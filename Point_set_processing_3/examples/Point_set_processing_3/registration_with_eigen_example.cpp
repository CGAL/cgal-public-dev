#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/draw_point_set_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/OpenGR/compute_registration_transformation.h>
#include <CGAL/OpenGR/register_point_sets.h>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <Eigen/SPQRSupport>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Vector_3                                     Vector;
typedef CGAL::Point_set_3<Point>                             PointSet;
typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef CGAL::Aff_transformation_3<Kernel>                   Transform;

typedef boost::associative_property_map< std::map<size_t, Point> >                          IndexMap;

typedef CGAL::Search_traits_3<Kernel>                                                       BaseTraits;
typedef CGAL::Search_traits_adapter<size_t, IndexMap, BaseTraits>                           Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>                                          Neighbor_search;
typedef Neighbor_search::Tree                                                               Tree;
typedef Tree::Splitter                                                                      Splitter;
typedef Neighbor_search::Distance                                                           Distance;

typedef Eigen::Triplet<double>                                                              T;

typedef std::unordered_map<size_t, size_t>                                                  Correspondence;

Mesh merge_meshes(Mesh a, Mesh b, const Correspondence& correspondence = {}) {
    // lots of redundant information bc of add_face
    // resizing the mesh before adding vertices and faces might be faster
    for (auto f : a.faces()) {
        std::array<Mesh::Vertex_index, 3> triangle;
        int i = 0;
        for (Mesh::Vertex_index v : a.vertices_around_face(a.halfedge(f))) {
            triangle[i] = b.add_vertex(a.point(v));
            if (correspondence.find(v) != correspondence.end()) {
                b.add_edge(v, triangle[i]);
			}
            ++i;
        }
        b.add_face(triangle[0], triangle[1], triangle[2]);
    }
    return b;
}

Mesh readToscaMesh(const std::string& filename) {
    Mesh mesh;
    std::string line;
    std::ifstream vertfile(filename + ".vert");
    while (std::getline(vertfile, line)) {
        std::istringstream iss(line);
        double x, y, z;
        iss >> x >> y >> z;
        Point p(x, y, z);
        Mesh::Vertex_index u = mesh.add_vertex(p);
    }
    vertfile.close();
    std::ifstream trifile(filename + ".tri");
    while (std::getline(trifile, line)) {
        std::istringstream iss(line);
        int u, v, w;
        iss >> u >> v >> w;
        mesh.add_face(Mesh::Vertex_index(u-1), Mesh::Vertex_index(v-1), Mesh::Vertex_index(w-1));
    }
    trifile.close();
    typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
    auto vnormals = mesh.add_property_map<vertex_descriptor, Vector>("v:normal", CGAL::NULL_VECTOR).first;
    CGAL::Polygon_mesh_processing::compute_vertex_normals(mesh, vnormals);
    return mesh;
}

PointSet Mesh2PointSet(const Mesh& mesh) {
    PointSet point_set;
    typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
    auto vertexnormals = mesh.property_map<vertex_descriptor, Vector>("v:normal").first;
    for (auto v : mesh.vertices()) {
        auto i = point_set.insert(mesh.point(v));
        point_set.normal(*i) = vertexnormals[v];
    }
    return point_set;
}

Mesh ModifyMeshWithPointSet(const PointSet& point_set, Mesh mesh) {
    for (auto v : mesh.vertices()) {
        mesh.point(v) = point_set.point(size_t(v));
    }
    return mesh;
}

Correspondence readCorrespondence(const std::string& filename) {
    Correspondence result;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        size_t p, q;
        iss >> p >> q;
        result.insert({ p, q });
    }
    file.close();
    return result;
}

// w1: point-to-point weight
// w2: point-to-plane weight
// w3: rigid transformation weight
std::pair<Transform, PointSet> rigid_registration(const PointSet& source, const PointSet& target,
    double w1 = 0.1, double w2 = 0.1, double w3 = 1.0,
    int max_iterations = 100, int num_threads = 1, const Correspondence& correspondence = {}) {
    Eigen::setNbThreads(num_threads);
    // resulting transformation is (R, t)
    Eigen::Matrix3d R_all = Eigen::Matrix3d::Identity();
    Eigen::Vector3d t_all = Eigen::Vector3d::Zero();
    size_t N = source.size();
    PointSet z = source;
    // build k-d tree for nearest neighbor search
    std::map<size_t, Point> index_map_base;
    for (auto it = target.begin(); it != target.end(); ++it) {
        index_map_base[*it] = target.point(*it);
    }
    IndexMap index_map(index_map_base);
    Distance distance(index_map);
    Tree tree(boost::counting_iterator<size_t>(0),
        boost::counting_iterator<size_t>(N),
        Splitter(),
        Traits(index_map)); // when called, returns the index of the nearest neighbor on target mesh
    // solver
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::Lower | Eigen::Upper> cg;
    cg.setMaxIterations(1000);
    cg.setTolerance(1e-6);
    // sparse coefficient matrix A
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(6 + 3 * N, 6 + 3 * N);
    std::vector<T> a(51 * N);
    size_t rolling = 0;
    // result vector b
    Eigen::VectorXd b(6 + 3 * N);
    double error = std::numeric_limits<double>::max();
    for (size_t iter = 0; iter < max_iterations; ++iter) {
        b.setZero();
        rolling = 0;
        auto start = std::chrono::high_resolution_clock::now();
        for (auto it = source.begin(); it != source.end(); ++it) {
            size_t i = *it;
            Point xp = source.point(*it);
            Eigen::Vector3d x(xp.x(), xp.y(), xp.z());
            // apply transformation to original x
            Eigen::Vector3d x_t = R_all * x + t_all;
            // search the nearest neighbor PI on the target mesh, n is the normal at PI
            Point zp = z.point(*it);
            size_t index = correspondence.find(i) != correspondence.end() ? correspondence.at(i) : Neighbor_search(tree, zp, 1, 0, true, distance).begin()->first;
            Point pnn = index_map[index];
            Eigen::Vector3d PI(pnn.x(), pnn.y(), pnn.z());
            Vector nnn = target.normal(index);
            Eigen::Vector3d n(nnn.x(), nnn.y(), nnn.z());
            // build A and b
            Eigen::Matrix3d n_matrix = n * n.transpose();
            Eigen::Matrix3d z_block = (1 + w1 / w3) * Eigen::Matrix3d::Identity() + w2 / w3 * n_matrix;
            Eigen::Matrix3d X_t;
            X_t << 0, -x_t(2), x_t(1),
                x_t(2), 0, -x_t(0),
                -x_t(1), x_t(0), 0;
            Eigen::Matrix3d XX_t(X_t * X_t);
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    a[rolling++] = T(j, k, -XX_t(j, k)); // A_rr
                    if (j != k) {
                        a[rolling++] = T(j, 3 + k, X_t(j, k)); // A_tr
                        a[rolling++] = T(3 + j, k, -X_t(j, k)); // A_rt
                        a[rolling++] = T(j, 6 + 3 * i + k, -X_t(j, k)); // A_zir
                        a[rolling++] = T(6 + 3 * i + j, k, X_t(j, k)); // A_rzj
                    } else {
                        a[rolling++] = T(3 + j, 3 + k, 1); // A_tt
                        a[rolling++] = T(3 + j, 6 + 3 * i + k, -1); // A_zit
                        a[rolling++] = T(6 + 3 * i + j, 3 + k, -1); // A_tzj
                    }
                    a[rolling++] = T(6 + 3 * i + j, 6 + 3 * i + k, z_block(j, k)); // A_zizj
                }
            }
            // b_r is zero
            b.segment<3>(3) -= x_t; // b_t
            b.segment<3>(6 + 3 * i) = (w1 / w3 * Eigen::Matrix3d::Identity() + w2 / w3 * n_matrix) * PI + x_t; // b_zi
        }
        A.setFromTriplets(a.begin(), a.end());
        if (iter == 0) {
            cg.analyzePattern(A);
        }
        cg.factorize(A);
        Eigen::VectorXd x = cg.solve(b);
        std::cout << "CG converged within " << cg.iterations() << " iterations. Error: " << cg.error() << "." << std::endl;
        auto stop = std::chrono::high_resolution_clock::now();
        std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << std::endl;
        // update r
        Eigen::Vector3d r(x[0], x[1], x[2]);
        auto roll = Eigen::AngleAxisd(x[0], Eigen::Vector3d::UnitX());
        auto pitch = Eigen::AngleAxisd(x[1], Eigen::Vector3d::UnitY());
        auto yaw = Eigen::AngleAxisd(x[2], Eigen::Vector3d::UnitZ());
        Eigen::Matrix3d R = (roll * pitch * yaw).matrix();
        R_all *= R;
        // update t
        Eigen::Vector3d t(x[3], x[4], x[5]);
        t_all += t;
        // update z
        for (auto it = z.begin(); it != z.end(); ++it) {
            size_t i = *it;
            Point p(x[6 + 3*i], x[6+ 3*i + 1], x[6 + 3*i+2]);
            z.point(*it) = p;
        }
        //// calculate error
        //double new_error = 0;
        //for (auto it = source.begin(); it != source.end(); ++it) {
        //    Point xp = source.point(*it);
        //    Eigen::Vector3d x(xp.x(), xp.y(), xp.z());
        //    Point zp = z.point(*it);
        //    Eigen::Vector3d ref_z(zp.x(), zp.y(), zp.z());
        //    size_t index = Neighbor_search(tree, zp, 1, 0, true, distance).begin()->first;
        //    Point pnn = index_map[index];
        //    Eigen::Vector3d PI(pnn.x(), pnn.y(), pnn.z());
        //    Vector nnn = target.normal(index);
        //    Eigen::Vector3d n(nnn.x(), nnn.y(), nnn.z());
        //    // point-to-point error
        //    new_error += w1 * (PI - ref_z).squaredNorm();
        //    // point-to-plane error
        //    new_error += w2 * std::pow(n.dot(PI - ref_z), 2);
        //    // rigid transformation error
        //    new_error += w3 * (R_all * x + t_all - ref_z).squaredNorm();
        //}
        //std::cout << "Iteration: " << iter << " Error: " << new_error << " \n" << R_all << "\n" << t_all << std::endl;
        //if (new_error < error) {
        //    error = new_error;
        //} else {
        //    break;
        //}
    }
    Transform transform(R_all(0, 0), R_all(0, 1), R_all(0, 2), t_all(0),
        R_all(1, 0), R_all(1, 1), R_all(1, 2), t_all(1),
        R_all(2, 0), R_all(2, 1), R_all(2, 2), t_all(2), 1);
    return std::make_pair(transform, z);
}







// This function implements various energy terms from paper (1) and (2) for non-rigid registration:
// - paper (1) https://taiya.github.io/pubs/tagliasacchi2016registration.pdf
// - paper (2) https://igl.ethz.ch/projects/ARAP/arap_web.pdf
// There are 5 energy terms in total:
// - point-to-point matching energy (equation 21 in paper (1))
// - point-to-plane matching energy (equation 22 in paper (1))
// - global rigidity energy (equation 15 in paper (1))
// - local rigidity energy (equation 16 in paper (1))
// - complementary as-rigid-as-possible energy (paper (2))
// The energy function is the sum of the above 5 energy terms.
// The energy function is minimized by solving a system of linear equations: Ax = b,
// where A is a sparse matrix and b is a vector.
// 
// Function parameters:
// - source: The surface we want to register. Since the function depends on connectivity information, it has to be a mesh.
// If you have a point set, see CGAL's surface reconstruction utilities to convert it to a mesh.
// - target: The surface we want to register the source to. It can be either a point set or a mesh.
// - w1: Weight of point-to-point matching energy.
// - w2: Weight of point-to-plane matching energy.
// - w3: Weight of global rigidity energy.
// - w4: Weight of local rigidity energy.
// - w5: Weight of complementary as-rigid-as-possible energy.
// - max_iter: Maximum number of iterations.
// - num_threads: Number of threads used for the solver (Eigen's LeastSquaresConjugateGradient).
// - correspondence: A map from source vertex indices to target vertex indices. Empty by default.
std::pair<Transform, Mesh> nonrigid_registration(const Mesh& source, const PointSet& target,
    double w1 = 1.0, double w2 = 1.0, double w3 = 0.1, double w4 = 10.0, double w5 = 1.0,
    size_t max_iter = 100, int num_threads = 1, const Correspondence& correspondence = {}) {
    
    Eigen::setNbThreads(num_threads);

    size_t N = source.number_of_vertices();
    size_t E = source.number_of_edges();

    // this will be the solution
    Mesh solution_mesh = source;
    // this will be the source mesh (and only transformed rigidly)
    Mesh source_mesh = source;

    // build k-d tree for nearest neighbor search
    std::map<size_t, Point> index_map_base;
    for (auto it = target.begin(); it != target.end(); ++it) {
        index_map_base[*it] = target.point(*it);
    }
    IndexMap index_map(index_map_base);
    Distance distance(index_map);
     // when called, returns the index of the nearest neighbor on target mesh
    Tree tree(boost::counting_iterator<size_t>(0),
              boost::counting_iterator<size_t>(N),
              Splitter(),
              Traits(index_map));

    // solver
    Eigen::LeastSquaresConjugateGradient< Eigen::SparseMatrix<double> > lscg; // other option: SPQR

    // for singular value decomposition in ARAP (see paper 2)
    Eigen::JacobiSVD<Eigen::Matrix3d> svd;

    // sparse coefficient matrix A
    Eigen::SparseMatrix<double> A(
        3 * N + // point-to-point matching energy
        N + // point-to-plane matching energy
        3 * N + // global rigidity energy
        3 * 2 * E + // local rigidity energy
        3 * N, // complementary as-rigid-as-possible energy
        3 * N + // updated positions
        6 + // global rotation and translation
        3 * N // local rotations
    );
    // we use an std::vector of triplets to construct A
    std::vector<T> a(21 * N + 15 * 2 * E);
    size_t rolling = 0;

    // preallocate constant elements of A
    size_t e = 0; // edge index
    for (auto it = source_mesh.vertices().begin(); it != source_mesh.vertices().end(); ++it) {
        size_t i = it->idx(); // source vertex index

        // if the corresponding point on target is given, increase the weight of the point
        bool is_corr = correspondence.find(i) != correspondence.end();
        double w1_corr = is_corr ? w1 * 1000 : w1;
        double w2_corr = is_corr ? w2 * 1000 : w2;

        // point-to-point matching energy - updated positions
        a[rolling++] = T(3 * i + 0, 3 * i + 0, -1 * w1_corr);
        a[rolling++] = T(3 * i + 1, 3 * i + 1, -1 * w1_corr);
        a[rolling++] = T(3 * i + 2, 3 * i + 2, -1 * w1_corr);

        // global rigidity energy - global translation
        a[rolling++] = T(3 * N + N + 3 * i + 0, 3 * N + 3, 1 * w3);
        a[rolling++] = T(3 * N + N + 3 * i + 1, 3 * N + 4, 1 * w3);
        a[rolling++] = T(3 * N + N + 3 * i + 2, 3 * N + 5, 1 * w3);

        // global rigidity energy - updated positions
        a[rolling++] = T(3 * N + N + 3 * i + 0, 3 * i + 0, -1 * w3);
        a[rolling++] = T(3 * N + N + 3 * i + 1, 3 * i + 1, -1 * w3);
        a[rolling++] = T(3 * N + N + 3 * i + 2, 3 * i + 2, -1 * w3);
        for (auto he : CGAL::halfedges_around_source(*it, source_mesh)) { // for each neighbor around source vertex
            auto v0 = CGAL::target(he, source_mesh);
            size_t j = v0.idx(); // target vertex index

            // local rigidity energy - updated positions
            a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 0, 3 * i + 0, -1 * w4);
            a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 1, 3 * i + 1, -1 * w4);
            a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 2, 3 * i + 2, -1 * w4);
            a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 0, 3 * j + 0, 1 * w4);
            a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 1, 3 * j + 1, 1 * w4);
            a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 2, 3 * j + 2, 1 * w4);

            // complementary as-rigid-as-possible energy (off-diagonal elements of graph Laplacian)
            a[rolling++] = T(3 * N + N + 3 * N + 3 * 2 * E + 3 * i + 0, 3 * j + 0, -1 * w5);
            a[rolling++] = T(3 * N + N + 3 * N + 3 * 2 * E + 3 * i + 1, 3 * j + 1, -1 * w5);
            a[rolling++] = T(3 * N + N + 3 * N + 3 * 2 * E + 3 * i + 2, 3 * j + 2, -1 * w5);

            ++e;
        }

        // complementary as-rigid-as-possible energy (diagonal elements of graph Laplacian)
        size_t diag = CGAL::halfedges_around_source(*it, source_mesh).size();
        a[rolling++] = T(3 * N + N + 3 * N + 3 * 2 * E + 3 * i + 0, 3 * i + 0, diag * w5);
        a[rolling++] = T(3 * N + N + 3 * N + 3 * 2 * E + 3 * i + 1, 3 * i + 1, diag * w5);
        a[rolling++] = T(3 * N + N + 3 * N + 3 * 2 * E + 3 * i + 2, 3 * i + 2, diag * w5);
    }

    // vector b
    Eigen::VectorXd b(
        3 * N + // point-to-point matching energy
        N + // point-to-plane matching energy
        3 * N + // global rigidity energy
        3 * 2 * E + // local rigidity energy
        3 * N // complementary as-rigid-as-possible energy
    );

    // iterative minimization
    for (size_t iter = 0; iter < max_iter; ++iter) {

        auto start = std::chrono::high_resolution_clock::now();

        double error = 0.0;

        // we collect the arap rotation matrices in an std::vector (see paper 2)
        std::vector<Eigen::Matrix3d> Rotations(N);

        b.setZero();
        // the first 12 * N + 9 * 2 * E triplets of A have constant values
        rolling = 12 * N + 9 * 2 * E;

        size_t e = 0; // edge index
        for (auto it = source_mesh.vertices().begin(); it != source_mesh.vertices().end(); ++it) {
            size_t i = it->idx(); // source vertex index

            // if the corresponding point on target is given, increase the weight of the point
            bool is_corr = correspondence.find(i) != correspondence.end();
            double w1_corr = is_corr ? w1 * 1000 : w1;
            double w2_corr = is_corr ? w2 * 1000 : w2;

            // solution_p is used to find better nearest neighbors on the target gradually in each iteration
            Point solution_p = solution_mesh.point(*it);
            // either use a given correspondence or find the index of nearest neighbor on target
            size_t index = is_corr ? correspondence.at(i) : Neighbor_search(tree, solution_p, 1, 0, true, distance).begin()->first;
            // position of nearest neighbor on target
            Point pnn = index_map[index];
            Eigen::Vector3d PI(pnn.x(), pnn.y(), pnn.z());
            error += CGAL::sqrt(CGAL::squared_distance(solution_p, pnn));
            // normal of nearest neighbor on target
            const Vector nnn = target.normal(index);
            const Eigen::Vector3d n(nnn.x(), nnn.y(), nnn.z());

            // position of point
            Point p = source_mesh.point(*it);
            Eigen::Vector3d x(p.x(), p.y(), p.z());

            // point-to-point matching energy - constants
            b[3 * i + 0] = -PI[0] * w1_corr;
            b[3 * i + 1] = -PI[1] * w1_corr;
            b[3 * i + 2] = -PI[2] * w1_corr;

            // point-to-plane matching energy - updated positions
            a[rolling++] = T(3 * N + i, 3 * i + 0, -n[0] * w2_corr);
            a[rolling++] = T(3 * N + i, 3 * i + 1, -n[1] * w2_corr);
            a[rolling++] = T(3 * N + i, 3 * i + 2, -n[2] * w2_corr);

            // point-to-plane matching energy - constants
            b[3 * N + i] = -n.dot(PI) * w2_corr;

            // global rigidity energy - global rotation
            a[rolling++] = T(3 * N + N + 3 * i + 0, 3 * N + 1, x[2] * w3);
            a[rolling++] = T(3 * N + N + 3 * i + 0, 3 * N + 2, -x[1] * w3);
            a[rolling++] = T(3 * N + N + 3 * i + 1, 3 * N + 0, -x[2] * w3);
            a[rolling++] = T(3 * N + N + 3 * i + 1, 3 * N + 2, x[0] * w3);
            a[rolling++] = T(3 * N + N + 3 * i + 2, 3 * N + 0, x[1] * w3);
            a[rolling++] = T(3 * N + N + 3 * i + 2, 3 * N + 1, -x[0] * w3);

            // local rigidity energy - constants
            b[3 * N + N + 3 * i + 0] = -x[0] * w3;
            b[3 * N + N + 3 * i + 1] = -x[1] * w3;
            b[3 * N + N + 3 * i + 2] = -x[2] * w3;

            // original position of point (for arap)
            Point p0 = source.point(*it);
            Eigen::Vector3d x0(p0.x(), p0.y(), p0.z());

            // covariance matrix for arap (see paper 2 equation 5)
            Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();

            for (auto he : CGAL::halfedges_around_source(*it, source_mesh)) { // for each neighbor around source vertex
                auto v0 = CGAL::target(he, source_mesh);
                size_t j = v0.idx(); // target vertex index

                // position of neighbor
                Point pk = source_mesh.point(v0);
                Eigen::Vector3d x_k(pk.x(), pk.y(), pk.z());

                // local rigidity energy - local rotations
                a[rolling++] = T(3 * N + N +  3 * N + 3 * e + 0, 3 * N + 6 + 3 * i + 1, (x[2] - x_k[2]) * w4);
                a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 0, 3 * N + 6 + 3 * i + 2, (x_k[1] - x[1]) * w4);
                a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 1, 3 * N + 6 + 3 * i + 0, (x_k[2] - x[2]) * w4);
                a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 1, 3 * N + 6 + 3 * i + 2, (x[0] - x_k[0]) * w4);
                a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 2, 3 * N + 6 + 3 * i + 0, (x[1] - x_k[1]) * w4);
                a[rolling++] = T(3 * N + N + 3 * N + 3 * e + 2, 3 * N + 6 + 3 * i + 1, (x_k[0] - x[0]) * w4);

                // local rigidity energy - constants
                b[3 * N + N + 3 * N + 3 * e + 1] = (x_k[1] - x[1]) * w4;
                b[3 * N + N + 3 * N + 3 * e + 2] = (x_k[2] - x[2]) * w4;
                b[3 * N + N + 3 * N + 3 * e + 0] = (x_k[0] - x[0]) * w4;

                ++e;

                // original position of neighbor (for arap)
                Point p0k = source.point(v0);
                Eigen::Vector3d x0_k(p0k.x(), p0k.y(), p0k.z());

                // for arap: (p_i - p_j) * (p_i' - p_j')^T (see paper 2 equation 5)
                covariance += (x0 - x0_k) * (x - x_k).transpose();
            }

            // rotation matrices for arap (see paper 2 equation 6)
            svd.compute(covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Rotations[i] = svd.matrixV() * svd.matrixU().transpose();
            // correct if it contains reflection
            double det = Rotations[i].determinant();
            if (det < 0) {
                Eigen::Matrix3d M = Eigen::Matrix3d::Identity();
            	M(2, 2) = -1;
                Rotations[i] = svd.matrixV() * M * svd.matrixU().transpose();
            }

        }

        // complementary as-rigid-as-possible energy - constants (see paper 2 equation 8)
        for (auto it = source_mesh.vertices().begin(); it != source_mesh.vertices().end(); ++it) {
            size_t i = it->idx();
            Point p0 = source.point(*it);
            Eigen::Vector3d x0(p0.x(), p0.y(), p0.z());
            for (auto he : CGAL::halfedges_around_source(*it, source_mesh)) {
                auto v0 = CGAL::target(he, source_mesh);
                size_t j = v0.idx();
                Point p0k = source.point(v0);
                Eigen::Vector3d x0_k(p0k.x(), p0k.y(), p0k.z());
                b.segment<3>(3 * N + N + 3 * N + 3 * 2 * E + 3 * i) = 0.5 * (Rotations[i] + Rotations[j]) * (x0 - x0_k) * w5;
            }
        }

        A.setFromTriplets(a.begin(), a.end()); // automatically calls makeCompressed()

        // solve Ax = b
        lscg.compute(A);
        Eigen::VectorXd solution = lscg.solve(b);
        // global rotation
        Eigen::Vector3d r(solution[3 * N + 0], solution[3 * N + 1], solution[3 * N + 2]);
        auto roll = Eigen::AngleAxisd(r[0], Eigen::Vector3d::UnitX());
        auto pitch = Eigen::AngleAxisd(r[1], Eigen::Vector3d::UnitY());
        auto yaw = Eigen::AngleAxisd(r[2], Eigen::Vector3d::UnitZ());
        Eigen::Matrix3d R = (roll * pitch * yaw).matrix();
        // global translation
        Eigen::Vector3d t(solution[3 * N + 3], solution[3 * N + 4], solution[3 * N + 5]);
        //std::cout << "R: " << std::endl << R << std::endl;
        //std::cout << "t: " << std::endl << t << std::endl;
        // updated positions
        
        for (auto it = source_mesh.vertices().begin(); it != source_mesh.vertices().end(); ++it) {
            size_t i = it->idx();
            solution_mesh.point(*it) = Point(solution[3 * i + 0], solution[3 * i + 1], solution[3 * i + 2]);
            // rigid transformation of source mesh
            Point p = source_mesh.point(*it);
            Eigen::Vector3d x(p.x(), p.y(), p.z());
            Eigen::Vector3d x_transformed = R * x + t;
            source_mesh.point(*it) = Point(x_transformed[0], x_transformed[1], x_transformed[2]);
        }

        auto stop = std::chrono::high_resolution_clock::now();
        
        std::cout << "Iteration " << iter << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms.";
		std::cout << " Chamfer distance from target: " << error << '.' << std::endl;
    }

    return std::make_pair(Transform(), solution_mesh);

}


int main(int argc, char* argv[]) {
    //const std::string filename1 = CGAL::data_file_path("../meshes/wolf1.ply");
    //Mesh mesh1;
    //CGAL::IO::read_PLY(std::ifstream(filename1), mesh1);
    //const std::string filename2 = CGAL::data_file_path("../meshes/wolf2.ply");
    //Mesh mesh2;
    //CGAL::IO::read_PLY(std::ifstream(filename2), mesh2);
    Mesh mesh1 = readToscaMesh("../meshes/toscahires-asci/wolf0");
    Mesh mesh2 = readToscaMesh("../meshes/toscahires-asci/wolf1");
    Correspondence correspondence = readCorrespondence("../meshes/correspondence.txt");
    CGAL::draw(merge_meshes(mesh1, mesh2, correspondence));
    
    //auto opengr_result = CGAL::OpenGR::compute_registration_transformation(Mesh2PointSet(mesh2), Mesh2PointSet(mesh1));
    //Transform opengr_transform = opengr_result.first;
    //CGAL::Polygon_mesh_processing::transform(opengr_transform, mesh1);
    //CGAL::draw(merge_meshes(mesh1, mesh2));

    auto rigid_result = rigid_registration(Mesh2PointSet(mesh1), Mesh2PointSet(mesh2), 0.1, 0.1, 1.0, 10, 4);
    Transform transform = rigid_result.first;
    CGAL::Polygon_mesh_processing::transform(transform, mesh1);
    //CGAL::Polygon_mesh_processing::transform(Transform(CGAL::Translation(), Vector(-40, 0, 0)), mesh1);
    //CGAL::draw(merge_meshes(mesh1, mesh2));
    auto nonrigid_result = nonrigid_registration(mesh1, Mesh2PointSet(mesh2), 1.0, 1.0, 0.1, 10.0, 0.0, 10, 4);
    Mesh z = nonrigid_result.second;
    CGAL::Polygon_mesh_processing::transform(Transform(CGAL::Translation(), Vector(-40, 0, 0)), z);
    CGAL::draw(merge_meshes(z, mesh2, correspondence));
    return EXIT_SUCCESS;
}
