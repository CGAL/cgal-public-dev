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
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>
#include <Eigen/KLUSupport>


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


Mesh merge_meshes(Mesh a, Mesh b) {
    for (auto f : a.faces()) {
        std::array<Mesh::Vertex_index, 3> triangle;
        int i = 0;
        for (Mesh::Vertex_index v : a.vertices_around_face(a.halfedge(f))) {
            triangle[i] = b.add_vertex(a.point(v));
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
        point_set.insert(mesh.point(v));
        point_set.add_normal_map(vertexnormals[v]);
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


// w1: point-to-point
// w2: point-to-plane
// w3: rigid transformation
// w4: as-rigid-as-possible transformation
std::pair<Transform, Mesh> nonrigid_registration(const Mesh& source, const PointSet& target,
    double w1 = 0.1, double w2 = 0.1, double w3 = 1.0, double w4 = 10.0,
    int max_iter = 100, int num_threads = 1, const Correspondence& correspondence = {}) {
    double avg_filling = 0.0, avg_factorizing = 0.0, avg_solving = 0.0;
    Eigen::setNbThreads(num_threads);
    size_t N = source.number_of_vertices();
    size_t E = source.number_of_edges();
    Mesh z = source;
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
    Eigen::BiCGSTAB< Eigen::SparseMatrix<double, Eigen::RowMajor> > cg;
    cg.setMaxIterations(1000);
    cg.setTolerance(1e-5);
    // sparse coefficient matrix A
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(6 + 2 * 3 * N, 6 + 2 * 3 * N);
    std::vector<T> a(72 * N + 30 * E);
    size_t rolling = 0;
    // result vector b
    Eigen::VectorXd b(6 + 2 * 3 * N);
    double error = std::numeric_limits<double>::max();
    //std::vector<double> weights;
    //   for (size_t i = 0; i < N; ++i) {
    //       if (rand() % 1 == 0) {
    //           weights.push_back(1.0);
    //       }
    //}
    for (size_t iter = 0; iter < max_iter; ++iter) {
        b.setZero();
        rolling = 0;
        auto start = std::chrono::high_resolution_clock::now();
        for (auto it = z.vertices().begin(); it != z.vertices().end(); ++it) {
            auto i = it->idx();
            bool is_corr = correspondence.find(i) != correspondence.end();
            double w1_corr = is_corr ? w1 * 1000 : w1;
            double w2_corr = is_corr ? w2 * 1000 : w2;
            Point zp = z.point(*it);
            Eigen::Vector3d x_t(zp.x(), zp.y(), zp.z());
            // search the nearest neighbor PI on the target mesh, n is the normal at PI
            size_t index = is_corr ? correspondence.at(i) : Neighbor_search(tree, zp, 1, 0, true, distance).begin()->first;
            Point pnn = index_map[index];
            Eigen::Vector3d PI(pnn.x(), pnn.y(), pnn.z());
            Vector nnn = target.normal(index);
            Eigen::Vector3d n(nnn.x(), nnn.y(), nnn.z());
            // build A and b
            Eigen::Matrix3d n_matrix = n * n.transpose();
            Eigen::Matrix3d z_diag_block = (w1_corr + w3) * Eigen::Matrix3d::Identity() + w2_corr * n_matrix;
            Eigen::Matrix3d X_t;
            X_t << 0, -x_t(2), x_t(1),
                x_t(2), 0, -x_t(0),
                -x_t(1), x_t(0), 0;
            Eigen::Matrix3d XX_t(X_t * X_t);
            Eigen::Matrix3d rirj_diag_block = Eigen::Matrix3d::Zero();
            Eigen::Matrix3d zirj_diag_block = Eigen::Matrix3d::Zero();
            Eigen::Matrix3d rizj_diag_block = Eigen::Matrix3d::Zero();
            for (auto he : CGAL::halfedges_around_target(*it, z)) {
                auto v0 = CGAL::source(he, z);
                Point xpk = z.point(v0);
                Eigen::Vector3d x_t_k(xpk.x(), xpk.y(), xpk.z());
                Eigen::Matrix3d X_t_k;
                X_t_k << 0, -x_t_k(2), x_t_k(1),
                    x_t_k(2), 0, -x_t_k(0),
                    -x_t_k(1), x_t_k(0), 0;
                rirj_diag_block -= (X_t_k - X_t) * (X_t_k - X_t);
                zirj_diag_block += X_t_k - X_t;
                rizj_diag_block -= w4 * (X_t_k - X_t);
                z_diag_block += 2 * w4 * Eigen::Matrix3d::Identity();
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        if (j != k) {
                            a[rolling++] = T(6 + 3 * size_t(v0) + j, 6 + 3 * N + 3 * i + k, X_t_k(j, k) - X_t(j, k)); // zirj off diag
                            a[rolling++] = T(6 + 3 * N + 3 * size_t(v0) + j, 6 + 3 * i + k, w4 * (X_t_k(j, k) - X_t(j, k))); // rizj off diag
                        } else {
                            a[rolling++] = T(6 + 3 * N + 3 * size_t(v0) + j, 6 + 3 * N + 3 * i + k, -2 * w4); // zizj off diag
                        }
                    }
                }
            }
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    a[rolling++] = T(j, k, -XX_t(j, k)); // A_rr
                    if (j != k) {
                        a[rolling++] = T(j, 3 + k, X_t(j, k)); // A_tr
                        a[rolling++] = T(3 + j, k, -X_t(j, k)); // A_rt
                        a[rolling++] = T(j, 6 + 3 * N + 3 * i + k, -X_t(j, k)); // A_zir
                        a[rolling++] = T(6 + 3 * N + 3 * i + j, k, w3 * X_t(j, k)); // A_rzj
                        a[rolling++] = T(6 + 3 * i + j, 6 + 3 * N + 3 * i + k, zirj_diag_block(j, k)); // A_zirj
                        a[rolling++] = T(6 + 3 * N + 3 * i + j, 6 + 3 * i + k, rizj_diag_block(j, k)); // A_rizj
                    } else {
                        a[rolling++] = T(3 + j, 3 + k, 1); // A_tt
                        a[rolling++] = T(3 + j, 6 + 3 * N + 3 * i + k, -1); // A_zit
                        a[rolling++] = T(6 + 3 * N + 3 * i + j, 3 + k, -w3); // A_tzj
                    }
                    a[rolling++] = T(6 + 3 * i + j, 6 + 3 * i + k, rirj_diag_block(j, k)); // A_rirj
                    a[rolling++] = T(6 + 3 * N + 3 * i + j, 6 + 3 * N + 3 * i + k, z_diag_block(j, k)); // A_zizj
                }
            }
            // b_r is zero
            b.segment<3>(3) -= x_t; // b_t
            // b_ri is zero
            b.segment<3>(6 + 3 * N + 3 * i) = (w1_corr * Eigen::Matrix3d::Identity() + w2_corr * n_matrix) * PI + w3 * x_t; // b_zi
            for (auto he : CGAL::halfedges_around_target(*it, z)) {
                auto v0 = CGAL::source(he, z);
                Point xpk = z.point(v0);
                Eigen::Vector3d x_t_k(xpk.x(), xpk.y(), xpk.z());
                b.segment<3>(6 + 3 * N + 3 * i) -= 2 * w4 * (x_t_k - x_t);
            }
        }
        A.setFromTriplets(a.begin(), a.end());
        auto stop = std::chrono::high_resolution_clock::now();
        double filling = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        if (iter == 0) {
            cg.analyzePattern(A);
        }
        start = std::chrono::high_resolution_clock::now();
        cg.factorize(A);
        stop = std::chrono::high_resolution_clock::now();
        double factorizing = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        start = std::chrono::high_resolution_clock::now();
        Eigen::VectorXd solution = cg.solve(b);
        stop = std::chrono::high_resolution_clock::now();
        double solving = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        std::cout << "Time for filling A: " << filling << std::endl;
        std::cout << "Time for factorizing A: " << factorizing << std::endl;
        std::cout << "Time for solving Ax=b: " << solving << std::endl;
        if (iter != 0) {
            avg_filling += filling;
            avg_factorizing += factorizing;
            avg_solving += solving;
        }
        // std::cout << "CG converged within " << cg.iterations() << " iterations. Error: " << cg.error() << "." << std::endl;
        // update z
        for (auto it = z.vertices().begin(); it != z.vertices().end(); ++it) {
            auto i = it->idx();
            Point p(solution[6 + 3 * N + 3 * i], solution[6 + 3 * N + 3 * i + 1], solution[6 + 3 * N + 3 * i + 2]);
            z.point(*it) = p;
        }
        std::cout << "Iteration: " << iter << std::endl;
    }
    std::cout << "Average time for filling A: " << avg_filling / (max_iter - 1) << std::endl;
    std::cout << "Average time for factorizing A: " << avg_factorizing / (max_iter - 1) << std::endl;
    std::cout << "Average time for solving Ax=b: " << avg_solving / (max_iter - 1) << std::endl;
    std::cout << "Number of threads: " << Eigen::nbThreads( ) << std::endl;
    Transform transform(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1);
    return std::make_pair(transform, z);
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
    CGAL::draw(merge_meshes(mesh1, mesh2));
    
    //auto opengr_result = CGAL::OpenGR::compute_registration_transformation(Mesh2PointSet(mesh2), Mesh2PointSet(mesh1));
    //Transform opengr_transform = opengr_result.first;
    //CGAL::Polygon_mesh_processing::transform(opengr_transform, mesh1);
    //CGAL::draw(merge_meshes(mesh1, mesh2));
    
    auto rigid_result = rigid_registration(Mesh2PointSet(mesh1), Mesh2PointSet(mesh2), 0.1, 0.1, 1.0, 10, 4, correspondence);
    Transform transform = rigid_result.first;
    CGAL::Polygon_mesh_processing::transform(transform, mesh1);
    //CGAL::Polygon_mesh_processing::transform(Transform(CGAL::Translation(), Vector(-40, 0, 0)), mesh1);
    //CGAL::draw(merge_meshes(mesh1, mesh2));
    auto nonrigid_result = nonrigid_registration(mesh1, Mesh2PointSet(mesh2), 0.1, 0.1, 1.0, 10.0, 10, 4, correspondence);
    Mesh z = nonrigid_result.second;
    CGAL::Polygon_mesh_processing::transform(Transform(CGAL::Translation(), Vector(-40, 0, 0)), z);
    CGAL::draw(merge_meshes(z, mesh2));
    return EXIT_SUCCESS;
}



//import random
//
//def sample(ratio):
//  with open('correspondence.txt', 'w') as file:
//    for i in range(0, 4344):
//      if random.random() < ratio:
//        file.write(str(i) + ' ' + str(i) + '\n')
//
//sample(0.5)
