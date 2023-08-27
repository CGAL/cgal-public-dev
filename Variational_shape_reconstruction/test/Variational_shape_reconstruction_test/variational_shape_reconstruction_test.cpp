#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/centroid.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Variational_shape_reconstruction.h>
#include <iostream>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Point_set_3/IO/XYZ.h>
// CGAL
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;


typedef Kernel::FT                  FT;
typedef Kernel::Point_3             Point;
typedef Kernel::Vector_3            Vector;
typedef Kernel::Triangle_3          Triangle;
typedef Kernel::Plane_3             Plane;

typedef Kernel::Point_2             Point_2;
typedef Kernel::Segment_2           Segment_2;

typedef CGAL::First_of_pair_property_map<std::pair<Point, Vector>>                     Point_map;
typedef CGAL::Second_of_pair_property_map<std::pair<Point, Vector>>                    Normal_map;
typedef CGAL::Point_set_3< Point, Vector > Pointset;

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <fstream>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>
// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
FT sphere_function (Point_3 p) {
  const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2+y2+z2-1;
}
FT torus_function (Point_3 p) {
  const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  const auto f = (x2+y2+z2-3*3-2*2); 
  return f*f - 4*3*3*(2*2-z2);
}
Pointset resize(Pointset pointset)
{
    std::vector<Point> m_points2;
    for( Pointset::const_iterator pointset_it = pointset.begin(); pointset_it != pointset.end(); ++ pointset_it )
    {
        const auto point = pointset.point(*pointset_it);
        m_points2.push_back(point);
    }

    auto m_bbox = CGAL::bbox_3(m_points2.begin(), m_points2.end());
    double dx = m_bbox.xmax() - m_bbox.xmin();
    float scaling = 100./dx;
    Point centroid = CGAL::centroid(m_points2.begin(), m_points2.end(),CGAL::Dimension_tag<0>());
    Aff_transformation trans(CGAL::TRANSLATION, Vector(-centroid.x(),-centroid.y(),-centroid.z()));
    Aff_transformation scale(CGAL::SCALING, scaling);

    m_points2.clear();
    for( Pointset::iterator pointset_it = pointset.begin(); pointset_it != pointset.end(); ++ pointset_it )
    {
        auto point = pointset.point(*pointset_it);
        pointset.point(*pointset_it)=scale(trans(point));
    }
    
    return pointset;
}
class MeshAnalysis
{
    public:
        MeshAnalysis(Polyhedron mesh)
        {
            m_mesh = mesh;
        }
        size_t number_connected_component()
        {
            return m_mesh.keep_largest_connected_components(1)+1;
  
        }
        int genus()
        {
            int F = m_mesh.size_of_facets();
            int V = m_mesh.size_of_vertices();
            int E = m_mesh.size_of_halfedges()/2;


            return  1 - (V-E+F)/2;
  
        }
        void test()
        {
              Tr tr;            // 3D-Delaunay triangulation
                C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
                // defining the surface
                Surface_3 surface(torus_function,             // pointer to function
                                    Sphere_3(CGAL::ORIGIN, 3.)); // bounding sphere
                // Note that "2." above is the *squared* radius of the bounding sphere!
                // defining meshing criteria
                CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                                    0.1,  // radius bound
                                                                    0.1); // distance bound
                // meshing surface
                CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
                Surface_mesh sm;
                CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);
                
                int E = sm.edges().size();
                int F = sm.faces().size();
                int V = sm.vertices().size();

               
                std::cout<<"gen "<<V-E+F<<" \n";
            }
    private:
     Polyhedron  m_mesh;
};

void test_point_set(const std::string fname)
{
    Pointset pointset;
    if (!CGAL::IO::read_XYZ( "../data/"+fname+".xyz",pointset))
    {
        std::cerr << "Error: cannot read file " << std::endl;
        return;
    } 
    const size_t generators = 30;
    const size_t steps = 10;
    const double split_threshold =0.01;
    const double distance_weight =0.000001;
    // reconstruction
    const double  dist_ratio = 0.001;
	const double  fitting = 0.43;
	const double  coverage = 0.27;
	const double  complexity = 0.3;

    size_t new_generators = generators; 
    size_t iteration = 0 ;
    //116.34
    pointset = resize(pointset);
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	qem::Variational_shape_reconstruction manager(pointset,generators,distance_weight,3,3);
        std::chrono::steady_clock::time_point begin_clustering = std::chrono::steady_clock::now();
    while(new_generators > 5 )
    {
        manager.region_growing(steps);
        new_generators = manager.guided_split_clusters(split_threshold, iteration++);
    }
    std::chrono::steady_clock::time_point end_clustering = std::chrono::steady_clock::now();
    std::cerr << "Clustering " << std::chrono::duration_cast<std::chrono::microseconds>(end_clustering - begin_clustering).count()/1000 << "[ms]" << std::endl;
    
    // Reconstruction
    manager.reconstruction(dist_ratio, fitting, coverage, complexity);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cerr << "Algo " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000 << "[ms]" << std::endl;


    Pointset point_cloud = manager.get_point_cloud_clustered();
    Polyhedron mesh = manager.get_reconstructed();

    //CGAL::IO::write_XYZ("clustering_"+fname+".txt", point_cloud);

    std::ofstream edge_file;
    edge_file.open("output/clustering_"+fname+".ply");

    edge_file << "ply\n"
                << "format ascii 1.0\n"
                << "element vertex " << point_cloud.size() << "\n"
                << "property float x\n"
                << "property float y\n"
                << "property float z\n"
                << "property uchar red\n"
                << "property uchar green\n"
                << "property uchar blue\n"
                << "end_header\n";
    std::vector<Point> points;
    for(Pointset::const_iterator it = point_cloud.begin(); it != point_cloud.end(); ++ it)
    {
        auto point = point_cloud.point(*it);
        points.push_back(point);
        edge_file << point.x() << " " << point.y() << " " << point.z() << " ";
        auto normal = point_cloud.normal(*it);
        edge_file << static_cast<int>(255*normal.x()) << " " << static_cast<int>(255*normal.y()) << " " << static_cast<int>(255*normal.z()) << std::endl;
    }

    edge_file.close();

    std::ofstream metrics_file;
    metrics_file.open("output/metrics_"+fname+".txt");
    if(!mesh.is_valid())
    {
        metrics_file<<"ERROR\n";

    }
    double dist = CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set	(mesh,points,0.0001);	
    MeshAnalysis mesh_analysis(mesh);
    metrics_file<<"genus "<<mesh_analysis.genus()<<"\n";
    metrics_file<<"number_connected_component "<<mesh_analysis.number_connected_component()<<"\n";
    metrics_file<<"is_valid "<<mesh.is_valid()<<"\n";
    metrics_file<<"is_pure_triangle "<<mesh.is_pure_triangle()<<"\n";
    metrics_file<<"is_closed "<<mesh.is_closed()<<"\n";
    metrics_file<<"distance "<<dist<<"\n";
    metrics_file.close();


    
    std::ofstream mesh_file;
    mesh_file.open("output/mesh_"+fname+".off");
    CGAL::write_off(mesh_file, mesh);
    mesh_file.close();
    manager.write_csv();
    
}

int main(int argc, char** argv)
{	
    const std::vector<std::string> files_to_test{
        //"double_sphere",
        "cubes",
        //"bones", //ok
        //"spheres", //ok
        //"piece_meca", //ok
        //"armjoin", // cluster ok not recons
        "guitar", // cluster ok not recons
        "g", // cluster ok not recons
        "capsule", //ok
        "hilbert_cube2_pds_100k",// cluster ok not recons
        "bunny_150k",// cluster ok not recons
        "joint",// ok
        "fertility",// cluster ok not recons
        "qtr_piston",// cluster ok recons weird
        "hand1" // ok
    };   
    int n = 1;
    
    for(int i = 0 ; i < n ; i ++)
    {
        std::cout<<"mesh "+files_to_test[i] +" \n";
        test_point_set(files_to_test[i]);
    }
	return 0;
}


