#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Variational_shape_reconstruction.h>
#include <iostream>

#include <CGAL/Point_set_3/IO/XYZ.h>
// CGAL
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Point_set_3.h>

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

    // reconstruction
    const double  dist_ratio = 0.001;
	const double  fitting = 0.43;
	const double  coverage = 0.27;
	const double  complexity = 0.3;

    size_t new_generators = generators; 
    size_t iteration = 0 ;

    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	qem::Variational_shape_reconstruction manager(pointset,generators);
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
    edge_file.open("clustering_"+fname+".ply");

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

    for(Pointset::const_iterator it = point_cloud.begin(); it != point_cloud.end(); ++ it)
    {
        auto point = point_cloud.point(*it);
        edge_file << point.x() << " " << point.y() << " " << point.z() << " ";
        auto normal = point_cloud.normal(*it);
        edge_file << static_cast<int>(255*normal.x()) << " " << static_cast<int>(255*normal.y()) << " " << static_cast<int>(255*normal.z()) << std::endl;
    }

    edge_file.close();

    if(!mesh.is_valid())
    {
        std::cout<<"ERROR";
    }
    std::ofstream mesh_file;
    mesh_file.open("mesh_"+fname+".off");
    CGAL::write_off(mesh_file, mesh);
    mesh_file.close();
    
}
int main(int argc, char** argv)
{	
    const std::vector<std::string> files_to_test{
        //"armjoin", // cluster ok not recons
        //"guitar", // cluster ok not recons
        "piece_meca", //ok
        //"g", // cluster ok not recons
        //"capsule", //ok
        //"hilbert_cube2_pds_100k",// cluster ok not recons
        //"bunny_150k",// cluster ok not recons
        //"joint",// ok
        //"fertility",// cluster ok not recons
        //"qtr_piston",// cluster ok recons weird
        "hand1" // ok

    };   
    int n = 1;
    for(int i = 0 ; i < n ; i ++)
    {
        test_point_set(files_to_test[i]);
    }
	return 0;
}


