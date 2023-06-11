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

int main(int argc, char** argv)
{	

    const std::string fname = argc > 1 ? argv[1] :  CGAL::data_file_path("data/hand1.xyz");
    Pointset pointset;
    if (!CGAL::IO::read_XYZ( fname,pointset))
    {
        std::cerr << "Error: cannot read file " << std::endl;
        return 0;
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
	return 0;
}