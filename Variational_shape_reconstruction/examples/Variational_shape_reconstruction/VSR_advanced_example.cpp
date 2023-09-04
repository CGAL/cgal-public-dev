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

int main()
{	

    Pointset pointset;
    if (!CGAL::IO::read_XYZ( "../data/guitar.xyz",pointset))
    {
        std::cerr << "Error: cannot read file " << std::endl;
        return EXIT_FAILURE;
    } 
    size_t generators = 30; 
    const size_t steps = 10;
    const double split_threshold = 10e-2;
    const double distance_weight = 10e-5;
    size_t iteration = 0 ;
	
    qem::Variational_shape_reconstruction manager(
        pointset,
        generators,
        distance_weight,
        qem::VERBOSE_LEVEL::HIGH,
        qem::INIT_QEM_GENERATOR::KMEANS_PLUSPLUS);
    while(generators > 5 )
    {
        
		manager.region_growing_and_update_poles(steps);
        generators = manager.guided_split_clusters(split_threshold, iteration++);
    }
    
    // Reconstruction parameters
    const double dist_ratio = 10e-3;
	const double fitting = 0.4;
	const double coverage = 0.3;
	const double complexity = 0.3;
	
    
    manager.reconstruction(dist_ratio, fitting, coverage, complexity,false);

	auto mesh = manager.get_reconstructed_mesh();
    std::ofstream mesh_file;
    mesh_file.open("mesh.off");
    CGAL::write_off(mesh_file, mesh);
    mesh_file.close();

	return EXIT_SUCCESS;
}


