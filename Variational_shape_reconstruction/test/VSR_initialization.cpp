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
void test_initialization()
{
    std::string fname = "spheres";
    Pointset pointset;
    if (!CGAL::IO::read_XYZ( "../data/"+fname+".xyz",pointset))
    {
        std::cerr << "Error: cannot read file " << std::endl;
        return;
    } 
    const size_t generators = 4;
    const size_t steps = 10;
    const double split_threshold = 10e-2;

    // reconstruction
    const double  dist_ratio = 10e-3;
	const double  fitting = 0.4;
	const double  coverage = 0.3;
	const double  complexity = 0.3;

    size_t new_generators = generators; 
    size_t iteration = 0 ;
    const double distance_weight = 10e-5;
    
    // for the test
    const int min_generator_count_per_component = 1; 

    
    qem::Variational_shape_reconstruction manager(
        pointset,
        generators,
        distance_weight,
        qem::VERBOSE_LEVEL::LOW,
        qem::INIT_QEM_GENERATOR::RANDOM);
    std::cout<< "TEST INITIALIZATION \n";
    auto generators_per_component = manager.get_generators_per_component();
    bool test_status = true;
    if(generators_per_component.size() < manager.get_component_count())
        test_status = false;
        
    for(auto generator_count : generators_per_component)
    {
        if( generator_count < min_generator_count_per_component)
            test_status = false; 
    }

    if(test_status)
        std::cout<< "test passed ! \n";
    else 
        std::cout<< "test failed ! \n"; 
        


}
int main(int argc, char** argv)
{	
    test_initialization();
	return 0;
}
