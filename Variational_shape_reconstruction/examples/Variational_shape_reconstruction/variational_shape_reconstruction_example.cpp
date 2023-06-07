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

int main(int argv, char **args)
{	
    std::ifstream fstream("../data/piece_meca.xyz");
    Pointset pointset;
    if (!fstream || !CGAL::IO::read_XYZ( fstream,pointset))
    {
        std::cerr << "Error: cannot read file " << std::endl;
        return 0;
    } 
    const size_t generators = 30;
    const size_t steps = 10;
    const double split_threshold =0.01;
    size_t new_generators = generators; 
    size_t iteration = 0 ;

	CGAL::Variational_shape_reconstruction manager;
	manager.initialize(pointset);
    
	manager.init_random_poles(generators);
    while(new_generators > 5 )
    {
        bool flag =1;
        // Clustering
        for(int i = 0 ; i < steps && flag;i++)
        {
            manager.region_growing(true);
            flag = manager.update_poles();
        }
        new_generators = manager.guided_split_clusters(split_threshold, iteration++);
    }

    
    // Reconstruction
	return 0;
}