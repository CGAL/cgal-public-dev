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

void test_generators_qem()
{
    std::string fname = "piece_meca";
    Pointset pointset;
    if (!CGAL::IO::read_XYZ( "../data/"+fname+".xyz",pointset))
    {
        std::cerr << "Error: cannot read file " << std::endl;
        return;
    } 
    const size_t generators = 100;
    const size_t steps = 10;
     const double distance_weight =0.00001;
	qem::Variational_shape_reconstruction manager(pointset,generators, distance_weight );
    manager.region_growing(steps);
    manager.write_csv();
    // todo check evolution generators
}
int main(int argc, char** argv)
{	
    test_generators_qem();
	return 0;
}
