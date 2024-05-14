// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Variational_shape_reconstruction.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Point_set_3.h>

#include <iostream>
#include <CGAL/Point_set_3/IO/XYZ.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::FT                  FT;
typedef Kernel::Plane_3             Plane;
typedef Kernel::Point_3             Point;
typedef Kernel::Vector_3            Vector;
typedef Kernel::Point_2             Point_2;
typedef Kernel::Triangle_3          Triangle;
typedef Kernel::Segment_2           Segment_2;

typedef CGAL::First_of_pair_property_map<std::pair<Point, Vector> > Point_map;
typedef CGAL::Second_of_pair_property_map<std::pair<Point, Vector> > Normal_map;
typedef CGAL::Point_set_3< Point, Vector > Pointset;

int main()
{	
    Pointset pointset;
    if (!CGAL::IO::read_XYZ("data/sphere.xyz", pointset))
    { 
        std::cerr << "Error: cannot read file " << std::endl;
        return EXIT_FAILURE;
    } 

    size_t generators = 5; 
    const size_t steps = 50;

    const double split_threshold = 10e-2;
    const double distance_weight = 10e-5;
	
    qem::Variational_shape_reconstruction vsr(
        pointset,
        generators,
        distance_weight,
        qem::VERBOSE_LEVEL::HIGH,
        qem::INIT_QEM_GENERATOR::RANDOM);

    vsr.clustering(steps, split_threshold);
    vsr.reconstruction();
	auto mesh = vsr.get_reconstructed_mesh();

    std::ofstream mesh_file;
    mesh_file.open("mesh.off");
    CGAL::write_off(mesh_file, mesh);
    mesh_file.close();

	return EXIT_SUCCESS;
}

