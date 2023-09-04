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
void test_sphere()
{
    std::string fname = "sphere";
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
    
    qem::Variational_shape_reconstruction manager(
        pointset,
        generators,
        distance_weight,
        qem::VERBOSE_LEVEL::HIGH,
        qem::INIT_QEM_GENERATOR::RANDOM);
    manager.region_growing(steps);
    manager.reconstruction(dist_ratio, fitting, coverage, complexity);


    Pointset point_cloud = manager.get_point_cloud_clustered();
    Polyhedron mesh = manager.get_reconstructed();
    
    std::ofstream mesh_file;
    mesh_file.open("output/mesh_"+fname+".off");
    CGAL::write_off(mesh_file, mesh);
    mesh_file.close();

    std::ifstream test_file;
    test_file.open("output/mesh_"+fname+".off");
    // generated with seed 27
    std::vector<std::string> ground_truth
    {
        "OFF",
        "4 4 0",
        "",
        "-0.484022 0.638033 1.01132",
        "1.29048 0.0272195 0.0375765",
        "-0.42351 -1.21661 0.0284324",
        "-0.391264 0.558783 -1.08751",
        "3  0 1 2",
        "3  3 1 0",
        "3  3 2 1",
        "3  0 2 3"
    };
    std::string line;
    int id=0;
    bool pass = true;
    while(getline(test_file, line)) {
     if(id <11&& line!=ground_truth[id])
     {
        std::cout<< "didn't pass : " <<line<<"\n";
        pass = false;
     }
     id++;
     
    }
    if(pass)
    {
        std::cout<<" It's a tetrahedron"<<"\n";
    }
    else
        std::cout<<" It's not a tetrahedron"<<"\n";
    test_file.close();


}
int main(int argc, char** argv)
{	
    test_sphere();
	return 0;
}
