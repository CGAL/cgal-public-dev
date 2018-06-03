#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/Writer_OFF.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygonal_surface_reconstruction.h>
#include <CGAL/Timer.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel		Kernel;

typedef Kernel::Point_3											Point;
typedef Kernel::Vector_3										Vector;
typedef std::pair<Point, Vector>								Point_with_normal;
typedef std::vector<Point_with_normal>							Point_list;
typedef CGAL::First_of_pair_property_map<Point_with_normal>		Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal>	Normal_map;


typedef	CGAL::Polygonal_surface_reconstruction<Kernel>			Polygonal_surface_reconstruction;
typedef CGAL::Surface_mesh<Point>								Surface_mesh;


/*
* The following example first extracts planes from the input point
* cloud and then reconstructs the surface model.
*/

int main()
{
	// Points with normals.
	Point_list points;

	// Loads point set from a file. 
	const std::string& input_file("data/cube.pwn");
    std::ifstream input_stream(input_file.c_str());
	std::cout << "Loading point cloud: " << input_file << "...";

	CGAL::Timer t;
	t.start();
    if (!input_stream ||
            !CGAL::read_xyz_points(input_stream,
		std::back_inserter(points),
		CGAL::parameters::point_map(Point_map()).normal_map(Normal_map())))
	{
		std::cerr << " Failed." << std::endl;
		return EXIT_FAILURE;
	}
	else
		std::cout << " Done. " << points.size() << " points. Time: " << t.time() << " sec." << std::endl;

	Surface_mesh model;
	Polygonal_surface_reconstruction algo;

	std::cout << "Reconstructing...";
	t.reset();
	if (!algo.reconstruct(points, model)) {
		std::cerr << " Failed." << std::endl;
		return EXIT_FAILURE;
	}

	const std::string& output_file("data/cube_result.off");
    std::ofstream output_stream(output_file.c_str());
    if (output_stream && CGAL::write_off(output_stream, model))
		std::cout << " Done. " << model.number_of_faces() << " faces. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;

	return EXIT_SUCCESS;
}
