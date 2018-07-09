#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygonal_surface_reconstruction.h>
#include <CGAL/IO/Writer_OFF.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/property_map.h>
#include <CGAL/Timer.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel		Kernel;
typedef Kernel::Point_3											Point;
typedef Kernel::Vector_3										Vector;
typedef CGAL::Surface_mesh<Point>								Surface_mesh;

// Point with normal, color and intensity
typedef CGAL::cpp11::tuple<Point, Vector, int> PNCI;
typedef CGAL::Nth_of_tuple_property_map<0, PNCI> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNCI> Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNCI> Color_map;
typedef CGAL::Nth_of_tuple_property_map<3, PNCI> Intensity_map;


/*
* The following example shows the reconstruction using user-provided
* planar segments.
*/

int main()
{
	const std::string& input_file("data/ball.ply");
	std::cout << "Loading point cloud: " << input_file << "...";
	CGAL::Timer t;
	t.start();

	std::ifstream input_stream(input_file.c_str());

	// Point with normal, color and intensity
	typedef CGAL::cpp11::tuple<Point, Vector, Color, int> PNCI;
	typedef CGAL::Nth_of_tuple_property_map<0, PNCI> Point_map;
	typedef CGAL::Nth_of_tuple_property_map<1, PNCI> Normal_map;
	typedef CGAL::Nth_of_tuple_property_map<2, PNCI> Color_map;
	typedef CGAL::Nth_of_tuple_property_map<3, PNCI> Intensity_map;


	if (!CGAL::read_point_set_with_segments(input_stream, point_set)) {
		std::cerr << " Failed." << std::endl;
		return EXIT_FAILURE;
	}
	else {
		std::cout << " Done. "
			<< point_set.number_of_points() << " points, "
			<< point_set.planar_segments().size() << " planar segments. Time: " << t.time() << " sec." << std::endl;
	}

	std::cout << "Generating candidate faces...";
	t.reset();

	Polygonal_surface_reconstruction algo(point_set);

	std::cout << " Done. Time: " << t.time() << " sec." << std::endl;

	Surface_mesh model;
	std::cout << "Optimizing...";
	t.reset();
	if (!algo.reconstruct(model)) {
		std::cerr << " Failed." << std::endl;
		return EXIT_FAILURE;
	}

	// save the mesh model
    const std::string& output_file("data/output/ball_result.off");
    std::ofstream output_stream(output_file.c_str());
    if (output_stream && CGAL::write_off(output_stream, model))
		std::cout << " Done. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;
	else {
        std::cerr << " Failed saving file." << std::endl;
        return EXIT_FAILURE;
    }

	return EXIT_SUCCESS;
}
