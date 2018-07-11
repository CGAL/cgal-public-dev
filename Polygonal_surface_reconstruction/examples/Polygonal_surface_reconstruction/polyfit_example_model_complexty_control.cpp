#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/Writer_OFF.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygonal_surface_reconstruction.h>
#include <CGAL/Timer.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel		Kernel;
typedef Kernel::Point_3											Point;
typedef Kernel::Vector_3										Vector;
typedef	CGAL::Polygonal_surface_reconstruction<Kernel>			Polygonal_surface_reconstruction;
typedef CGAL::Surface_mesh<Point>								Surface_mesh;

// Point with normal, and plane index
typedef boost::tuple<Point, Vector, int>						PNI;
typedef CGAL::Nth_of_tuple_property_map<0, PNI>					Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNI>					Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNI>					Plane_index_map;


/*
* The following example shows how to control the model complexity by
* increasing the influence of the model complexity term.
* In this example, the intermediate results from plane extraction and 
* candidate generation are cached and reused.
*/

int main()
{
	const std::string& input_file("data/building-ascii.ply");
	std::ifstream input_stream(input_file.c_str());

	std::vector<PNI> points; // store points

	std::cout << "Loading point cloud: " << input_file << "...";
	CGAL::Timer t;
	t.start();

	if (!input_stream ||
		!CGAL::read_ply_points_with_properties(
			input_stream,
			std::back_inserter(points),
			CGAL::make_ply_point_reader(Point_map()),
			CGAL::make_ply_normal_reader(Normal_map()),
			std::make_pair(Plane_index_map(), CGAL::PLY_property<int>("segment_index"))))
	{
		std::cerr << "Error: cannot read file " << input_file << std::endl;
		return EXIT_FAILURE;
	}
	else
		std::cout << " Done. " << points.size() << " points. Time: " << t.time() << " sec." << std::endl;

	//////////////////////////////////////////////////////////////////////////

	std::cout << "Generating candidate faces...";
	t.reset();

	Polygonal_surface_reconstruction algo(
		points,
		Point_map(),
		Normal_map(),
		Plane_index_map()
	);

	std::cout << " Done. Time: " << t.time() << " sec." << std::endl;

	//////////////////////////////////////////////////////////////////////////

	// reconstruction with complexity control
	
	// model 1: more details
	Surface_mesh model;

	std::cout << "Reconstructing with complexity = 0.2...";
	t.reset();
	if (!algo.reconstruct(model, 0.43, 0.27, 0.2)) {
		std::cerr << " Failed: " << algo.error_message() << std::endl;
		return EXIT_FAILURE;
	}
	else {
       const std::string& output_file = "data/output/building_result_complexity-0.2.off";
       std::ofstream output_stream(output_file.c_str());
       if (output_stream && CGAL::write_off(output_stream, model))
			std::cout << " Done. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;
		else {
           std::cerr << " Failed saving file." << std::endl;
           return EXIT_FAILURE;
       }
	}

	// model 2: less details
	std::cout << "Reconstructing with complexity = 0.5...";
	t.reset();
	if (!algo.reconstruct(model, 0.43, 0.27, 0.5)) {
		std::cerr << " Failed: " << algo.error_message() << std::endl;
		return EXIT_FAILURE;
	}
	else {
       const std::string& output_file = "data/output/building_result_complexity-0.5.off";
       std::ofstream output_stream(output_file.c_str());
       if (output_stream && CGAL::write_off(output_stream, model))
			std::cout << " Done. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;
		else {
           std::cerr << " Failed saving file." << std::endl;
           return EXIT_FAILURE;
       }
	}

	return EXIT_SUCCESS;
}
