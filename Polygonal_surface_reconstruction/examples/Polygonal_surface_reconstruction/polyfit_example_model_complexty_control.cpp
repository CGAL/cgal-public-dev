#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/algo/point_set_with_planes.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygonal_surface_reconstruction.h>
#include <CGAL/IO/read_point_set_with_segments.h>
#include <CGAL/Timer.h>

#include <fstream>


//typedef CGAL::Exact_predicates_inexact_constructions_kernel		Kernel;

//typedef CGAL::Point_set_with_planes<Kernel>					Point_set_with_planes;
//typedef	CGAL::Polygonal_surface_reconstruction<Kernel>			Polygonal_surface_reconstruction;

//typedef Kernel::Point_3											Point;
//typedef CGAL::Surface_mesh<Point>								Surface_mesh;


/*
* The following example shows how to control the model complexity by
* increasing the influence of the model complexity term.
* In this example, the intermediate results from plane extraction and 
* candidate generation are cached and reused.
*/

int main()
{
//	Point_set_with_planes point_set;

//	const std::string& input_file("data/building.vg");
//	std::cout << "Loading point cloud: " << input_file << "...";
//	CGAL::Timer t;
//	t.start();

//	std::ifstream input_stream(input_file.c_str());
//	if (!CGAL::read_point_set_with_segments(input_stream, point_set)) {
//		std::cerr << " Failed." << std::endl;
//		return EXIT_FAILURE;
//	}
//	else {
//		std::cout << " Done. "
//			<< point_set.number_of_points() << " points, "
//			<< point_set.planar_segments().size() << " planar segments. Time: " << t.time() << " sec." << std::endl;
//	}

//	std::cout << "Generating candidate faces...";
//	t.reset();

//	Polygonal_surface_reconstruction algo(point_set);

//	std::cout << " Done. Time: " << t.time() << " sec." << std::endl;

//	// reconstruction with complexity control
	
//	// model 1: more details
//	Surface_mesh model;
//	std::cout << "Optimizing with complexity = 0.2...";
//	t.reset();
//	if (!algo.reconstruct(model, 0.43, 0.27, 0.2)) {
//		std::cerr << " Failed." << std::endl;
//		return EXIT_FAILURE;
//	}
//	else {
//        const std::string& output_file = "data/output/building_result_complexity-0.2.off";
//        std::ofstream output_stream(output_file.c_str());
//        if (output_stream && CGAL::write_off(output_stream, model))
//			std::cout << " Done. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;
//		else {
//            std::cerr << " Failed saving file." << std::endl;
//            return EXIT_FAILURE;
//        }
//	}

//	// model 2: less details
//	std::cout << "Optimizing with complexity = 0.5...";
//	t.reset();
//	if (!algo.reconstruct(model, 0.43, 0.27, 0.5)) {
//		std::cerr << " Failed." << std::endl;
//		return EXIT_FAILURE;
//	}
//	else {
//        const std::string& output_file = "data/output/building_result_complexity-0.5.off";
//        std::ofstream output_stream(output_file.c_str());
//        if (output_stream && CGAL::write_off(output_stream, model))
//			std::cout << " Done. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;
//		else {
//            std::cerr << " Failed saving file." << std::endl;
//            return EXIT_FAILURE;
//        }
//	}

	return EXIT_SUCCESS;
}
