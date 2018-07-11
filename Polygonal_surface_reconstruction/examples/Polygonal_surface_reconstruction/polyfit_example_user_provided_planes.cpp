#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Writer_OFF.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/property_map.h>
#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygonal_surface_reconstruction.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel		Kernel;
typedef Kernel::Point_3											Point;
typedef Kernel::Vector_3										Vector;
typedef	CGAL::Polygonal_surface_reconstruction<Kernel>			Polygonal_surface_reconstruction;
typedef CGAL::Surface_mesh<Point>								Surface_mesh;

// Point with normal, and plane index
typedef CGAL::cpp11::tuple<Point, Vector, int>					PNI;
typedef CGAL::Nth_of_tuple_property_map<0, PNI>					Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNI>					Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNI>					Plane_index_map;

/*
* The following example shows the reconstruction using user-provided
* planar segments stored in PLY format. In the PLY format, a property
* named "segment_index" stores the plane index for each point (-1 if 
* the point is not assigned to a plane). 
*/

int main()
{
	const std::string& input_file("data/ball-ascii.ply");
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

	//// Display points read
	//for (std::size_t i = 0; i < points.size(); ++i)
	//{
	//	const Point& p = get<0>(points[i]);
	//	const Vector& n = get<1>(points[i]);
	//	int idx = get<2>(points[i]);
	//	std::cerr << "Point (" << p << ") with normal (" << n	<< "), and plane index " << idx << std::endl;
	//}

	std::cout << "Generating candidate faces...";
	t.reset();

	Polygonal_surface_reconstruction algo(
		points,
		CGAL::parameters::point_map(Point_map()).normal_map(Normal_map()).
		plane_index_map(Plane_index_map())
	);

	std::cout << " Done. Time: " << t.time() << " sec." << std::endl;

	Surface_mesh model;
	std::cout << "Optimizing...";
	t.reset();
	if (!algo.reconstruct(model)) {
		std::cerr << " Failed: " << algo.error_message() << std::endl;
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
