#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/Writer_OFF.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygonal_surface_reconstruction.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/Timer.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel			Kernel;

typedef Kernel::Point_3												Point;
typedef Kernel::Vector_3											Vector;
typedef std::pair<Point, Vector>									Point_with_normal;
typedef std::vector<Point_with_normal>								Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>			Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal>		Normal_map;

typedef CGAL::Shape_detection_3::Shape_detection_traits<Kernel,	Pwn_vector, Point_map, Normal_map>	Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>			Efficient_ransac;
typedef CGAL::Shape_detection_3::Plane<Traits>						Plane;

typedef CGAL::Shape_detection_3::Plane_map<Traits>					Plane_map;
typedef CGAL::Shape_detection_3::Point_to_shape_index_map<Traits>	Point_to_plane_index_map;

typedef	CGAL::Polygonal_surface_reconstruction<Kernel>				Polygonal_surface_reconstruction;
typedef CGAL::Surface_mesh<Point>									Surface_mesh;

/*
* The following example first extracts planes from the input point
* cloud and then reconstructs the surface model.
*/

int main()
{
	// Points with normals.
	Pwn_vector points;

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

	// Shape detection
	Efficient_ransac ransac;
	ransac.set_input(points);
	ransac.add_shape_factory<Plane>();

	std::cout << "Extracting planes...";
	t.reset();
	ransac.detect();

	Efficient_ransac::Plane_range planes = ransac.planes();
	std::size_t num_planes = planes.size();
	std::cout << " Done. " << num_planes << " planes extracted. Time: " << t.time() << " sec." << std::endl;

	//////////////////////////////////////////////////////////////////////////

	std::cout << "Generating candidate faces...";
	t.reset();

	Polygonal_surface_reconstruction algo(
		points,
		CGAL::parameters::point_map(Point_map()).normal_map(Normal_map()).
		plane_map(Plane_map()).plane_index_map(Point_to_plane_index_map(points, planes))
	);

	std::cout << " Done. Time: " << t.time() << " sec." << std::endl;

	//////////////////////////////////////////////////////////////////////////

	Surface_mesh model;

	std::cout << "Reconstructing...";
	t.reset();

	if (!algo.reconstruct(model)) {
		std::cerr << " Failed: " << algo.error_message() << std::endl;
		return EXIT_FAILURE;
	}

    const std::string& output_file("data/output/cube_result.off");
    std::ofstream output_stream(output_file.c_str());
    if (output_stream && CGAL::write_off(output_stream, model))
		std::cout << " Done. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;
	else {
        std::cerr << " Failed saving file." << std::endl;
        return EXIT_FAILURE;
    }

	return EXIT_SUCCESS;
}
