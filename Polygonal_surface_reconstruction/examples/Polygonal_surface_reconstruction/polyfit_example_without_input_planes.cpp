#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/Writer_OFF.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygonal_surface_reconstruction.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/algo/point_set_with_segments.h>
#include <CGAL/Timer.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel		Kernel;

typedef Kernel::Point_3											Point;
typedef Kernel::Vector_3										Vector;
typedef std::pair<Point, Vector>								Point_with_normal;
typedef std::vector<Point_with_normal>							Point_list;
typedef CGAL::First_of_pair_property_map<Point_with_normal>		Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal>	Normal_map;

typedef CGAL::Planar_segment<Kernel>							Planar_segment;
typedef CGAL::Point_set_with_segments<Kernel>					Point_set_with_segments;

typedef	CGAL::Polygonal_surface_reconstruction<Kernel>			Polygonal_surface_reconstruction;
typedef CGAL::Surface_mesh<Point>								Surface_mesh;


/*
* The following example first extracts planes from the input point
* cloud and then reconstructs the surface model.
*/

bool extract_planes(
	const Point_list& points,
	Point_set_with_segments& segments)
{
	typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
	typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
	typedef CGAL::Shape_detection_3::Shape_detection_traits
		<Kernel, Point_list, Point_map, Normal_map>              Traits;
	typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>    Efficient_ransac;
	typedef CGAL::Shape_detection_3::Plane<Traits>               Plane;

	// Instantiates the Efficient_ransac engine.
	Efficient_ransac ransac;

	// Provides the input data.
	// Why does RANSAC require a non-const input?
	Point_list* point_set = const_cast<Point_list*>(&points);
	ransac.set_input(*point_set);

	// Registers planar shapes via template method.
	ransac.template add_shape_factory<Plane>();

	// Detects registered shapes with default parameters.
	ransac.detect();

	const typename Efficient_ransac::Shape_range& shapes = ransac.shapes();

	// update the output.

	// clean first
	segments.clear();

	// upload points
	std::size_t num = points.size();
	for (std::size_t i = 0; i < num; ++i)
		segments.insert(points[i].first);

	// ignore the colors
	// ignore the normals

	// now the planar segments

	typename Efficient_ransac::Shape_range::const_iterator it = shapes.begin();
	for (; it != shapes.end(); ++it) {
		boost::shared_ptr<typename Efficient_ransac::Shape> shape = *it;
		const std::vector<std::size_t>& indices = (*it)->indices_of_assigned_points();
		Planar_segment* s = new Planar_segment;
		s->set_point_set(&segments);
		s->insert(s->end(), indices.begin(), indices.end());
		s->fit_supporting_plane();
		segments.planar_segments().push_back(s);
	}

	return !shapes.empty();
}



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

	Point_set_with_segments point_set_with_planes;
	
	std::cout << "Extracting planes...";
	t.reset();

	extract_planes(points, point_set_with_planes);

	std::size_t num_planes = point_set_with_planes.planar_segments().size();
	std::cout << " Done. " << num_planes << " planes extracted. Time: " << t.time() << " sec." << std::endl;

	std::cout << "Reconstructing...";
	t.reset();

	Polygonal_surface_reconstruction algo(point_set_with_planes);

	Surface_mesh model;
	if (!algo.reconstruct(model)) {
		std::cerr << " Failed." << std::endl;
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
