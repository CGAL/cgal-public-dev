#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygonal_surface_reconstruction.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel		Kernel;

typedef Kernel::Point_3											Point;
typedef Kernel::Vector_3										Vector;
typedef std::pair<Point, Vector>								Point_with_normal;
typedef std::vector<Point_with_normal>							Point_list;
typedef CGAL::First_of_pair_property_map<Point_with_normal>		Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal>	Normal_map;

typedef CGAL::Polygonal_surface_reconstruction_traits<Kernel, Point_list, Point_map, Normal_map>	Traits;
typedef	CGAL::Polygonal_surface_reconstruction<Traits>			Polygonal_surface_reconstruction;

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
	std::ifstream stream("data/ball.pwn");
	if (!stream || !CGAL::read_xyz_points(stream,
		std::back_inserter(points),
		CGAL::parameters::point_map(Point_map()).normal_map(Normal_map())))
	{
		std::cerr << "Error: cannot read file data/ball.pwn" << std::endl;
		return EXIT_FAILURE;
	}

	Surface_mesh output_mesh;

	Polygonal_surface_reconstruction algo;
	if (!algo.reconstruct(points, output_mesh))
	{
		std::cerr << "reconstruction failed" << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
