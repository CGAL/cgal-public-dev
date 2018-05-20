#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Point_set_with_segments.h>
#include <CGAL/Polygonal_surface_reconstruction.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel		Kernel;
typedef CGAL::Point_set_with_segments<Kernel>					Point_set_with_segments;

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
* The following example shows the reconstruction using user-provided
* planar segments.
*/

int main()
{
	const std::string& data_file("data/apartment.vg");
	Point_set_with_segments point_set;

	if (!point_set.read(data_file)) {
		std::cerr << "failed reading data file \'" << data_file << "\'" << std::endl;
		return 1;
	}

	Polygonal_surface_reconstruction algo;

	Surface_mesh candidate_faces;
	if (!algo.generate_candidate_faces(point_set, candidate_faces)) {
		std::cerr << "failed to generate candidate faces" << std::endl;
		return 1;
	}

	Surface_mesh output_mesh;
	if (!algo.select_faces(candidate_faces, output_mesh)) {
		std::cerr << "failed in face selection" << std::endl;
		return 1;
	}

	// save the mesh model
	// ...

	return 0;
}
