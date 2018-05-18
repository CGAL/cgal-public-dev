#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygonal_surface_reconstruction.h>


typedef CGAL::Simple_cartesian<double>			Kernel;
typedef Kernel::Vector_3						Vector;

typedef	CGAL::Polygonal_surface_reconstruction<Kernel>				Polygonal_surface_reconstruction;
typedef Polygonal_surface_reconstruction::Point_set					Point_set;
typedef Point_set::Property_map<Vector>								Normal_map;
typedef Polygonal_surface_reconstruction::Point_set_with_segments	Point_set_with_segments;
typedef Polygonal_surface_reconstruction::Surface_mesh				Surface_mesh;


/*
 * The following example first extracts planes from the input point 
 * cloud and then reconstructs the surface model.
 */

int main()
{
	Point_set point_set;
	Normal_map normal_map;

	// read point set from file
	// ...

	Surface_mesh output_mesh;

	Polygonal_surface_reconstruction algo;
	if (!algo.reconstruct(point_set, normal_map, output_mesh)) {
		std::cout << "reconstruction failed" << std::endl;
		return 1;
	}

	return 0;
}
