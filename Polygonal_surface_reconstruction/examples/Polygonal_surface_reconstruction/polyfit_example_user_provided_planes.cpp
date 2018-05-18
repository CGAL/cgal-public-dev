#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygonal_surface_reconstruction.h>


typedef CGAL::Simple_cartesian<double>								Kernel;

typedef	CGAL::Polygonal_surface_reconstruction<Kernel>				Polygonal_surface_reconstruction;
typedef Polygonal_surface_reconstruction::Point_set_with_segments	Point_set_with_segments;
typedef Polygonal_surface_reconstruction::Surface_mesh				Surface_mesh;


/*
* The following example shows the reconstruction using user-provided
* planar segments.
*/

int main()
{
	Point_set_with_segments point_set;
	// read point set and planar segments from file
	// ...

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
