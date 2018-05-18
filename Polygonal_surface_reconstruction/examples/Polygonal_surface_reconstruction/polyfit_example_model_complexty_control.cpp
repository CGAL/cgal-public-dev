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
* The following example shows how to control the model complexity by
* tuning the weight of the model complexity term.
* This example uses CGAL::Polygonal_surface_reconstruction() class,
* so the plane extraction and candidate generation only need to be
* done once.
*/

int main()
{
	Point_set point_set;
	Normal_map normal_map;
	// read point set from file
	// ...


	Point_set_with_segments segments;
	Polygonal_surface_reconstruction algo;
	if (!algo.extract_planes(point_set, normal_map, segments)) {
		std::cerr << "failed to extract planar segments" << std::endl;
		return 1;
	}
	
	Surface_mesh candidate_faces;
	if (!algo.generate_candidate_faces(segments, candidate_faces)) {
		std::cerr << "failed to generate candidate faces" << std::endl;
		return 1;
	}

	// increasing the weight of the model complexity term results in less detailed 3D models.
	
	// model 1: most details
	Surface_mesh model;
	if (!algo.select_faces(candidate_faces, model, 0.43, 0.27, 0.20)) {
		std::cerr << "failed in face selection" << std::endl;
		return 1;
	}
	else {
		std::string file_name = "reconstruction_1.obj";
		// save model...
	}

	if (!algo.select_faces(candidate_faces, model, 0.43, 0.27, 0.40)) {
		std::cerr << "failed in face selection" << std::endl;
		return 1;
	}
	else {
		std::string file_name = "reconstruction_2.obj";
		// save model...
	}

	if (!algo.select_faces(candidate_faces, model, 0.43, 0.27, 0.60)) {
		std::cerr << "failed in face selection" << std::endl;
		return 1;
	}
	else {
		std::string file_name = "reconstruction_3.obj";
		// save model...
	}

	// model 4: least details
	if (!algo.select_faces(candidate_faces, model, 0.43, 0.27, 0.90)) {
		std::cerr << "failed in face selection" << std::endl;
		return 1;
	}
	else {
		std::string file_name = "reconstruction_4.obj";
		// save model...
	}

	return 0;
}
