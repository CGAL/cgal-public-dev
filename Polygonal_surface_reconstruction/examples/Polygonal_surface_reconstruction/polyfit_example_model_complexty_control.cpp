#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_with_segments.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygonal_surface_reconstruction.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel		Kernel;

typedef CGAL::Point_set_with_segments<Kernel>					Point_set_with_segments;
typedef	CGAL::Polygonal_surface_reconstruction<Kernel>			Polygonal_surface_reconstruction;

typedef Kernel::Point_3											Point;
typedef CGAL::Surface_mesh<Point>								Surface_mesh;


/*
* The following example shows how to control the model complexity by
* tuning the weight of the model complexity term.
* This example uses CGAL::Polygonal_surface_reconstruction() class,
* so the plane extraction and candidate generation only need to be
* done once.
*/

int main()
{
	const std::string& data_file("data/building.vg");
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

	// model 2: less details
	if (!algo.select_faces(candidate_faces, model, 0.43, 0.27, 0.40)) {
		std::cerr << "failed in face selection" << std::endl;
		return 1;
	}
	else {
		std::string file_name = "reconstruction_2.obj";
		// save model...
	}

	// model 3: even less details
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
