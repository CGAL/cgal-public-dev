#include "scene.h"

#include <iostream>
#include <fstream>

#include <QFileInfo>
#include <QTextStream>
#include <QInputDialog>
#include <QtCore/qglobal.h>

#include <QGLViewer/qglviewer.h>

#include "function_sphere.h"
#include "function_torus.h"
#include "function_ellipsoid.h"

#include <CGAL/Eigen_solver_traits.h>

Scene::Scene(): tr(), m_c2t3(tr)
{
	// view options
	m_view_mesh = true;
	m_view_edges = true;
	m_view_vertices = true;
}

Scene::~Scene()
{
}

void Scene::render()
{
	/*if (m_c2t3.nb_vertices() == 0)
		return;

	::glDisable(GL_LIGHTING);

	if(m_view_vertices)
		m_c2t3.render_vertices(3.0f, 0, 0, 0);

	if(m_view_edges)
		m_c2t3.render_edges(1.0f, 128, 128, 128);
*/
	//if(m_view_mesh)
		//m_c2t3.render_cells(0.0);


}


void Scene::mesh_sphere(const FT angle,
	const FT sizing,
	const FT approximation)
{
	// Define domain (careful: Sphere_3 constructor uses squared radius)
	/*Sphere bounding_sphere(CGAL::ORIGIN, 1.5 * 1.5);
	Mesh_domain domain(function_sphere, bounding_sphere); // unit sphere

	// set mesh criteria
	Mesh_criteria criteria(facet_angle = angle,
		facet_size = sizing,
		facet_distance = approximation,
		cell_radius_edge = 2.0,
		cell_size = sizing);

	// Mesh generation
	m_c2t3.clear();
	std::cout << "Meshing a unit sphere...";
	m_c2t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
	std::cout << "done" << std::endl;*/

}


void Scene::mesh_torus(const FT angle,
	const FT sizing,
	const FT approximation)
{
}

void Scene::mesh_ellipsoid(const FT angle,
	const FT sizing,
	const FT approximation)
{
}

void Scene::read_xyz(QString filename){
	FT sm_angle = 20.0; // Min triangle angle (degrees).
	FT sm_radius = 100; // Max triangle size w.r.t. point set average spacing.
	FT sm_distance = 0.5; // Approximation error w.r.t. point set average spacing.

	// Accumulated errors
	int accumulated_fatal_err = EXIT_SUCCESS;
	CGAL::Timer task_timer; task_timer.start();
	PointList points;
	std::string input_filename = filename.toStdString();
	// If OFF file format
	std::cerr << "Open " << input_filename << " for reading..." << std::endl;
	std::string extension = input_filename.substr(input_filename.find_last_of('.'));
	if (extension == ".xyz" || extension == ".XYZ" ||
					 extension == ".pwn" || extension == ".PWN")
	{
		// Reads the point set file in points[].
		// Note: read_xyz_points_and_normals() requires an iterator over points
		// + property maps to access each point's position and normal.
		// The position property map can be omitted here as we use iterators over Point_3 elements.
		std::ifstream stream(input_filename.c_str());
		if (!stream ||
				!CGAL::read_xyz_points(
															stream,
															std::back_inserter(points),
															CGAL::parameters::normal_map
															(CGAL::make_normal_of_point_with_normal_map(PointList::value_type()))
															))
		{
			std::cerr << "Error: cannot read file " << input_filename << std::endl;
			accumulated_fatal_err = EXIT_FAILURE;
			//continue;
		}
	}
	else
	{
		std::cerr << "Error: cannot read file " << input_filename << std::endl;
		accumulated_fatal_err = EXIT_FAILURE;
		//continue;
	}
	// Prints status
	std::size_t memory = CGAL::Memory_sizer().virtual_size();
	std::size_t nb_points = points.size();
	std::cerr << "Reads file " << input_filename << ": " << nb_points << " points, "
																											<< task_timer.time() << " seconds, "
																											<< (memory>>20) << " Mb allocated"
																											<< std::endl;
	task_timer.reset();

	//***************************************
	// Checks requirements
	//***************************************

	if (nb_points == 0)
	{
		std::cerr << "Error: empty point set" << std::endl;
		accumulated_fatal_err = EXIT_FAILURE;
		//continue;
	}

	bool points_have_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
	if ( ! points_have_normals )
	{
		std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
		// this is not a bug => do not set accumulated_fatal_err
		//continue;
	}

	//CGAL::Timer reconstruction_timer; reconstruction_timer.start();
	CGAL::Timer reconstruction_timer; reconstruction_timer.start();

	//***************************************
	// Computes implicit function
	//***************************************

	std::cerr << "Computes Poisson implicit function...\n";

	// Creates implicit function from the read points.
	// Note: this method requires an iterator over points
	// + property maps to access each point's position and normal.
	// The position property map can be omitted here as we use iterators over Point_3 elements.
	Poisson_reconstruction_function function(
														points.begin(), points.end(),
														CGAL::make_normal_of_point_with_normal_map(PointList::value_type())
														);

	// Computes the Poisson indicator function f()
	// at each vertex of the triangulation.
	//typedef CGAL::Eigen_sparse_matrix<double>::EigenType EigenMatrix;

	//CGAL::Eigen_solver_traits< Eigen::BiCGSTAB<EigenMatrix> > solver;
	if ( ! function.compute_implicit_function() )
	{
		std::cerr << "Error: cannot compute implicit function" << std::endl;
		accumulated_fatal_err = EXIT_FAILURE;
		//continue;
	}

	// Prints status
	std::cerr << "Total implicit function (triangulation+refinement+solver): " << task_timer.time() << " seconds\n";
	task_timer.reset();

	//***************************************
	// Surface mesh generation
	//***************************************

	std::cerr << "Surface meshing...\n";

	// Computes average spacing
	FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points, 6 /* knn = 1 ring */);

	// Gets one point inside the implicit surface
	Point inner_point = function.get_inner_point();
	FT inner_point_value = function(inner_point);
	if(inner_point_value >= 0.0)
	{
		std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
		accumulated_fatal_err = EXIT_FAILURE;
		//continue;
	}

	// Gets implicit function's radius
	Sphere bsphere = function.bounding_sphere();
	FT radius = std::sqrt(bsphere.squared_radius());

	// Defines the implicit surface: requires defining a
	// conservative bounding sphere centered at inner point.
	FT sm_sphere_radius = 5.0 * radius;
	FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
	Surface_3 surface(function,
										Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
										sm_dichotomy_error/sm_sphere_radius);

	// Defines surface mesh generation criteria
	CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
																											sm_radius*average_spacing,  // Max triangle size
																											sm_distance*average_spacing); // Approximation error

	CGAL_TRACE_STREAM << "  make_surface_mesh(sphere center=("<<inner_point << "),\n"
										<< "                    sphere radius="<<sm_sphere_radius<<",\n"
										<< "                    angle="<<sm_angle << " degrees,\n"
										<< "                    triangle size="<<sm_radius<<" * average spacing="<<sm_radius*average_spacing<<",\n"
										<< "                    distance="<<sm_distance<<" * average spacing="<<sm_distance*average_spacing<<",\n"
										<< "                    dichotomy = distance/"<<sm_distance*average_spacing/sm_dichotomy_error<<",\n"
										<< "                    Manifold_with_boundary_tag)\n";

	// Generates surface mesh with manifold option
	STr tr; // 3D Delaunay triangulation for surface mesh generation
	C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
	CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
													surface,                              // implicit surface
													criteria,                             // meshing criteria
													CGAL::Manifold_with_boundary_tag());  // require manifold mesh

	// Prints status
	/*long*/ memory = CGAL::Memory_sizer().virtual_size();
	std::cerr << "Surface meshing: " << task_timer.time() << " seconds, "
																	 << tr.number_of_vertices() << " output vertices, "
																	 << (memory>>20) << " Mb allocated"
																	 << std::endl;
	task_timer.reset();

	if(tr.number_of_vertices() == 0) {
		accumulated_fatal_err = EXIT_FAILURE;
		//continue;
	}

	// Converts to polyhedron
	Polyhedron output_mesh;
	CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, output_mesh);

	// Prints total reconstruction duration
	std::cerr << "Total reconstruction (implicit function + meshing): " << reconstruction_timer.time() << " seconds\n";



	std::cerr << std::endl;

	// Returns accumulated fatal error
	std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
}
