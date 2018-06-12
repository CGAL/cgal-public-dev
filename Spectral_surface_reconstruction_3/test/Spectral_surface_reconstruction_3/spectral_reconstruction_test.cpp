// spectral_reconstruction_test.cpp

//----------------------------------------------------------
// Test the Spectral Delaunay Reconstruction method:
// For each input point set or mesh's set of vertices, reconstruct a surface.
// No output.
//----------------------------------------------------------
// spectral_reconstruction_test mesh1.off point_set2.xyz...

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/trace.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Spectral_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>

#include <deque>
#include <cstdlib>
#include <fstream>
#include <math.h>

#include <boost/foreach.hpp>
#include "boost/program_options.hpp"

#include <CGAL/disable_warnings.h>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::deque<Point_with_normal> PointList;

// polyhedron
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

// Spectral implicit function
typedef CGAL::Spectral_reconstruction_function<Kernel> Spectral_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Spectral_reconstruction_function> Surface_3;


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Test the Spectral Delaunay Reconstruction method" << std::endl;

  //***************************************
  // decode parameters
  //***************************************

  namespace po = boost::program_options;
  po::options_description desc("Options");
  // Spectral options
  desc.add_options()
      ("help,h", "Display this help message")
      ("input,i", po::value<std::vector<std::string> >(), "Input files")
      ("output,o", po::value<std::string>()->default_value("out.off"), "The suffix of the output files")
      ("bilaplacian,b", po::value<double>()->default_value(0.5), "The global bilaplacian coefficient")
      ("laplacian,l", po::value<double>()->default_value(1.), "The global laplacian coefficient")
      ("ratio,r", po::value<double>()->default_value(10.), "The largest eigenvalue of the tensor C")
      ("fitting,f", po::value<double>()->default_value(1.), "The data fitting term")
      ("sm_angle,a", po::value<double>()->default_value(20.), "The min triangle angle (degrees).")
      ("sm_radius,s", po::value<double>()->default_value(2.), "The max triangle size w.r.t. point set average spacing.")
      ("sm_distance,d", po::value<double>()->default_value(1), "The approximation error w.r.t. point set average spacing.");

  // Parse input files
  po::positional_options_description p;
  p.add("input", -1);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);

  if(vm.count("help")){
    std::cerr << "For each input point set or mesh's set of vertices, reconstruct a surface.\n";
    std::cerr << "\n\n";
    std::cerr << desc;
    return EXIT_FAILURE;
  }

  std::vector<std::string> files;

  if(vm.count("input")){
    files = vm["input"].as<std::vector<std::string> >();
    for(std::string file : files){
        std::cerr << "Input file: " << file << std::endl;
    }
  }
  else{
    std::cerr << "No input file specified." << std::endl;
    return EXIT_FAILURE;
  }

  double sm_angle = vm["sm_angle"].as<double>();
  double sm_radius = vm["sm_radius"].as<double>();
  double sm_distance = vm["sm_distance"].as<double>();

  double laplacian = vm["laplacian"].as<double>();
  double bilaplacian = vm["bilaplacian"].as<double>();
  double ratio = vm["ratio"].as<double>();
  double fitting = vm["fitting"].as<double>();

  std::string outfile = vm["output"].as<std::string>();

  // Accumulated errors
  int accumulated_fatal_err = EXIT_SUCCESS;

  // Process each input file
  for (int i = 1; i <= files.size(); i++)
  {
    CGAL::Timer task_timer; task_timer.start();

    std::cerr << std::endl;

    //***************************************
    // Loads mesh/point set
    //***************************************

    // File name is:
    std::string input_filename  = files[i - 1];

    PointList points;

    // If OFF file format
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;
    std::string extension = input_filename.substr(input_filename.find_last_of('.'));
    if (extension == ".off" || extension == ".OFF")
    {
      // Reads the mesh file in a polyhedron
      std::ifstream stream(input_filename.c_str());
      Polyhedron input_mesh;
      CGAL::scan_OFF(stream, input_mesh, true /* verbose */);
      if(!stream || !input_mesh.is_valid() || input_mesh.empty())
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        accumulated_fatal_err = EXIT_FAILURE;
        continue;
      }

      // Converts Polyhedron vertices to point set.
      // Computes vertices normal from connectivity.
      BOOST_FOREACH(boost::graph_traits<Polyhedron>::vertex_descriptor v, 
                    vertices(input_mesh)){
        const Point& p = v->point();
        Vector n = CGAL::Polygon_mesh_processing::compute_vertex_normal(v,input_mesh);
        points.push_back(Point_with_normal(p,n));
      }
    }
    // If XYZ file format
    else if (extension == ".xyz" || extension == ".XYZ" ||
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
        continue;
      }
    }
    else
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
      continue;
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
      continue;
    }

    bool points_have_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
    if ( ! points_have_normals )
    {
      std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
      // this is not a bug => do not set accumulated_fatal_err
      continue;
    }

    CGAL::Timer reconstruction_timer; reconstruction_timer.start();

    //***************************************
    // Computes implicit function
    //***************************************

    std::cerr << "Computes Spectral implicit function...\n";

    // Creates implicit function from the read points.
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    Spectral_reconstruction_function function(
                              points.begin(), points.end(),
                              CGAL::make_normal_of_point_with_normal_map(PointList::value_type())
                              , CGAL::Default_property_map<PointList::iterator, FT>(1.)
                              );
  
    // Computes the Spectral indicator function f()
    // at each vertex of the triangulation.
    if ( ! function.compute_implicit_function(bilaplacian, laplacian, fitting, ratio) )
    {
      std::cerr << "Error: cannot compute implicit function" << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
     continue;
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
      continue;
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
      continue;
    }

    // Converts to polyhedron
    Polyhedron output_mesh;
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, output_mesh);

    // Prints total reconstruction duration
    std::cerr << "Total reconstruction (implicit function + meshing): " << reconstruction_timer.time() << " seconds\n";

    // Output the 3D complex to an OFF file. 
    std::ofstream out(std::to_string(i) + "_" + outfile); 
    out << output_mesh;

  } // for each input file

  std::cerr << std::endl;

  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;

  return 0;
}
