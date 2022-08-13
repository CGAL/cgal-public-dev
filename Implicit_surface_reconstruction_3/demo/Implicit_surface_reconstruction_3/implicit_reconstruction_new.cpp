// implicit_reconstruction_test.cpp

//----------------------------------------------------------
// Test the Implicit Delaunay Reconstruction method:
// For each input point set or mesh's set of vertices, reconstruct a surface.
// No output.
//----------------------------------------------------------
// implicit_reconstruction_test -h

#undef min
#undef max

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Implicit_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>


#include <list>
#include <cstdlib>
#include <fstream>
#include <math.h>

#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include "boost/program_options.hpp"

#include <CGAL/disable_warnings.h>

//Mesh
#include <CGAL/Surface_mesh.h>





// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
//typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef std::pair<Point, Vector> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::list<Point_with_normal> PointList;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// polyhedron
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

// Spectral implicit function
typedef CGAL::Implicit_reconstruction_function<Kernel, PointList, Normal_map> Implicit_reconstruction_function;

//typedef CGAL::Implicit_reconstruction_function<Kernel, PointList, CGAL::Identity_property_map<Point_with_normal> > Implicit_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Implicit_reconstruction_function> Surface_3;

//Mesh
typedef CGAL::Surface_mesh<Point>                           Mesh;

// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Test the Implicit Delaunay Reconstruction method" << std::endl;

  //***************************************
  // decode parameters
  //***************************************

  namespace po = boost::program_options;
  po::options_description desc("Options");
  // Implicit options
  desc.add_options()
      ("help,h", "Display this help message")
      ("input,i", po::value<std::vector<std::string> >(), "Input files")
      ("output,o", po::value<std::string>()->default_value("out.off"), "The suffix of the output files")
      ("bilaplacian,b", po::value<double>()->default_value(100.), "The global bilaplacian coefficient")
      ("laplacian,l", po::value<double>()->default_value(1), "The global laplacian coefficient")
      ("hessian", po::value<double>()->default_value(1e-3), "The global hessian coefficient")
      ("ratio,r", po::value<double>()->default_value(10.), "The largest eigenvalue of the tensor C")
      ("fitting,f", po::value<double>()->default_value(10), "The data fitting term")
      ("size,s", po::value<int>()->default_value(500), "The number of slice")
      ("x", po::value<double>()->default_value(0), "The chosen x coordinate")
      ("isovalue,a", po::value<double>()->default_value(0.), "The isovalue to extract")
      ("octree,c", po::bool_switch()->default_value(false), "Use Octree Refinement / Delaunay Refinement")
      ("octree_debug,d", po::bool_switch()->default_value(false), "Octree Refinement Debug Mode ")
      ("poisson", po::bool_switch()->default_value(false), "Use spectral (false) / poisson (true)")
      ("ssd", po::bool_switch()->default_value(false), "Use spectral (false) / ssd (true)")
      ("vals,v", po::bool_switch()->default_value(false), "Save function value for all points in a ply file (true/false)")
      ("marching,m", po::bool_switch()->default_value(false), "Use marching tet to reconstruct surface")
      ("sm_angle", po::value<double>()->default_value(20.), "The min triangle angle (degrees).")
      ("sm_radius", po::value<double>()->default_value(100.), "The max triangle size w.r.t. point set average spacing.")
      ("sm_distance", po::value<double>()->default_value(0.25), "The approximation error w.r.t. point set average spacing.")
      ("scale", po::value<double>()->default_value(1), "scaling for testing.");;

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
  double hessian = vm["hessian"].as<double>();
  double ratio = vm["ratio"].as<double>();
  double fitting = vm["fitting"].as<double>();
  double isovalue = vm["isovalue"].as<double>();
  double scale = vm["scale"].as<double>();

  bool flag_poisson = vm["poisson"].as<bool>();
  bool flag_ssd = vm["ssd"].as<bool>();
  //bool flag_vals = vm["vals"].as<bool>();
  bool flag_marching = vm["marching"].as<bool>();
  bool flag_octree = vm["octree"].as<bool>();
  bool flag_octree_debug = vm["octree_debug"].as<bool>();

  int size = vm["size"].as<int>();
  double x = vm["x"].as<double>();

  std::string outfile = vm["output"].as<std::string>();

  // Accumulated errors
  int accumulated_fatal_err = EXIT_SUCCESS;

  // Process each input file
  for (int i = 1; i <= (int)files.size(); i++)
  {
    CGAL::Timer task_timer; task_timer.start();

    std::cerr << std::endl;

    //***************************************
    // Loads mesh/point set
    //***************************************

    // File name is:
    std::string input_filename  = files[i - 1];

    PointList points_buf;

    // If OFF file format
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;

    if(!CGAL::IO::read_points(input_filename, std::back_inserter(points_buf),
                              CGAL::parameters::point_map(Point_map())
                                               .normal_map(Normal_map())))
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      accumulated_fatal_err = EXIT_FAILURE;
      continue;
    }

    // scaling
    PointList points;
    auto it = std::back_inserter(points);
    
    for (auto p = points_buf.begin(); p != points_buf.end(); p++)
    {
        Point pos = p->first;
        *it = std::make_pair(Point(pos.x() * scale, pos.y() * scale, pos.z() * scale), p->second);
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

    bool points_have_normals = (points.begin()->second != CGAL::NULL_VECTOR);
    if ( ! points_have_normals )
    {
      std::cerr << "Input point set not supported: this reconstruction method requires unoriented normals" << std::endl;
      // this is not a bug => do not set accumulated_fatal_err
      continue;
    }

    CGAL::Timer reconstruction_timer; reconstruction_timer.start();

    //***************************************
    // Computes implicit function
    //***************************************

    std::cerr << "Initializing class object...\n";

    // Creates empty implicit function
    Implicit_reconstruction_function function;
    
	// from the read points.
	// Note: this method requires an iterator over points
	// + property maps to access each point's position and normal.
	// The position property map can be omitted here as we use iterators over Point_3 elements.
	function.initialize_point_map(points, Point_map(), Normal_map(), flag_octree, flag_octree_debug);

    std::cerr << "Initialization: " << reconstruction_timer.time() << " seconds\n";
    
    // Computes the Implicit indicator function f()
    // at each vertex of the triangulation.
    std::cerr << "Computes implicit function...\n";

    if(flag_poisson)
	{
      if (! function.compute_poisson_implicit_function_new(fitting)){
        std::cerr << "Error: cannot compute implicit function" << std::endl;
        accumulated_fatal_err = EXIT_FAILURE;
        continue;
      }
    }
    else if (flag_ssd) {
        if (! function.compute_ssd_implicit_function(fitting, laplacian, hessian) )
        {
            std::cerr << "Error: cannot compute implicit function" << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            continue;
        }
    }
    else {
      if (! function.compute_spectral_implicit_function_test(bilaplacian, laplacian, fitting) )
      {
        std::cerr << "Error: cannot compute implicit function" << std::endl;
        accumulated_fatal_err = EXIT_FAILURE;
        continue;
      }
    }

    // Prints status
    std::cerr << "Total implicit function (triangulation+refinement+solver): " << task_timer.time() << " seconds\n";
    task_timer.reset();

    if(flag_marching)
    {
      Mesh output_mesh;
      std::cerr << "Marching tets..." << std::endl;
      std::string curr_outfile(std::to_string(i) + "_" + outfile);
      function.marching_tetrahedra(isovalue, output_mesh);
      std::ofstream out(curr_outfile);
      out << output_mesh;
    }
    else{
      Point inner_point = function.get_inner_point();
      Sphere bsphere = function.bounding_sphere();
      double radius = std::sqrt(bsphere.squared_radius());

      // Computes average spacing
      double spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points, 6, CGAL::parameters::point_map(Point_map()) /* knn = 1 ring */);
      double sm_sphere_radius = 5.0 * radius;
      double sm_dichotomy_error = sm_distance * spacing / 1000.0;
      
      Surface_3 surface(function,
                        Sphere (inner_point, sm_sphere_radius * sm_sphere_radius),
                        sm_dichotomy_error / sm_sphere_radius);

      CGAL::Surface_mesh_default_criteria_3<STr> criteria (sm_angle,
                                                          sm_radius * spacing,
                                                          sm_distance * spacing / 10);

      STr tr;
      C2t3 c2t3(tr);
      
      CGAL::make_surface_mesh(c2t3,
                              surface,
                              criteria,
                              CGAL::Manifold_with_boundary_tag());

      // saves reconstructed surface mesh
      std::string curr_outfile(std::to_string(i) + "_" + outfile);
      std::ofstream out(curr_outfile);
      Polyhedron output_mesh;
      CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, output_mesh);
      out << output_mesh;
    }
    
    std::string f_outfile(std::to_string(i) + "_fvalue.ply");
    function.draw_xslice_function(size, x, 0, f_outfile);

  } // for each input file
  
  std::cerr << std::endl;

  // Returns accumulated fatal error
  std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
  return accumulated_fatal_err;

  return 0;
}
