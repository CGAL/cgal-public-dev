// rename isr_test_reconstruction_utils.h
#ifndef ISR_TEST_UTIL_RECONSTRUCTION_H
#define ISR_TEST_UTIL_RECONSTRUCTION_H

// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include <iostream>

#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <Eigen/Core>

//Mesh
#include <CGAL/Surface_mesh.h>

//Mesher
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polyhedron_3.h>

//Reconstruction
#include <CGAL/Implicit_reconstruction_function.h>

//PMP
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

//Boost
#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>

//file includes & utils
#include "isr_test_types.h"
#include "isr_test_param_class.h"
#include "isr_test_io_utils.h"

namespace PMP = CGAL::Polygon_mesh_processing;

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

//property maps
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// Spectral implicit function
typedef CGAL::Implicit_reconstruction_function<Kernel, PwnList, Normal_map> Implicit_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Implicit_reconstruction_function> Surface_3;


// ----------------------------------------------------------------------------


bool surface_mesh_reconstruction(const Param &p, PwnList &pwnl, Mesh &m) // surface_mesh_reconstruction
{
  static int i = 0;
  ++i;

  //COMPUTES IMPLICIT FUNCTION
  Implicit_reconstruction_function function;
  function.initialize_point_map(pwnl, Point_map(), Normal_map(), p.octree, 0);

  if(p.poisson)
  {
    if (! function.compute_poisson_implicit_function()){
      std::cerr << "Error: cannot compute implicit function" << std::endl;
      return false;
    }
  }
  else if (p.spectral)
  {
    double fitting = 0.1;
    double ratio = 10.0;
    double bilaplacian = 1.0;
    double laplacian = 0.1; /*default values*/
    if (! function.compute_spectral_implicit_function(fitting, ratio, bilaplacian, laplacian) )
    {
      std::cerr << "Error: cannot compute implicit function" << std::endl;
      return false;
    }
  }

  //BUILDS RECONSTRUCTION
  if(p.march_tets)
  {
    double isovalue = 0.0;
    if(!function.marching_tetrahedra(isovalue, m)) {
      std::cerr << "Error : Mesh is not 2-manifold" << std::endl;
            std::string curr_outfile(std::to_string(i) + ".off");
            std::ofstream out(curr_outfile);
            out << m;
      return false;
    }
  }
  else if (p.make_sm)
  {
    Point inner_point = function.get_inner_point();
    Sphere bsphere = function.bounding_sphere();
    double radius = std::sqrt(bsphere.squared_radius());

    // Computes average spacing
    double spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(pwnl, 6, CGAL::parameters::point_map(Point_map()) /* knn = 1 ring */);
    double sm_sphere_radius = 5.0 * radius;
    double sm_distance = 0.25;
    double sm_dichotomy_error = sm_distance * spacing / 1000.0;
    double sm_angle = 20.0;
    double sm_radius = 100.0;
    
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
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, m);
  }
  
  if (m.is_empty()) {
    std::cout << "Error : Mesh is empty" << std::endl;
    return false;
  }
  
  if (!m.is_valid()) {
    std::cout << "Error : Mesh is not valid" << std::endl;
  }

      std::string curr_outfile(std::to_string(i) + ".off");
      std::ofstream out(curr_outfile);
      out << m;

  return true;
}

#endif //ISR_TEST_UTIL_RECONSTRUCTION_H
