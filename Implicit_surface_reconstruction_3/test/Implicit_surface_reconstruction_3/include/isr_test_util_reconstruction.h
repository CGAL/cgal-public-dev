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
#include "isr_test_util_file_reading.h"
#include "isr_test_util_mesh_validity.h"

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

//boost
typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;


// ----------------------------------------------------------------------------


bool read_file(const std::string &in_file, PwnList &pwnl, Mesh &m)
{
  bool success = true;

  //READS INPUT FILE
  if(is_mesh(in_file))
  {
    if (read_input_mesh_file(in_file, m)) 
    {
      Mesh::Property_map<vertex_descriptor, Vector> vnormals_pm = m.add_property_map<vertex_descriptor, Vector>
                                                              ("v:normals", CGAL::NULL_VECTOR).first;
      BOOST_FOREACH(vertex_descriptor v, m.vertices()) {
        const Point& p = m.point(v);
        Vector n = PMP::compute_vertex_normal(v , m , vnormals_pm);
        pwnl.push_back(std::make_pair(p, n));
      }
    }
    else
      success = false;
  }
  else if (is_point_set(in_file))
  {
    if (! read_input_point_set_file(in_file, pwnl)) {
      success = false;
    }
  }
  else
  {
    std::cerr << "Error: file not supported" << in_file << std::endl;
    success = false;
  }

  return success;
}


bool mesh_reconstruction(const std::string &in_file, const Param &p, PwnList &pwnl, Mesh &m)
{
  bool success = true;

  //READS INPUT FILE
  if(!read_file(in_file, pwnl, m)) {
    success = false;
    return (success);
  }

  //CHECK REQUIREMENTS
  std::size_t nb_points = pwnl.size();

  if (nb_points == 0)
  {
    std::cerr << "Error: empty point set" << std::endl;
    success = false;
    return (success);
  }

  bool points_have_normals = (pwnl.begin()->second != CGAL::NULL_VECTOR);
  if ( ! points_have_normals )
  {
    std::cerr << "Input point set not supported: this reconstruction method requires unoriented normals" << std::endl;
    // this is not a bug => do not set success
    return true;
  }

  //COMPUTES IMPLICIT FUNCTION
  Implicit_reconstruction_function function;
  if (p.octree)
    function.initialize_point_map(pwnl, Point_map(), Normal_map(), 1, 0); /*dernier argument = octree debug peut etre ; mettre variable?*/
  else if (p.del_ref)
    function.initialize_point_map(pwnl, Point_map(), Normal_map(), 0, 0);

  if(p.poisson)
  {
    if (! function.compute_poisson_implicit_function()){
      std::cerr << "Error: cannot compute implicit function" << std::endl;
      success = false;
      return (success);
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
      success = false;
      return (success);
    }
  }

  //BUILDS RECONSTRUCTION
  if(p.march_tets)
  {
    double isovalue = 0.0;
    function.marching_tetrahedra(isovalue, m);
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
  
  if (! is_valid(m))
    success = false;

  return (success);
}

#endif //ISR_TEST_UTIL_RECONSTRUCTION_H