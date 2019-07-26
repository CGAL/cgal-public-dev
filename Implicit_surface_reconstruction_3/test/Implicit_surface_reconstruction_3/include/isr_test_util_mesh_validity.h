#ifndef ISR_TEST_UTIL_MESH_VALIDITY_H
#define ISR_TEST_UTIL_MESH_VALIDITY_H

//includes
#include "include/isr_test_types.h"

#include <boost/foreach.hpp>

#include <CGAL/Polygon_mesh_processing/repair.h>

namespace PMP = CGAL::Polygon_mesh_processing;

//types


bool is_valid(const Mesh &m) /*changer la partie manifold une fois que Tong a fini sa partie*/
{
  bool validity_res = true;
  bool valid = m.is_valid();
  bool empty = m.is_empty();
  bool manifold = true;

/*  BOOST_FOREACH(boost::graph_traits<Mesh>::vertex_descriptor vd, m.vertices()) {
    if (PMP::is_non_manifold_vertex(vd,m)) {
      manifold = false;
      break;
    }
  }*/ 

  if (empty) {
    std::cout << "Error : reconstructed mesh is empty" << std::endl;
    validity_res = false;
  }

/*  if (!manifold) {
    std::cout << "Error : reconstructed mesh is not 2-manifold" << std::endl;
    validity_res = false;
  }*/

  if (!valid) {
    std::cout << "Error : reconstructed mesh is not valid" << std::endl;
    validity_res = false;
  }

  return (validity_res);
}



#endif //ISR_TEST_UTIL_MESH_VALIDITY_H