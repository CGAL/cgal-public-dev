// rename file isr_test_bbox_utils.h
#ifndef ISR_TEST_UTIL_BBOX_H
#define ISR_TEST_UTIL_BBOX_H

// ----------------------------------------------------------------------------
// Includes
// ----------------------------------------------------------------------------

#include "isr_test_types.h"

//Mesh
#include <CGAL/Surface_mesh.h> // maybe dont need

//Bounding box
#include <CGAL/bounding_box.h>

//boost
#include <boost/function_output_iterator.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

double util_bb_diag(PwnList pwnl) // bb_diag_length
{
  typedef typename PwnList::value_type PwnList_t;
  boost::function<Point(PwnList_t&)> pwn_it_to_point_it = boost::bind(&PwnList_t::first, _1);
  Kernel::Iso_cuboid_3 c3 = CGAL::bounding_box(boost::make_transform_iterator(pwnl.begin(), pwn_it_to_point_it), 
                                               boost::make_transform_iterator(pwnl.end(), pwn_it_to_point_it));

  return CGAL::sqrt(CGAL::squared_distance(c3.min(), c3.max()));
}

#endif //ISR_TEST_UTIL_BBOX_H
