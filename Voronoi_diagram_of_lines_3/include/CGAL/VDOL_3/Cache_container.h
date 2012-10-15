// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id:
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//
//

#ifndef CGAL_VDL3_CACHE_CONTAINER_H
#define CGAL_VDL3_CACHE_CONTAINER_H

#include <CGAL/VDOL_3/CGAL_SNAP_SVCET_3_TYPEDEFS.h>
#include <CGAL/tuple.h>

// The following is required only by MSVC, but harmless for other.
#include <boost/tuple/tuple_comparison.hpp>

namespace CGAL {
namespace VDOL_3 {

template <typename SVCET_3>
struct Cache_container{  
  CGAL_SNAP_SVCET_3_TYPEDEFS;

  /*! Cache for trisectors begin */
  typedef CGAL::cpp0x::tuple<Sign,Poly_int_3,Poly_int_3> Trisector_cache_key;
  typedef std::list<CGAL::Object>                        Trisector_cache_data;
  typedef std::map<Trisector_cache_key,Trisector_cache_data>  Trisector_cache;
  mutable Trisector_cache m_trisector_cache;
  /*! Cache for trisectors end */

//   // KEY: SVCET_3::State, ABOVE/BELOW, arc id, two lines 
//   typedef CGAL::cpp0x::tuple<Sign,Sign,int,Poly_int_3,Poly_int_3> Compare_at_arc_cache_key;
//   typedef boost::optional<std::pair<CGAL::Sign,X_monotone_curve_2> > Compare_at_arc_cache_data;
//   typedef std::map<Compare_at_arc_cache_key,Compare_at_arc_cache_data>  
//   Compare_at_arc_cache;
//   mutable Compare_at_arc_cache m_compare_at_arc_cache; 

  // projection cache
  typedef std::pair<Poly_int_3, Poly_int_3> Projection_cache_key; 
  typedef boost::optional<Poly_int_2> Projection_cache_data;  
  typedef std::map<Projection_cache_key, Projection_cache_data> 
  Projection_cache;
  mutable Projection_cache m_projection_cache;
  
  // projection cache
  typedef std::pair<Poly_int_3, Poly_int_3> Projection_sqff_cache_key; 
  typedef boost::optional<std::list<std::pair<Poly_int_2,int> > > Projection_sqff_cache_data; 
  typedef std::map<Projection_sqff_cache_key, Projection_sqff_cache_data>  
  Projection_sqff_cache;
  mutable Projection_sqff_cache m_projection_sqff_cache;

  // swap_2 cache 
  typedef Poly_int_2 Swap_2_cache_key;
  typedef boost::optional<Poly_int_2> Swap_2_cache_data;
  typedef std::map<Swap_2_cache_key, Swap_2_cache_data> Swap_2_cache;
  mutable Swap_2_cache m_swap_2_cache;

  // solve_1 cache
  typedef Poly_int_1 Solve_1_cache_key;
  typedef boost::optional<std::list<Coordinate_1> > Solve_1_cache_data;
  typedef std::map<Solve_1_cache_key, Solve_1_cache_data> Solve_1_cache;
  mutable Solve_1_cache m_solve_1_cache;
  

//   typedef std::pair<FT,FT> Rational_point_above_or_below_arc_cache_key;
//   typedef boost::optional<std::pair<FT,FT> > Rational_point_above_or_below_arc_cache_data;
//   typedef std::map<
//             Rational_point_above_or_below_arc_cache_key, 
//             Rational_point_above_or_below_arc_cache_data> 
//           Rational_point_above_or_below_arc_cache;
//   mutable Rational_point_above_or_below_arc_cache m_rational_point_above_or_below_arc_cache;
  
};

} // namespace VDOL_3
} //namespace CGAL

#endif // CGAL_VDL3_CACHE_CONTAINER_H
