// Copyright (c) 2010  Tel-Aviv University (Israel).
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
// $Id: $
// 
//
// Author(s)     : Asaf Porat          <asafpor1@post.tau.ac.il>

#ifndef LINES_THROUGH_SEGMENTS_TRAITS_3_H
#define LINES_THROUGH_SEGMENTS_TRAITS_3_H

#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
/*! \file
*************************************************************
* The following class represents traits for the computation of line through
* segments.
*
* Input:
*
* 1. Algebraic number type -  The returned lines are algebraic.
* 2. Rational number type - The input segments are rational.
* 3. Traits class for a 2D arrangement on plane. 
* 4. Traits class for a 2D arrangement on sphere. 
*
*************************************************************
*/

namespace CGAL {

template <typename Alg_kernel_, 
          typename Rational_kernel_,
          typename Traits_arr_on_plane_2_ =  
          CGAL::Arr_conic_traits_2<Rational_kernel_, 
                                   Alg_kernel_, 
                                   CGAL::CORE_algebraic_number_traits>,
          typename Traits_arr_on_sphere_2_ = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Rational_kernel_> >
class Lines_through_segments_traits_3
{
public:
  typedef Alg_kernel_                                  Alg_kernel;
  typedef Rational_kernel_                             Rational_kernel;
  typedef typename Rational_kernel::Segment_3          Rational_segment_3;
  typedef CGAL::Arr_consolidated_curve_data_traits_2<
    Traits_arr_on_plane_2_,const Rational_segment_3*>  Traits_arr_on_plane_2;
  typedef CGAL::Arr_consolidated_curve_data_traits_2<
    Traits_arr_on_sphere_2_,const Rational_segment_3*> Traits_arr_on_sphere_2;
};
 
} //namespace CGAL

#endif /*LINES_THROUGH_SEGMENTS_TRAITS_3_H*/
