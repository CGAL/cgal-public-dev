// Copyright (c) 2005-2009 Tel-Aviv University (Israel).
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

// SVCET_3 = Single_voronoi_cell_envelope_traits_3


#ifndef CGAL_VDL3_SVCET_33_CONSTRUCT_PROJECTED_BOUNDARY_2
#define CGAL_VDL3_SVCET_33_CONSTRUCT_PROJECTED_BOUNDARY_2

#include <CGAL/config.h>
#include <CGAL/VDOL_3/SVCET_3_functors/SVCET_3_functor_base.h>


namespace CGAL { 
namespace VDOL_3 {
namespace SVCET_3_functors {

template < class SVCET_3 >
class Construct_projected_boundary_2 : public 
Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > {

public:
  //! first template parameter
  typedef SVCET_3 Single_voronoi_cell_envelope_traits_3;
  //! the base type 
  typedef Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > Base;
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
  Construct_projected_boundary_2(const SVCET_3 *traits): Base(traits) {};
  
public:
  // insert into the OutputIterator all the (2d) curves of the boundary of
  // the vertical projection of the surface on the xy-plane
  // the OutputIterator value type is X_monotone_curve_2
  template <class OutputIterator>
  OutputIterator
  operator()(const Xy_monotone_surface_3& s, OutputIterator o)
  {
    // Just take the boundaries from the surface and make them CGAL::Object

    typedef typename Xy_monotone_surface_3::Boundary_container
                                                            Boundary_container;
    typedef typename Boundary_container::const_iterator     Boundary_const_iter;

    for (Boundary_const_iter it = s.boundary().begin();
         it != s.boundary().end(); ++it)
      *o++ = make_object(*it);

    return o;
  }  
};

} // namespace SVCET_3_functors
} // namespace VDOL_3 
} // namespace CGAL

#endif // CGAL_VDL3_SVCET_33_CONSTRUCT_PROJECTED_BOUNDARY_2
