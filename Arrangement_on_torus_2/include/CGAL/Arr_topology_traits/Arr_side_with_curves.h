// Copyright (c) 2009 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>

#ifndef CGAL_ARR_SIDE_WITH_CURVES_H
#define CGAL_ARR_SIDE_WITH_CURVES_H 1

/*!\file include/CGAL/Arr_topology_traits_2/Arr_side_with_curves.h
 * \brief Class to collect all specific operations for a closed side
 * or pair of identified sides.
 */

#include <CGAL/config.h>

namespace CGAL {


/*!\brief Class template that provides operations
 * related to a closed side or a pair of opposite identified sides
 */
template < class TopologyTraits_2, 
           Arr_parameter_space pss,
           Arr_parameter_space psl >
class Arr_side_with_curves {
  
public:
  
#if 0 // TODO add mpl assert
  BOOST_MPL_ASSERT(pss == psl || 
                   pss == ARR_LEFT_BOUNDARY && psl == ARR_RIGHT_BOUNDARY ||
                   pss == ARR_BOTTOM_BOUNDARY && psl == ARR_TOP_BOUNDARY);
#endif  

  //! this instance template parameter
  typedef TopologyTraits_2 Topology_traits_2;

  ///! \name The geometry-traits types.
  //!@{
  
  //! the adapted geometric traits
  typedef typename Topology_traits_2::Traits_adaptor_2 Traits_adaptor_2;

  //! point type
  typedef typename Traits_adaptor_2::Point_2 Point_2;

  //! curve type
  typedef typename Traits_adaptor_2::X_monotone_curve_2 X_monotone_curve_2;
  
  //!@}

  //! \name The DCEL types.
  //!@{
  typedef typename Topology_traits_2::Dcel           Dcel;
  typedef typename Dcel::Size                       Size;
  typedef typename Dcel::Vertex                     Vertex;
  typedef typename Dcel::Halfedge                   Halfedge;
  typedef typename Dcel::Face                       Face;
  typedef typename Dcel::Outer_ccb                  Outer_ccb;
  typedef typename Face::Outer_ccb_const_iterator   Outer_ccb_const_iterator;
  typedef typename Dcel::Inner_ccb                  Inner_ccb;
  typedef typename Dcel::Isolated_vertex            Isolated_vertex;
  
  //}@
  
protected:
  
  struct Less_point_2 {
    
    /*! Construct default */
    Less_point_2() : 
      _m_traits(NULL) {
    }
    
    /*! Construct */    
    Less_point_2(Traits_adaptor_2 * traits) : 
      _m_traits(traits) {
    }
    
    Traits_adaptor_2 * _m_traits;
    
    bool operator()(const Point_2& p1, const Point_2& p2) const {
      if (pss == ARR_LEFT_BOUNDARY || pss == ARR_RIGHT_BOUNDARY) {
        // compare on left and/or right side
        return (_m_traits->compare_y_on_boundary_2_object()(
                    p1, p2
                ) == CGAL::SMALLER);
      }
      // else
      assert(pss == ARR_BOTTOM_BOUNDARY || pss == ARR_TOP_BOUNDARY);
      // compare on bottom and or top side
      return (_m_traits->compare_x_on_boundary_2_object()(
                  p1, p2
              ) == CGAL::SMALLER);
    }
  };
  
  //! type of curve of identification
  typedef std::map< Point_2, Vertex*, Less_point_2 > Points_on_side;
  
  // check Vertex_less
  struct Less_vertex {
    bool operator() (Vertex *v1, Vertex *v2) {
      return &(*v1) < &(*v2);
    }
  };
  
  typedef std::map< Vertex*, typename Points_on_side::iterator, Less_vertex >
  Vertices_on_side;
  
  
public:
  
  //!\name Constructors
  //!@{
  
  Arr_side_with_curves(Traits_adaptor_2 *traits) : 
    _m_traits(traits),
    _m_points_on_side(Points_on_side(Less_point_2(_m_traits))) {
  }
  
  //!@}
  
  // TODO place boundary vertex
  
  // TODO locate ...
  
protected:
  
  //! map to associate points with vertices
  Points_on_side _m_points_on_side;
  
  //! map to associate vertices with position (iterator) in previous map
  Vertices_on_side _m_vertices_on_sides;
  
  Traits_adaptor_2 *_m_traits;
  
};
  
} //namespace CGAL

#endif // CGAL_ARR_DUPIN_CYCLIDE_TOPOLOGY_TRAITS_2_IMPL_H
// EOF
