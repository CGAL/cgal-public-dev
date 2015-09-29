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
// Author(s): Ophir Setter          <ophir.setter@cs.tau.ac.il>
//
//

#ifndef CGAL_ALGEBRAIC_FARTHEST_POINT_FARTHEST_SITE_TRAITS_2_H
#define CGAL_ALGEBRAIC_FARTHEST_POINT_FARTHEST_SITE_TRAITS_2_H

/*! \file
 * This is the file containing a traits class to create Voronoi diagrams
 * with respect to the furthest apollonius distance function. This distance
 * function is used to calculate the minimum width annulus of disks (see 
 * below.)
 */

#include <CGAL/basic.h>

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Envelope_voronoi_traits_2/Algebraic_apollonius_traits_base_2.h>

namespace CGAL {

/*!
 * \class Algebraic_farthest_point_farthest_site_traits_2 We define a distance function
 * from a disk to as the longest distance from the given point to a point on
 * on the disk. We define the Vornoi diagram with respect to that distance 
 * funtion - "furthest" distance function. 
 * In short, the distance that this traits class supports is:
 *        $d(x) = \sqrt{(p_x - c_x)^2 + (p_y - c_y)^2} + r
 *
 * Where x = (p_x, p_y), (c_x, c-y) is the center of the circle and r is its
 * radius.
 *        
 * We can use this Voronoi diagram together with the regular apollonius diagram
 * to compute the minimum width annulus that contains the circles.
 */  

template <class CurvedKernel_>
class Algebraic_farthest_point_farthest_site_traits_2 
: public Algebraic_apollonius_traits_base_2< 
                                             CurvedKernel_, 
                                             Algebraic_farthest_point_farthest_site_traits_2
                                             <
                                               CurvedKernel_
                                             >,
                                             Envelope_voronoi_traits_2::
                                              _Apollonius_disk< 
                                               typename CurvedKernel_::
                                               Curve_kernel_2::
                                               Coordinate_1::
                                               Rational 
                                               >,
                                             CGAL::Field_with_sqrt_tag
                                           >
{
  typedef CurvedKernel_                            Curved_kernel_2;
  typedef Algebraic_farthest_point_farthest_site_traits_2< 
    Curved_kernel_2 >                              Self;
  typedef Envelope_voronoi_traits_2::_Apollonius_disk< 
    typename CurvedKernel_::
    Curve_kernel_2::
    Coordinate_1::
    Rational >                                     Site_2;

  
  typedef Algebraic_apollonius_traits_base_2< 
    Curved_kernel_2, Self, Site_2, 
    CGAL::Field_with_sqrt_tag >                    Base;

public:

  typedef typename Base::Xy_monotone_surface_3     Xy_monotone_surface_3;


  //! Compute the distance from a point to a site.
  /*! Compute the distance from a point to a site.
    \param px The x-value of the point.
    \param py The y-value of the point.
    \param site The site we compute the distance to.
    \return The distance from the point to this site.
  */
  template <typename NT>
    static NT distance(const NT &px,
                       const NT &py, const Site_2& site)
  {
    const NT cx = site.center().first;
    const NT cy = site.center().second;
    
    NT x = px - cx;
    NT y = py - cy;
    
    return sqrt(x*x + y*y) + NT(site.r());
  }

  class Construct_projected_intersections_2
  {
  protected:
    const Self * m_traits;

  public:
    Construct_projected_intersections_2(const Self * traits)
      : m_traits(traits)
    {}

    template <class OutputIterator>
      OutputIterator operator()(const Xy_monotone_surface_3& s1,
                                const Xy_monotone_surface_3& s2,
                                OutputIterator o) const
    {
      // The bisector of this traits function is similar to the
      // bisector of the apollonius diagram:
      //
      //        ||x - p_1|| + r_1 = ||x - p_2|| + r_2
      //
      //                    ||
      //                    \/
      //
      //        ||x - p_1|| - r_2 = ||x - p_2||  r_1
      //
      // So we just need to swap the radii.
      return m_traits->apollonius_bisector(s1.center().first,
                                           s1.center().second,
                                           s2.r(),
                                           s2.center().first,
                                           s2.center().second,
                                           s1.r(),
                                           *m_traits,
                                           o);
    }
  };

  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const
  {
    return Construct_projected_intersections_2(this);
  }
};


} //namespace CGAL

#endif //  CGAL_ALGEBRAIC_FARTHEST_POINT_FARTHEST_SITE_TRAITS_2_H
