// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_SMS_3_SINGLE_MOLD_TRANSLATIONAL_CASTING_TOP_FACETS_H
#define CGAL_SMS_3_SINGLE_MOLD_TRANSLATIONAL_CASTING_TOP_Facets_H

#include <iostream>
#include <list>
#include <vector>

#include <boost/type_traits/is_same.hpp>

#include <CGAL/Kernel_traits.h>
#include <CGAL/enum.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/property_map.h>

#include <CGAL/Set_movable_separability_3/NormalsDirection.h>
#include <CGAL/Set_movable_separability_3/check_upper_facet.h>
#include <CGAL/Set_movable_separability_3/coveringset_finder.h>
#include <CGAL/Set_movable_separability_3/lp_wrapper.h>

namespace CGAL {

enum geom_traits_t { geom_traits };
enum all_default_t { all_default };

namespace Set_movable_separability_3 {
namespace Single_mold_translational_casting {

template <typename Polyhedron, typename NamedParameters>
class Get_kernel {
  typedef typename boost::property_map<Polyhedron,
                                       boost::vertex_point_t>::const_type
    Propert_map;
  typedef typename boost::property_traits<Propert_map>::value_type      Point;

  typedef typename CGAL::Kernel_traits<Point>::Kernel
    Default_kernel;

public:
  typedef typename CGAL::internal_np::Lookup_named_param_def<CGAL::geom_traits_t,
                                                       NamedParameters,
                                                       Default_kernel>::type
    type;
};

template <typename Polyhedron,  typename NamedParameters,
          typename OutputIterator>
OutputIterator top_facets_impl(const Polyhedron& polyhedron,
                               const NamedParameters& np,
                               OutputIterator oi, boost::false_type)
{
  return oi;
}

//return single direction for each facet
template <typename Polyhedron,  typename NamedParameters,
          typename OutputIterator>
OutputIterator top_facets_impl(const Polyhedron& polyhedron,
                               const NamedParameters& np,
                               OutputIterator oi, boost::true_type)
{
  typedef typename Get_kernel<Polyhedron, NamedParameters>::type    Kernel;

  typedef typename Kernel::Direction_3                         	    Direction_3;

  std::vector< std::pair<typename Kernel::Direction_3, std::vector<typename Polyhedron::Facet_const_iterator>>>
  polyhedronNormals =  internal::findDirections<Polyhedron,Kernel>(polyhedron);

  int outLength(0);
   unsigned int outIndexs[6];

   // 2 hemisphere that the intersection of their boundary is a point that
   // wasn't covered
   Direction_3 outDirection;
   bool outDirectionExists;

   outLength =
     internal::findCoveringSet<Kernel>(polyhedronNormals, outIndexs,
                                       &outDirection,&outDirectionExists);
   while (outLength--) {
     std::pair<bool, Direction_3> tmp =
       internal::checkUpperFacet<Kernel>(polyhedronNormals,
                                         outIndexs[outLength]);
     if(tmp.first) {
       *oi++ = std::make_pair(polyhedronNormals[outIndexs[outLength]].second,
                              tmp.second);
     }
   }
   if (outDirectionExists) {
     std::pair<bool, unsigned int> topFacet =
       internal::checkDirection<Kernel>(polyhedronNormals, outDirection);
     if (topFacet.first) {
       *oi++ = std::make_pair(polyhedronNormals[topFacet.second].second,
                              outDirection);

     }
     topFacet =
       internal::checkDirection<Kernel>(polyhedronNormals, -outDirection);
     if(topFacet.first) {
       *oi++ = std::make_pair(polyhedronNormals[topFacet.second].second,
                              -outDirection);
     }
   }

  return oi;
}

/*! \fn OutputIterator find_top_facets(const Polyhedron& polyhedron, OutputIterator oi)
 * \param[in] polyhedron the input polyhedron.
 * \param[out] oi the output iterator. Its value type is a pair, where
 *             (i) the first element in the pair identifies a valid top face
 *                 represented by its index the type of which is convertible to
                   `boost::graph_traits<Polyhedron>::face_descriptor`, and
 *             (ii) the second element is a closed spherical patch of pull-out
 *                  3D directions represented as a sequence of the extreme
 *                  directions in the patch of type `Kernel::Direction_3`.
 * \return the past-the-end iterator of the output container.
 * \pre `polyhedron` must be non-degenerate (has at least 4 vertices and 6
 *      edges), simple, and does not have neighboring coplanar facets.
 */
template <typename Polyhedron,  typename NamedParameters,
          typename OutputIterator, typename DirectionType>
OutputIterator top_facets(const Polyhedron& polyhedron,
                          const NamedParameters& np, OutputIterator oi)
{
  typedef typename Get_kernel<Polyhedron, NamedParameters>::type Kernel;
  typedef typename Kernel::Direction_3                           Direction_3;

  //! \todo consider using CGAL::is_same_or_derived instaed of boost::is_same
  typedef typename boost::is_same<DirectionType, Direction_3>::type
    Is_direction;

  return top_facets_impl(polyhedron, np, oi, Is_direction());
}

template <typename Polyhedron,  typename NamedParameters,
          typename OutputIterator>
OutputIterator top_facets(const Polyhedron& polyhedron,
                          const NamedParameters& np, OutputIterator oi)
{
  typedef typename value_type_traits<OutputIterator>::type       Value_type;
  typedef typename Value_type::second_type                       Direction_type;

  typedef typename Get_kernel<Polyhedron, NamedParameters>::type Kernel;
  typedef typename Kernel::Direction_3                           Direction_3;

  //! \todo consider using CGAL::is_same_or_derived instead of boost::is_same
  typedef typename boost::is_same<Direction_type, Direction_3>::type
    Is_direction;

  return top_facets_impl(polyhedron, np, oi, Is_direction());
}

} // end of namespace Single_mold_translational_casting
} // end of namespace Set_movable_separability_3
} // end of namespace CGAL

#endif
