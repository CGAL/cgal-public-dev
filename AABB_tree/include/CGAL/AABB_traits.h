// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
//
// Author(s) : Stéphane Tayeb, Pierre Alliez, Camille Wormser
//

#ifndef CGAL_AABB_TRAITS_H_
#define CGAL_AABB_TRAITS_H_

#include <CGAL/AABB_traits_d.h>
/// \file AABB_traits.h

namespace CGAL {

/// \addtogroup PkgAABB_tree
/// @{

/// This traits class handles any type of 2D/3D geometric
/// primitives provided that the proper intersection tests and
/// constructions are implemented. It handles points, rays, lines and
/// segments as query types for intersection detection and
/// computations, and it handles points as query type for distance
/// queries. The geometry specific functions are derived from the
/// base class AABB_traits_d.
/// \cgalModels AABBTraits
/// \tparam GeomTraits must  be a model of the concept `AABBGeomTraits_2` or `AABBGeomTraits_3`,
/// and provide the geometric types as well as the intersection tests and computations.
/// \tparam Primitive provide the type of primitives stored in the `AABB_tree`.
///   It is a model of the concept `AABBPrimitive` or `AABBPrimitiveWithSharedData`.
///
/// \sa `AABBTraits`
/// \sa `AABB_tree`
/// \sa `AABBPrimitive`
/// \sa `AABBPrimitiveWithSharedData`
template<typename GeomTraits, typename Primitive>
class AABB_traits:
  public AABB_traits_d<GeomTraits,Primitive,Primitive::Datum::Ambient_dimension::value>
{
public:
  /// Default constructor.
  AABB_traits() { };

   /*
   * Computes the bounding box of a set of primitives
   * @param first an iterator on the first primitive
   * @param beyond an iterator on the past-the-end primitive
   * @return the bounding box of the primitives of the iterator range
   */
  class Compute_bbox {
    const AT& m_traits;
  public:
    Compute_bbox(const AT& traits)
      :m_traits (traits) {}

    template<typename ConstPrimitiveIterator>
    typename AT::Bounding_box operator()(ConstPrimitiveIterator first,
                                         ConstPrimitiveIterator beyond) const
      {
        typename AT::Bounding_box bbox = internal::Primitive_helper<AT>::get_datum(*first,m_traits).bbox();
        for(++first; first != beyond; ++first)
        {
          bbox = bbox + internal::Primitive_helper<AT>::get_datum(*first,m_traits).bbox();
        }
        return bbox;
      }
  };

  Compute_bbox compute_bbox_object() const {return Compute_bbox(*this);}


};  // end class AABB_traits


}  // end namespace CGAL

#endif // CGAL_AABB_TRAITS_H_
