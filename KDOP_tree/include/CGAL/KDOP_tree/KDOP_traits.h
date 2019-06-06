// Copyright (c) 2019  University of Cambridge (UK), GeometryFactory (France)
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s) : Xiao Xiao, Fehmi Cirak, Andreas Fabri

#ifndef CGAL_KDOP_TREE_KDOP_TRAITS_H_
#define CGAL_KDOP_TREE_KDOP_TRAITS_H_

#include <CGAL/Default.h>
#include <CGAL/KDOP_tree/KDOP_kdop.h>

namespace CGAL {
namespace KDOP_tree {

/// \addtogroup PkgKDOPTree
/// @{

  /**
   * Traits class
   * \todo Add KDOP_traits_base class.
   */

/// \tparam GeomTraits must  be a model of the concept \ref KDOPGeomTraits,
/// and provide the geometric types as well as the intersection tests and computations.
/// \tparam KdopPrimitive provide the type of primitives stored in the KDOP_tree.
///   It is a model of the concept `KDOPPrimitive`.
///
/// \tparam KdopMap must be a model of `ReadablePropertyMap` that has as key type a primitive id,
///                 and as value type a `Kdop`.
///                 If the type is `Default` the `Datum` must have the
///                 member function `kdop()` that returns the k-dop of the primitive.
///
/// If the argument `GeomTraits` is a model of the concept \ref
/// KdopRayIntersectionGeomTraits, this class is also a model of \ref
/// KdopRayIntersectionTraits.
///
/// \sa `KDOPTraits`
/// \sa `KDOP_tree`
/// \sa `KDOPPrimitive`

  template<typename GeomTraits, typename KdopPrimitive, typename KdopMap = Default>
  class KDOP_traits
  {
  public:

    /// \name Types
    /// @{

    /// Type of geometry traits (kernel)
    typedef GeomTraits Geom_traits;

    /// Type of fild number of the kernel
    typedef typename GeomTraits::FT FT;

    /// Type of k-dop traits
    typedef KDOP_traits<GeomTraits, KDOPPrimitive, KDOPMap> KT;

    /// Type of primitives
    typedef KDOPPrimitive Primitive;

    /// 3D point and Primitive Id type
    typedef typename std::pair<typename GeomTraits::Point_3, typename Primitive::Id> Point_and_primitive_id;

    /// Type of intersection result
    /// \todo definition of the structure
    template<typename Query>
    struct Intersection_and_primitive_id {
      //TODO definition of the structure
    };

    /// Type of 3D point
    typedef typename GeomTraits::Point_3 Point_3;

    /// Type of k-dop
    typedef typename CGAL::KDOP_tree::KDOP_kdop Kdop;

    KdopMap kdm;

    /// @}

    /// \name Constructor
    /// @{

    /// Default constructor
    KDOP_traits() { }

    /// Constructor with given k-dop map
    KDOP_traits(KdopMap kdm)
      : kdm(kdm)
    {}

    /// @}

    /**
     * Split a range of primitives defined by [first, beyond).
     *
     * @param first iterator on the first element
     * @param beyond iterator on the past-the-end element
     *
     * \todo Split the primitives with an octree or a binary tree.
     *
     */
    class Split_primitives
    {
      //TODO split the primitives recursively.
    };

    Split_primitives split_primitives_object() const {return Split_primitives(*this);}

    /**
     * Compute the k-dop of a set of primitives.
     *
     * @param first iterator on the first primitive
     * @param beyond iterator on the past-the-end primitive
     *
     * @return the k-dop of the primitives within the iterator range
     *
     * \todo Recursively compute the kdops of nodes in the tree, including the
     * union operation to obtain the k-dop of a node from its children.
     *
     */
    class Compute_kdop
    {
      //TODO compute support heights, define operator
    };

    Compute_kdop compute_kdop_object() const {return Compute_kdop(*this);}

    /**
     * Check if the query intersects the primitive
     *
     * @param q query object
     * @param pr primitive
     *
     * @return a bool result
     *
     * \todo Define operators to check intersection with k-dops.
     *
     */
    class Do_intersect
    {
      //TODO define operators to check intersection
    };

    Do_intersect do_intersect_object() const {return Do_intersect(*this);}

    /**
     * Compute the intersection between a query and a primitive
     *
     * @param query query object
     * @param primitive primitive
     *
     * @return the intersection result
     *
     * \todo Define operators to compute intersection.
     *
     */
    class Intersection
    {
      //TODO define operator to compute intersection
    };

    Intersection intersection_object() const {return Intersection(*this);}

  private:
    /**
     * Compute the k-dop of a primitive
     *
     * @param pr primitive
     *
     * @return the k-dop of the primitive \c pr
     *
     */
    template<typename PM>
    Kdop compute_kdop(const Primitive& pr, const PM&) const
    {
      return get(kdm, pr.id());
    }

  };

/// @}

}
}



#endif // CGAL_KDOP_TREE_KDOP_TRAITS_H_
