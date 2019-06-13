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

#ifndef CGAL_KDOP_TREE_INTERNAL_KDOP_RAY_INTERSECTION_H_
#define CGAL_KDOP_TREE_INTERNAL_KDOP_RAY_INTERSECTION_H_

#include <functional>
#include <boost/optional.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/variant/apply_visitor.hpp>

#if BOOST_VERSION >= 105000
#  if defined(BOOST_MSVC)
#    pragma warning(push)
#    pragma warning(disable: 4996)
#  endif
#  include <boost/heap/priority_queue.hpp>
#  if defined(BOOST_MSVC)
#    pragma warning(pop)
#  endif
#else
#  include <queue>
#endif

#include <CGAL/assertions.h>

/// \file KDOP_ray_intersection.h

namespace CGAL {
namespace KDOP_tree {
namespace internal {

/// \addtogroup PkgKDOPTree
/// @{

/*!
 * Class to compute intersection with rays.
 *
 */

  template<typename KDOPTree, typename SkipFunctor>
  class KDOP_ray_intersection
  {
    typedef typename KDOPTree::KDOP_traits KDOP_traits;
    typedef typename KDOP_traits::Ray_3 Ray;
    typedef typename KDOPTree::template Intersection_and_primitive_id<Ray>::Type Ray_intersection_and_primitive_id;
    typedef typename Ray_intersection_and_primitive_id::first_type Ray_intersection;

  public:
    /// \name Constructor
    /// @{

    /// Constructor
    /// \param tree the k-dop tree object
    KDOP_ray_intersection(const KDOPTree& tree) : tree_(tree) { }

    /// @}

    /// \name Functions
    /// @{

    /// Return intersection points between the ray and the k-dop tree
    /// \param query the query
    /// \param skip skip function
    /// \return intersection points
    /// \todo Add code based on the splitting method.
    boost::optional< Ray_intersection_and_primitive_id >
    ray_intersection(const Ray& query, SkipFunctor skip) const {

      //TODO function body

      boost::optional< Ray_intersection_and_primitive_id > intersection, p;

      return p;
    }

    /// @}

  private:
    const KDOPTree& tree_;
    typedef typename KDOPTree::Point Point;
    typedef typename KDOPTree::FT FT;
    typedef typename KDOPTree::Node Node;
    typedef typename KDOPTree::size_type size_type;

  }; // end KDOP_ray_intersection class

  //TODO "first_intersection" function

  //TODO "first_intersected_primitive" function

  /// @}

} // namespace internal
} // namespace KDOP_tree
} // namespace CGAL



#endif // CGAL_KDOP_TREE_INTERNAL_KDOP_RAY_INTERSECTION_H_
