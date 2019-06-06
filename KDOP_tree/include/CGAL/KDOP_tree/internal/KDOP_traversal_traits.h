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

#ifndef CGAL_KDOP_TREE_INTERNAL_KDOP_TRAVERSAL_TRAITS_H_
#define CGAL_KDOP_TREE_INTERNAL_KDOP_TRAVERSAL_TRAITS_H_

#include <CGAL/KDOP_tree/internal/KDOP_node.h>

#include <boost/optional.hpp>

/// \file KDOP_traversal_traits.h

namespace CGAL {
namespace KDOP_tree {

/// \addtogroup PkgKDOPTree

  template<typename ValueType, typename IntegralType>
  class Counting_output_iterator
  {
    typedef Counting_output_iterator<ValueType, IntegralType> Self;
    IntegralType* i;

  public:
    Counting_output_iterator(IntegralType* i_) : i(i_) { };

    struct Proxy {
      Proxy& operator = (const ValueType&) { return *this; }
    };

    Proxy operator * () { return Proxy(); }

    Self& operator ++ () { ++*i; return *this; }

    Self& operator ++ (int) { ++*i; return *this; }

  };

  // the following are traits classes for traversal computation

  /*!
   * @class First_intersection_traits
   * Traits used in the k-dop tree traversal to get first intersection.
   * \todo Add member functions, for example computing intersection, check intersection, return intersection, etc.
   */
  template<typename KDOPTraits, typename Query>
  class First_intersection_traits
  {
    //TODO add member functions
  };

  /*!
   * @class Listing_intersection_traits
   * Traits used in the k-dop tree traversal to get intersections.
   * \todo Add member functions.
   */
  template<typename KDOPTraits, typename Query, typename OutputIterator>
  class Listing_intersection_traits
  {
    //TODO add member functions
  };

  /*!
   * @class Listing_primitive_traits
   * Traits used in the k-dop tree traversal to get intersected primitives.
   * \todo Add member functions.
   */
  template<typename KDOPTraits, typename Query, typename OutputIterator>
  class Listing_primitive_traits
  {
    //TODO add member functions
  };

  /*!
   * @class First_primitive_traits
   * Traits used in the k-dop tree traversal to get first intersected primitive.
   * \todo Add member functions.
   */
  template<typename KDOPTraits, typename Query>
  class First_primitive_traits
  {
    //TODO add member functions
  };

  /*!
   * @class Do_intersect_traits
   * Traits used in the k-dop tree traversal to check intersection.
   * \todo Add member functions.
   */
  template<typename KDOPTraits, typename Query>
  class Do_intersect_traits
  {
    //TODO add member functions
  };

  /*!
   * @class Compute_kdop_traits
   * Traits used in the k-dop tree traversal to compute k-dops.
   * \todo Add member functions.
   */
  template<typename KDOPTraits>
  class Compute_kdop_traits
  {

  };

  /// @}

} // namespace KDOP_tree
} // namespace CGAL


#endif // CGAL_KDOP_TREE_INTERNAL_KDOP_TRAVERSAL_TRAITS_H_
