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

namespace CGAL {
namespace KDOP_tree {

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
  /**
   * @class First_intersection_traits
   */
  template<typename KDOPTraits, typename Query>
  class First_intersection_traits
  {
    //TODO add member functions
  };

  template<typename KDOPTraits, typename Query, typename OutputIterator>
  class Listing_intersection_traits
  {
    //TODO add member functions
  };

  template<typename KDOPTraits, typename Query, typename OutputIterator>
  class Listing_primitive_traits
  {
    //TODO add member functions
  };

  template<typename KDOPTraits, typename Query>
  class First_primitive_traits
  {
    //TODO add member functions
  };

  template<typename KDOPTraits, typename Query>
  class Do_intersect_traits
  {
    //TODO add member functions
  };

  template<typename KDOPTraits>
  class Projection_traits
  {
    //TODO add member functions
  };

}
}


#endif // CGAL_KDOP_TREE_INTERNAL_KDOP_TRAVERSAL_TRAITS_H_
