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
   *
   */

//TODO add KDOP_traits_base class

  template<typename GeomTraits, typename KDOPPrimitive, typename KDOPMap = Default>
  class KDOP_traits
  {
  public:
    typedef GeomTraits Geom_traits;

    typedef KDOP_traits<GeomTraits, KDOPPrimitive, KDOPMap> KT;

    typedef typename GeomTraits::FT FT;

    typedef KDOPPrimitive Primitive;

    typedef typename std::pair<typename GeomTraits::Point_3, typename Primitive::Id> Point_and_primitive_id;

    //TODO define Intersection_and_primitive_id struct

    typedef typename GeomTraits::Point_3 Point_3;

    typedef typename CGAL::KDOP_tree::KDOP_kdop Kdop;

    KDOPMap kdm;

    // Default constructor
    KDOP_traits() { }

    KDOP_traits(KDOPMap kdm)
      : kdm(kdm)
    {}

    /**
     * Sorts the range defined by [first, beyond).
     *
     * @param first iterator on the first element
     * @param beyond iterator on the past-the-end element
     *
     */
    class Split_primitives
    {
      //TODO split primitives using bbox?
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
