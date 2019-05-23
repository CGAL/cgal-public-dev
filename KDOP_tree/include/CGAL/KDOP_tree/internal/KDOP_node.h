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

#ifndef CGAL_KDOP_TREE_INTERNAL_KDOP_NODE_H_
#define CGAL_KDOP_TREE_INTERNAL_KDOP_NODE_H_

#include <CGAL/Profile_counter.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/intersections.h>

#include <vector>

namespace CGAL {
namespace KDOP_tree {

  template<typename KDOPTraits>
  class KDOP_node
  {
  public:
    typedef typename KDOPTraits::Kdop Kdop;

    // constructor
    KDOP_node()
      : m_kdop()
      , m_p_left_child(NULL)
      , m_p_right_child(NULL) { };

    // Non-virtual destructor
    ~KDOP_node() { };

    // return the kdop of the node
    const Kdop& kdop() const { return m_kdop; }

    template<typename ConstPrimitiveIterator>
    void expand(ConstPrimitiveIterator first,
                ConstPrimitiveIterator beyond,
                const std::size_t range,
                const KDOPTraits&);

    template<class Traversal_traits, class Query>
    void traversal(const Query& query,
                   Traversal_traits& traits,
                   const std::size_t nb_primitives) const;

  private:
    typedef KDOPTraits KDOP_traits;
    typedef KDOP_node<KDOP_traits> Node;
    typedef typename KDOP_traits::Primitive Primitive;

  public:
    const Node& left_child() const { return *static_cast<Node*>(m_p_left_child); }
    const Node& right_child() const { return *static_cast<Node*>(m_p_right_child); }
    const Primitive& left_data() const { return *static_cast<Primitive*>(m_p_left_child); }
    const Primitive& right_data() const { return *static_cast<Primitive*>(m_p_right_child); }

  private:
    Node& left_child() { return *static_cast<Node*>(m_p_left_child); }
    Node& right_child() { return *static_cast<Node*>(m_p_right_child); }
    Primitive& left_data() { return *static_cast<Primitive*>(m_p_left_child); }
    Primitive& right_data() { return *static_cast<Primitive*>(m_p_right_child); }

  private:
    // node kdop
    Kdop m_kdop;

    // children nodes, either pointing towards children (if children are not leaves),
    // or pointing toward input primitives (if children are leaves)
    void *m_p_left_child;
    void *m_p_right_child;

  }; // end class KDOP_node

  template<typename Tr>
  template<typename ConstPrimitiveIterator>
  void
  KDOP_node<Tr>::expand(ConstPrimitiveIterator first,
                        ConstPrimitiveIterator beyond,
                        const std::size_t range,
                        const Tr& traits)
  {
    m_kdop = traits.compute_bbox_object()(first, beyond);

    //TODO sort primitives as AABB does?

    switch(range)
    {
    case 2:
      m_p_left_child = &(*first);
      m_p_right_child = &(*(++first));
      break;
    case 3:
      m_p_left_child = &(*first);
      m_p_right_child = static_cast<Node*>(this) + 1;
      right_child().expand(first + 1, beyond, 2, traits);
      break;
    default:
      const std::size_t new_range = range/2;
      m_p_left_child = static_cast<Node*>(this) + 1;
      m_p_right_child = static_cast<Node*>(this) + new_range;
      left_child().expand(first, first + new_range, new_range, traits);
      right_child().expand(first + new_range, beyond, range - new_range, traits);
    }
  }

  template<typename Tr>
  template<class Traversal_traits, class Query>
  void
  KDOP_node<Tr>::traversal(const Query& query,
                           Traversal_traits& traits,
                           const std::size_t nb_primitives) const
  {
    // recursive traversal
    switch(nb_primitives)
    {
    case 2:
      traits.intersection(query, left_data());
      if ( traits.go_further() ) {
        traits.intersection(query, right_data());
      }
      break;
    case 3:
      traits.intersection(query, left_data());
      if ( traits.go_further() && traits.do_intersect(query, right_child()) ) {
        right_child().traversal(query, traits, 2);
      }
      break;
    default:
      if ( traits.do_intersect(query, left_child()) ) {
        left_child().traversal(query, traits, nb_primitives/2);
        if ( traits.go_further() && traits.do_intersect(query, right_child()) ) {
          right_child().traversal(query, traits, nb_primitives - nb_primitives/2);
        }
      }
      else if ( traits.do_intersect(query, right_child()) ) {
        right_child().traversal(query, traits, nb_primitives - nb_primitives/2);
      }
    }
  }

} // end namespace KDOP
} // end namespace CGAL

#endif // CGAL_KDOP_TREE_INTERNAL_KDOP_NODE_H_
