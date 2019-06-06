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

/// \file KDOP_node.h

namespace CGAL {
namespace KDOP_tree {

/// \addtogroup PkgKDOPTree
/// @{

/*! K-dop tree node class
 * \tparam KDOPTraits is a model of the concept \ref KDOPTraits.
 */

  template<typename KDOPTraits>
  class KDOP_node
  {
  public:

    /// \name Types
    /// @{

    /// Type of k-dop
    typedef typename KDOPTraits::Kdop Kdop;

    /// @}

    /// \name Constructor
    /// @{

    /// Null constructor
    KDOP_node()
      : m_kdop()
      , m_p_left_child(NULL)
      , m_p_right_child(NULL) { };

    /// Non-virtual destructor
    ~KDOP_node() { };

    /// @}

    /// \name Functions
    /// @{

    /// return the kdop of the node
    const Kdop& kdop() const { return m_kdop; }

    /*!
     * @brief Build the tree by recursive expansion.
     * @param first the first primitive to insert
     * @param last the last primitive to insert
     * @param range the number of primitive of the range
     *
     * [first,last[ is the range of primitives to be added to the tree.
     *
     * \todo Add the recursive code without computing k-dops; consider to create
     * an octree or a binary tree.
     */
    template<typename ConstPrimitiveIterator>
    void expand(ConstPrimitiveIterator first,
                ConstPrimitiveIterator beyond,
                const std::size_t range,
                const KDOPTraits&);

    /*!
     * @brief General traversal query
     * @param query the query
     * @param traits the traversal traits that define the traversal behaviour
     * @param nb_primitives the number of primitive
     *
     * General traversal query. The traits class allows using it for the various
     * traversal methods we need: listing, counting, computing k-dops, detecting
     * intersections.
     *
     * \todo Add recursive code to traverse the tree.
     */
    template<typename Traversal_traits, typename Query>
    void traversal(const Query& query,
                   Traversal_traits& traits,
                   const std::size_t nb_primitives) const;

    /// @}

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
    //todo may not need to compute k-dop in this process
    m_kdop = traits.compute_kdop_object()(first, beyond);

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
  template<typename Traversal_traits, typename Query>
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

  /// @}

} // end namespace KDOP
} // end namespace CGAL

#endif // CGAL_KDOP_TREE_INTERNAL_KDOP_NODE_H_
