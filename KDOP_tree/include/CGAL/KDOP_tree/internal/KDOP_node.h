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
namespace internal {

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

    /// Type of bounding box
    typedef typename KDOPTraits::Bounding_box Bounding_box;

    /// Type of k-dop
    typedef typename KDOPTraits::Kdop Kdop;
    typedef typename Kdop::Vec_direction Vec_direction;
    typedef typename Kdop::Array_height Array_height;

    const static unsigned int num_directions = Kdop::num_directions;

    /// @}

    /// \name Constructor
    /// @{

    /// Null constructor
    KDOP_node()
      : m_support_heights()
      , m_p_left_child(NULL)
      , m_p_right_child(NULL) { };

    /// Non-virtual destructor
    ~KDOP_node() { };

    /// @}

    /// \name Functions
    /// @{

    /// return support heights of the node
    const Array_height& support_heights() const { return m_support_heights; }

    /// return k-DOP of the node
    Kdop kdop() const {
      Kdop kdop(m_support_heights);
      return kdop;
    }

    /*!
     * @brief Build the tree by recursive expansion.
     * @param first the first primitive to insert
     * @param last the last primitive to insert
     * @param range the number of primitive of the range
     *
     * [first,last[ is the range of primitives to be added to the tree.
     *
     */
    template<typename ConstPrimitiveIterator>
    void expand(ConstPrimitiveIterator first,
                ConstPrimitiveIterator beyond,
                const std::size_t range,
                const KDOPTraits&);

    template<typename Traversal_traits>
    void kdop_traversal(Traversal_traits& traits,
                        const std::size_t nb_primitives,
                        const Vec_direction& directions);

    void union_support_heights(const Array_height& left_height,
                               const Array_height& right_height,
                               Array_height& height_union);

    /*!
     * @brief General traversal query
     * @param query_pair the query and its k-dop
     * @param traits the traversal traits that define the traversal behaviour
     * @param nb_primitives the number of primitive
     *
     * General traversal query. The traits class allows using it for the various
     * traversal methods we need: listing, counting, computing k-dops, detecting
     * intersections.
     *
     */
    template<typename Traversal_traits, typename QueryPair>
    void traversal(const QueryPair& query_pair,
                   Traversal_traits& traits,
                   const std::size_t nb_primitives) const;

    template<typename Traversal_traits>
    void kdop_heights(const Traversal_traits& traits,
                      const std::size_t nb_primitives,
                      const Vec_direction& directions,
                      std::vector< Array_height >& heights);

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
    // node bounding box
    Bounding_box m_bbox;

    // node support heights
    Array_height m_support_heights;

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
    // binary splitting as AABB
    m_bbox = traits.compute_bbox_object()(first, beyond);

    traits.split_primitives_object()(first, beyond, m_bbox);

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
  void
  KDOP_node<Tr>::union_support_heights(const Array_height& left_heights,
                                       const Array_height& right_heights,
                                       Array_height& heights_union)
  {
    for (int i = 0; i < num_directions/2; ++i) { // consider half the number of directions
      if (left_heights[i] >= right_heights[i]) heights_union[i] = left_heights[i];
      else heights_union[i] = right_heights[i];

      // opposite direction
      if (left_heights[i + num_directions/2] >= right_heights[i + num_directions/2]) {
        heights_union[i + num_directions/2] = left_heights[i + num_directions/2];
      }
      else {
        heights_union[i + num_directions/2] = right_heights[i + num_directions/2];
      }
    }
  }

  template<typename Tr>
  template<typename Traversal_traits>
  void
  KDOP_node<Tr>::kdop_traversal(Traversal_traits& traits,
                                const std::size_t nb_primitives,
                                const Vec_direction& directions)
  {
    // recursive traversal
    switch(nb_primitives)
    {
    case 2:
    {
      Kdop left_leaf_kdop = traits.compute_kdop(left_data(), directions);
      Kdop right_leaf_kdop = traits.compute_kdop(right_data(), directions);

      const Array_height& left_support_heights = left_leaf_kdop.support_heights();
      const Array_height& right_support_heights = right_leaf_kdop.support_heights();

      // union of support heights of two children
      this->union_support_heights(left_support_heights, right_support_heights,
                                  m_support_heights);
    }
    break;
    case 3:
    {
      Kdop left_leaf_kdop = traits.compute_kdop(left_data(), directions);

      const Array_height& left_support_heights = left_leaf_kdop.support_heights();

      right_child().kdop_traversal(traits, 2, directions);

      const Array_height& right_support_heights = right_child().support_heights();

      // union of support heights of two children
      this->union_support_heights(left_support_heights, right_support_heights,
                                  m_support_heights);
    }
    break;
    default:
      left_child().kdop_traversal(traits, nb_primitives/2, directions);
      right_child().kdop_traversal(traits, nb_primitives - nb_primitives/2, directions);

      const Array_height& left_support_heights = left_child().support_heights();
      const Array_height& right_support_heights = right_child().support_heights();

      this->union_support_heights(left_support_heights, right_support_heights,
                                  m_support_heights);
    }
  }

  template<typename Tr>
  template<typename Traversal_traits, typename QueryPair>
  void
  KDOP_node<Tr>::traversal(const QueryPair& query_pair,
                           Traversal_traits& traits,
                           const std::size_t nb_primitives) const
  {
    // recursive traversal
    switch(nb_primitives)
    {
    case 2:
      traits.intersection(query_pair.first, left_data());
      if ( traits.go_further() ) {
        traits.intersection(query_pair.first, right_data());
      }
      break;
    case 3:
      traits.intersection(query_pair.first, left_data());
      if ( traits.go_further() && traits.do_intersect(query_pair.first, query_pair.second, right_child()) ) {
        right_child().traversal(query_pair, traits, 2);
      }
      break;
    default:
      if ( traits.do_intersect(query_pair.first, query_pair.second, left_child()) ) {
        left_child().traversal(query_pair, traits, nb_primitives/2);
        if ( traits.go_further() && traits.do_intersect(query_pair.first, query_pair.second, right_child()) ) {
          right_child().traversal(query_pair, traits, nb_primitives - nb_primitives/2);
        }
      }
      else if ( traits.do_intersect(query_pair.first, query_pair.second, right_child()) ) {
        right_child().traversal(query_pair, traits, nb_primitives - nb_primitives/2);
      }
    }
  }

  template<typename Tr>
  template<typename Traversal_traits>
  void
  KDOP_node<Tr>::kdop_heights(const Traversal_traits& traits,
                              const std::size_t nb_primitives,
                              const Vec_direction& directions,
                              std::vector< Array_height >& heights)
  {
    switch(nb_primitives)
    {
    case 2:
    {
      // do nothing;
    }
      break;
    case 3:
    {
      const Array_height& kdop_heights = this->support_heights();
      heights.push_back(kdop_heights);
    }
      break;
    default:
      const Array_height& kdop_heights = this->support_heights();
      heights.push_back(kdop_heights);

      left_child().kdop_heights(traits, nb_primitives/2, directions, heights);
      right_child().kdop_heights(traits, nb_primitives - nb_primitives/2, directions, heights);
    }
  }

  /// @}

} // end namespace internal
} // end namespace KDOP
} // end namespace CGAL

#endif // CGAL_KDOP_TREE_INTERNAL_KDOP_NODE_H_
