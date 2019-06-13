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

    /// @}

    /// \name Constructor
    /// @{

    /// Null constructor
    KDOP_node()
      : m_kdop()
      , m_p_left_child(NULL)
      , m_p_right_child(NULL)
      , m_directions()
      , m_direction_number() { };

    /// Non-virtual destructor
    ~KDOP_node() { };

    /// @}

    /// \name Functions
    /// @{

    /// store the kdop of the node
    void set_kdop(const Kdop& kdop) { m_kdop = kdop; }

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

    template<typename Traversal_traits>
    void kdop_traversal(Traversal_traits& traits,
                        const std::size_t nb_primitives,
                        const Vec_direction& directions,
                        const int direction_number);

    void union_support_heights(const std::vector<double>& left_height,
                               const std::vector<double>& right_height,
                               std::vector<double>& height_union,
                               const int direction_number);

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
    // node bounding box
    Bounding_box m_bbox;

    // node kdop
    Kdop m_kdop;
    Vec_direction m_directions;

    int m_direction_number;

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

#ifdef TEST_
  template<typename Tr>
  template<typename ConstPrimitiveIterator>
  void
  KDOP_node<Tr>::expand(ConstPrimitiveIterator first,
                        ConstPrimitiveIterator beyond,
                        const std::size_t range,
                        const Tr& traits)
  {
    // binary splitting as AABB
    traits.split_primitives_object()(first, beyond);

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
#endif

  template<typename Tr>
  void
  KDOP_node<Tr>::union_support_heights(const std::vector<double>& left_heights,
                                       const std::vector<double>& right_heights,
                                       std::vector<double>& heights_union,
                                       const int direction_number)
  {
    for (int i = 0; i < direction_number; ++i) {
      double left_height = left_heights[i];
      double right_height = right_heights[i];

      if (left_height >= right_height) heights_union.push_back(left_height);
      else heights_union.push_back(right_height);
    }
  }

  template<typename Tr>
  template<typename Traversal_traits>
  void
  KDOP_node<Tr>::kdop_traversal(Traversal_traits& traits,
                                const std::size_t nb_primitives,
                                const Vec_direction& directions,
                                const int direction_number
                                )
  {
    // recursive traversal
    switch(nb_primitives)
    {
    case 2:
    {
      Kdop left_leaf_kdop = traits.compute_kdop(left_data(), directions, direction_number);

      Kdop right_leaf_kdop = traits.compute_kdop(right_data(), directions, direction_number);

      std::vector<double> left_support_heights = left_leaf_kdop.give_support_heights();
      std::vector<double> right_support_heights = right_leaf_kdop.give_support_heights();

      // union of support heights of two children
      std::vector<double> union_support_heights; // union of support heights in all directions

      this->union_support_heights(left_support_heights, right_support_heights,
                                  union_support_heights, direction_number);

      Kdop kdop(directions);

      kdop.set_support_heights(union_support_heights);

      this->set_kdop(kdop);

    }
      break;
    default:
      left_child().kdop_traversal(traits, nb_primitives/2, directions, direction_number);
      right_child().kdop_traversal(traits, nb_primitives - nb_primitives/2, directions, direction_number);

      Kdop left_kdop = left_child().kdop();
      Kdop right_kdop = right_child().kdop();

      std::vector<double> left_support_heights = left_kdop.give_support_heights();
      std::vector<double> right_support_heights = right_kdop.give_support_heights();

      // union of support heights of two children
      std::vector<double> union_support_heights; // union of support heights in all directions

      this->union_support_heights(left_support_heights, right_support_heights,
                                  union_support_heights, direction_number);

      typename Kdop::Vec_direction vec_direction = left_kdop.give_directions();

      Kdop kdop(vec_direction);

      kdop.set_support_heights(union_support_heights);

      this->set_kdop(kdop);

      std::cout << "union support heights: " << std::endl;
      for (int i = 0; i < direction_number; ++i) {
        std::cout << union_support_heights[i] << std::endl;
      }
      std::cout << std::endl;

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

} // end namespace internal
} // end namespace KDOP
} // end namespace CGAL

#endif // CGAL_KDOP_TREE_INTERNAL_KDOP_NODE_H_
