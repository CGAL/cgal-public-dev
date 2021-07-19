// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Camille Wormser, Pierre Alliez, Stephane Tayeb

#ifndef CGAL_AABB_NODE_H
#define CGAL_AABB_NODE_H

#include <CGAL/license/AABB_tree.h>


#include <CGAL/Profile_counter.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <vector>
#include <CGAL/Iterator_range.h>

namespace CGAL {

/**
 * @class AABB_node
 *
 *
 */
  template<typename AABBTraits>
  class AABB_node {
  private:
    typedef AABB_node<AABBTraits> Self;

  public:
    typedef typename AABBTraits::Bounding_box Bounding_box;

    static const size_t N = 2;

    /// Constructor
    AABB_node()
            : m_bbox(), m_contents(static_cast<Primitive *>(nullptr)) {};

    AABB_node(Self &&node) = default;

    // Disabled copy constructor & assignment operator
    AABB_node(const Self &src) = delete;

    Self &operator=(const Self &src) = delete;

    /// Returns the bounding box of the node
    const Bounding_box &bbox() const { return m_bbox; }

    /**
     * @brief General traversal query
     * @param query the query
     * @param traits the traversal traits that define the traversal behaviour
     * @param nb_primitives the number of primitive
     *
     * General traversal query. The traits class allows using it for the various
     * traversal methods we need: listing, counting, detecting intersections,
     * drawing the boxes.
     */
    template<class Traversal_traits, class Query>
    void traversal(const Query &query,
                   Traversal_traits &traits,
                   const std::size_t nb_primitives) const;

  private:
    typedef AABBTraits AABB_traits;
    typedef AABB_node<AABB_traits> Node;
    typedef typename AABB_traits::Primitive Primitive;


  public:
    /// Helper functions

    template<class T>
    void set_children(T *children) {
      m_contents = children;
    }

    void set_bbox(const Bounding_box &bbox) {
      m_bbox = bbox;
    }

    std::array<Node, N> &children() {
      return *static_cast<std::array<Node, N> *>(m_contents);
    }

    const std::array<Node, N> &children() const {
      return *static_cast<std::array<Node, N> *>(m_contents);
    }

    Primitive &data() { return *static_cast<Primitive *>(m_contents); }

    const Primitive &data() const { return *static_cast<Primitive *>(m_contents); }

    // TODO This is inefficient, perhaps nodes should know their own size?
    std::size_t num_primitives(const Node &child, std::size_t total_count) const {
      // Assumes that primitives are distributed as evenly as possible between children
      // When there are leftover primitives, they go to earlier nodes first
      std::size_t child_index = &child - children().begin();
      return num_primitives(child_index, total_count);
    }

    std::size_t num_primitives(std::size_t index, std::size_t total_count) const {
      // Assumes that primitives are distributed as evenly as possible between children
      // When there are leftover primitives, they go to earlier nodes first
      return total_count / N + (index < total_count % N);
    }

  private:
    /// node bounding box
    Bounding_box m_bbox;

    // Points to either children or a primitive
    void *m_contents;

  };  // end class AABB_node


  template<typename Tr>
  template<class Traversal_traits, class Query>
  void
  AABB_node<Tr>::traversal(const Query &query,
                           Traversal_traits &traits,
                           const std::size_t nb_primitives) const {

    // This is a Depth-first traversal

    if (nb_primitives == 1) {

      // If there's only one primitive, we've reached a leaf node
      traits.intersection(query, data());

    } else {

      // Otherwise, recursively search the child nodes
      for (const auto &child : children()){
        if (traits.go_further() && traits.do_intersect(query, child))
          child.traversal(query, traits, num_primitives(child, nb_primitives));
      }
    }
  }

} // end namespace CGAL

#endif // CGAL_AABB_NODE_H
