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

#ifndef KDOP_TREE_INCLUDE_CGAL_KDOP_TREE_INTERNAL_KDOP_RAY_INTERSECTION_H_
#define KDOP_TREE_INCLUDE_CGAL_KDOP_TREE_INTERNAL_KDOP_RAY_INTERSECTION_H_

#include <CGAL/license/AABB_tree.h>

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

namespace CGAL {
namespace KDOP_tree {
namespace internal {

  template<typename KDOPTree, typename SkipFunctor>
  class KDOP_ray_intersection
  {
    typedef typename KDOPTree::KDOP_traits KDOP_traits;
    typedef typename KDOP_traits::Ray_3 Ray;
    typedef typename KDOPTree::template Intersection_and_primitive_id<Ray>::Type Ray_intersection_and_primitive_id;
    typedef typename Ray_intersection_and_primitive_id::first_type Ray_intersection;

    typedef typename KDOPTree::FT FT;
    typedef typename KDOPTree::Point Point;
    typedef typename KDOPTree::Node Node;
    typedef typename KDOPTree::size_type size_type;

    typedef typename KDOP_traits::Kdop Kdop;

    typedef typename std::pair<Ray, Kdop> RayPair;

  public:
    KDOP_ray_intersection(const KDOPTree& tree) : tree_(tree) {}

    boost::optional< Ray_intersection_and_primitive_id >
    ray_intersection(const RayPair& query_pair, SkipFunctor skip) const {

      boost::optional< Ray_intersection_and_primitive_id > intersection, p;

      typename KDOP_traits::Intersection intersection_obj = tree_.traits().intersection_object();
      typename KDOP_traits::Intersection_distance intersection_distance_obj = tree_.traits().intersection_distance_object();

      typedef
#if BOOST_VERSION >= 105000
          boost::heap::priority_queue< Node_ptr_with_ft, boost::heap::compare< std::greater<Node_ptr_with_ft> > >
#else
      std::priority_queue< Node_ptr_with_ft>
#endif
      Heap_type;

      Heap_type pq;

      as_ray_param_visitor param_visitor = as_ray_param_visitor( &(query_pair.first) );

      FT t = (std::numeric_limits<double>::max)();

      // start with the root node
      pq.push(Node_ptr_with_ft(tree_.root_node(), 0, tree_.size()));

      while ( !pq.empty() && pq.top().value < t ) {
        Node_ptr_with_ft current = pq.top();
        pq.pop();

        switch(current.nb_primitives)
        {
        case 2:
        {
          // left child
          if ( !skip(current.node->left_data().id()) ) {
            intersection = intersection_obj(query_pair.first, current.node->left_data());
            if (intersection) {
              FT ray_distance = boost::apply_visitor(param_visitor, intersection->first);
              if (ray_distance < t) {
                t = ray_distance;
                p = intersection;
              }
            }
          }

          // right child
          if ( !skip(current.node->right_data().id()) ) {
            intersection = intersection_obj(query_pair.first, current.node->right_data());
            if (intersection) {
              FT ray_distance = boost::apply_visitor(param_visitor, intersection->first);
              if (ray_distance < t) {
                t = ray_distance;
                p = intersection;
              }
            }
          }
          break;
        }
        case 3:
        {
          // left child
          if ( !skip(current.node->left_data().id()) ) {
            intersection = intersection_obj(query_pair.first, current.node->left_data());
            if (intersection) {
              FT ray_distance = boost::apply_visitor(param_visitor, intersection->first);
              if (ray_distance < t) {
                t = ray_distance;
                p = intersection;
              }
            }
          }

          // right child
          const Node* child = &(current.node->right_child());
          boost::optional<FT> dist = intersection_distance_obj(query_pair.second, child->support_heights());
          if (dist) {
            pq.push(Node_ptr_with_ft(child, *dist, 2));
          }
          break;
        }
        default:
        {
          const Node* child_left = &(current.node->left_child());
          boost::optional<FT> dist = intersection_distance_obj(query_pair.second, child_left->support_heights());
          if (dist) {
            pq.push(Node_ptr_with_ft(child_left, *dist, current.nb_primitives/2));
          }

          const Node* child_right = &(current.node->right_child());
          dist = intersection_distance_obj(query_pair.second, child_right->support_heights());
          if (dist) {
            pq.push(Node_ptr_with_ft(child_right, *dist, current.nb_primitives - current.nb_primitives/2));
          }
          break;
        }
        }
      }

      return p;
    }

  private:
    const KDOPTree& tree_;

    struct Node_ptr_with_ft {
      Node_ptr_with_ft(const Node* node, const FT& value, size_type nb_primitives)
      : node(node), nb_primitives(nb_primitives), value(value) {}
      const Node* node;
      size_type nb_primitives;
      FT value;
#if BOOST_VERSION >= 105000
      bool operator<(const Node_ptr_with_ft& other) const { return value < other.value; }
      bool operator>(const Node_ptr_with_ft& other) const { return value > other.value; }
#else
      bool operator>(const Node_ptr_with_ft& other) const { return value < other.value; }
      bool operator<(const Node_ptr_with_ft& other) const { return value > other.value; }

#endif
    };

    struct as_ray_param_visitor {
      typedef FT result_type;
      as_ray_param_visitor(const Ray* ray)
      : ray_(ray), max_i_(0)
      {
        typename KDOP_traits::Geom_traits::Vector_3 v = ray->to_vector();
        for (int i = 1; i < 3; ++i)
          if( CGAL::abs(v[i]) > CGAL::abs(v[max_i_]) )
            max_i_ = i;
      }

      template<typename T>
      FT operator()(const T& s)
      {
        // intersection is a segment, returns the min relative distance
        // of its endpoints
        FT r1 = this->operator()(s[0]);
        FT r2 = this->operator()(s[1]);
        return (std::min)(r1, r2);
      }

      FT operator()(const Point& point) {
        typename KDOP_traits::Geom_traits::Vector_3 x(ray_->source(), point);
        typename KDOP_traits::Geom_traits::Vector_3 v = ray_->to_vector();

        return x[max_i_] / v[max_i_];
      }

      const Ray* ray_;
      int max_i_;
    };

  };

} // namespace internal
} // namespace KDOP_tree
} // namespace CGAL

#endif /* KDOP_TREE_INCLUDE_CGAL_KDOP_TREE_INTERNAL_KDOP_RAY_INTERSECTION_H_ */
