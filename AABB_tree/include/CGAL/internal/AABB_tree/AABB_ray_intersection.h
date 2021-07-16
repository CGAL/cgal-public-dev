// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Philipp Moeller
//

#ifndef CGAL_AABB_RAY_INTERSECTION_H
#define CGAL_AABB_RAY_INTERSECTION_H

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

  template<typename AABBTree, typename SkipFunctor>
  class AABB_ray_intersection {
    typedef typename AABBTree::AABB_traits AABB_traits;
    typedef typename AABB_traits::Ray_3 Ray;
    typedef typename AABBTree::template Intersection_and_primitive_id<Ray>::Type Ray_intersection_and_primitive_id;
    typedef typename Ray_intersection_and_primitive_id::first_type Ray_intersection;
  public:
    AABB_ray_intersection(const AABBTree &tree) : tree_(tree) {}

    boost::optional<Ray_intersection_and_primitive_id>
    ray_intersection(const Ray &query, SkipFunctor skip) const {

      // A heap is used to examine the nearest nodes first
#if BOOST_VERSION >= 105000
      typedef boost::heap::priority_queue<Node_ptr_with_ft, boost::heap::compare<std::greater<Node_ptr_with_ft> > > Heap_type;
#else
      typedef std::priority_queue< Node_ptr_with_ft> Heap_type;
#endif

      // This is our priority queue
      Heap_type pq;

      // Create the functor we'll use to check for intersections
      typename AABB_traits::Intersection intersection_obj =
              tree_.traits().intersection_object();

      // This functor determines distance, in addition to intersections
      typename AABB_traits::Intersection_distance intersection_distance_obj =
              tree_.traits().intersection_distance_object();

      // This functor estimates the distance between its ray and any line segment
      // TODO Is this a good description?
      as_ray_param_visitor param_visitor = as_ray_param_visitor(&query);

      // Our result, and an intermediate used during calculation
      // TODO Perhaps the separate intermediate value could be eliminated?
      boost::optional<Ray_intersection_and_primitive_id> intersection, p;

      // "t" begins with an effectively infinite value
      // this is not the right way to do it, but using numeric_limits<FT>::{max,infinity} will not work with Epeck.
      // TODO std::numeric_limits can have custom specializations, maybe we should create one for Epeck!
      // https://stackoverflow.com/questions/16122912/is-it-ok-to-specialize-stdnumeric-limitst-for-user-defined-number-like-class
      FT t = (std::numeric_limits<double>::max)();

      // The root node should be evaluated first
      pq.push(Node_ptr_with_ft(tree_.root_node(), 0, tree_.size()));

      // Continue to evaluate nodes until the heap is empty, or the closest node can't be closer than t
      while (!pq.empty() && pq.top().value < t) {

        // Retrieve the next node to be evaluated
        Node_ptr_with_ft current = pq.top();
        pq.pop();

        if (current.nb_primitives == 1) {
          // If this node is a leaf, directly check for intersection

          // Only check this node if we're not told to skip it
          if (!skip(current.node->data().id())) {

            // Check if it intersects with the query ray, using the functor we created before
            intersection = intersection_obj(query, current.node->data());
            if (intersection) {

              // If it did intersect, use our distance heuristic to estimate how close the first intersection is
              FT ray_distance = boost::apply_visitor(param_visitor, intersection->first);

              // If this intersection is closer than the previous best, update the closest intersection
              if (ray_distance < t) {
                t = ray_distance;
                p = intersection;
              }
            }
          }

        } else {
          // If the node has children, add each of them to the queue to be checked later
          for (const auto &child : current.node->children()) {
            boost::optional<FT> dist = intersection_distance_obj(query, child.bbox());
            if (dist) pq.push(Node_ptr_with_ft(&child, *dist, current.node->num_primitives(child, current.nb_primitives)));
          }
        }
      }

      // p should now hold the closest intersection, and its associated primitive ID
      return p;
    }

  private:
    const AABBTree &tree_;
    typedef typename AABBTree::Point Point;
    typedef typename AABBTree::FT FT;
    typedef typename AABBTree::Node Node;
    typedef typename AABBTree::size_type size_type;

    struct Node_ptr_with_ft {
      Node_ptr_with_ft(const Node *node, const FT &value, size_type nb_primitives)
              : node(node), nb_primitives(nb_primitives), value(value) {}

      const Node *node;
      size_type nb_primitives;
      FT value;
#if BOOST_VERSION >= 105000

      bool operator<(const Node_ptr_with_ft &other) const { return value < other.value; }

      bool operator>(const Node_ptr_with_ft &other) const { return value > other.value; }

#else
      bool operator>(const Node_ptr_with_ft& other) const { return value < other.value; }
      bool operator<(const Node_ptr_with_ft& other) const { return value > other.value; }

#endif
    };

    struct as_ray_param_visitor {
      typedef FT result_type;

      as_ray_param_visitor(const Ray *ray)
              : ray(ray), max_i(0) {
        typename AABB_traits::Geom_traits::Vector_3 v = ray->to_vector();
        for (int i = 1; i < 3; ++i)
          if (CGAL::abs(v[i]) > CGAL::abs(v[max_i]))
            max_i = i;
      }

      template<typename T>
      FT operator()(const T &s) {
        // intersection is a segment, returns the min relative distance
        // of its endpoints
        FT r1 = this->operator()(s[0]);
        FT r2 = this->operator()(s[1]);
        return (std::min)(r1, r2);
      }

      FT operator()(const Point &point) {
        typename AABB_traits::Geom_traits::Vector_3 x(ray->source(), point);
        typename AABB_traits::Geom_traits::Vector_3 v = ray->to_vector();

        return x[max_i] / v[max_i];
      }

      const Ray *ray;
      int max_i;
    };
  };

  template<typename AABBTraits>
  template<typename Ray, typename SkipFunctor>
  boost::optional<typename AABB_tree<AABBTraits>::template Intersection_and_primitive_id<Ray>::Type>
  AABB_tree<AABBTraits>::first_intersection(const Ray &query,
                                            const SkipFunctor &skip) const {
    CGAL_static_assertion_msg((boost::is_same<Ray, typename AABBTraits::Ray_3>::value),
                              "Ray and Ray_3 must be the same type");

    switch (size()) // copy-paste from AABB_tree::traversal
    {
      case 0: // Tree empty, nothing to intersect
        break;
      case 1: // Tree has 1 node, intersect directly
        return traits().intersection_object()(query, singleton_data());
      default: // Tree has >= 2 nodes
        if (traits().do_intersect_object()(query, root_node()->bbox())) {
          AABB_ray_intersection<AABB_tree<AABBTraits>, SkipFunctor> ri(*this);
          return ri.ray_intersection(query, skip);
        } else {
          // but we don't hit the root
          break;
        }
    }
    return boost::none;
  }

  template<typename AABBTraits>
  template<typename Ray, typename SkipFunctor>
  boost::optional<typename AABB_tree<AABBTraits>::Primitive_id>
  AABB_tree<AABBTraits>::first_intersected_primitive(const Ray &query,
                                                     const SkipFunctor &skip) const {
    boost::optional<
            typename AABB_tree<AABBTraits>::
            template Intersection_and_primitive_id<Ray>::Type> res =
            first_intersection(query, skip);
    if ((bool) res)
      return boost::make_optional(res->second);
    return boost::none;
  }

}

#endif /* CGAL_AABB_RAY_INTERSECTION_H */
