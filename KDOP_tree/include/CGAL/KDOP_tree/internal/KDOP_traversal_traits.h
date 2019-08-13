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
   * \ingroup PkgKDOPTreeTraitsClasses
   * @class Any_intersection_traits
   * Traits used in the k-dop tree traversal to get the first intersection obtained in the tree traversal,
   * which is different from `First_intersection_traits` which is to get the intersection closest to the
   * source of a ray.
   */
  template<typename KDOPTraits, typename Query>
  class Any_intersection_traits
  {
    typedef typename KDOPTraits::Primitive Primitive;
    typedef typename KDOPTraits::Kdop Kdop;

    typedef CGAL::KDOP_tree::internal::KDOP_node<KDOPTraits> Node;

  public:
    typedef boost::optional< typename KDOPTraits::template Intersection_and_primitive_id<Query>::Type > Result_type;

    Any_intersection_traits(const KDOPTraits& traits)
      : m_result(), m_traits(traits) {}

    bool go_further() const { return !m_result; }

    void intersection(const Query& query, const Primitive& primitive)
    {
      m_result = m_traits.intersection_object()(query, primitive);
    }

    bool do_intersect(const Query& query, const Kdop& kdop_query, const Node& node) const
    {
      return m_traits.do_intersect_object()(query, kdop_query, node.support_heights());
    }

    Result_type result() const { return m_result; }

    bool is_intersection_found() const { return m_result; }

  private:
    Result_type m_result;
    const KDOPTraits& m_traits;
  };

  /*!
   * \ingroup PkgKDOPTreeTraitsClasses
   * @class Any_primitive_traits
   * Traits used in the k-dop tree traversal to get first intersected primitive.
   */
  template<typename KDOPTraits, typename Query>
  class Any_primitive_traits
  {
    typedef typename KDOPTraits::Primitive Primitive;
    typedef typename Primitive::Id Primitive_id;

    typedef typename KDOPTraits::Kdop Kdop;

    typedef CGAL::KDOP_tree::internal::KDOP_node<KDOPTraits> Node;

  public:
    Any_primitive_traits(const KDOPTraits& traits)
  : m_is_found(false)
  , m_result()
  , m_traits(traits) {}

    bool go_further() const { return !m_is_found; }

    void intersection(const Query& query, const Primitive& primitive)
    {
      if ( m_traits.do_intersect_object()(query, primitive) ) {
        m_result = boost::optional<Primitive_id>(primitive.id());
        m_is_found = true;
      }
    }

    bool do_intersect(const Query& query, const Kdop& kdop_query, const Node& node) const
    {
      return m_traits.do_intersect_object()(query, kdop_query, node.support_heights());
    }

    boost::optional<Primitive_id> result() const { return m_result; }

    bool is_intersection_found() const { return m_is_found; }

  private:
    bool m_is_found;
    boost::optional<Primitive_id> m_result;
    const KDOPTraits& m_traits;
  };

  /*!
   * \ingroup PkgKDOPTreeTraitsClasses
   * @class Listing_intersection_traits
   * Traits used in the k-dop tree traversal to get intersections.
   */
  template<typename KDOPTraits, typename Query, typename OutputIterator>
  class Listing_intersection_traits
  {
    typedef typename KDOPTraits::Primitive Primitive;
    typedef typename KDOPTraits::Kdop Kdop;

    typedef CGAL::KDOP_tree::internal::KDOP_node<KDOPTraits> Node;

  public:
    Listing_intersection_traits(OutputIterator out_it, const KDOPTraits& traits)
      : m_out_it(out_it), m_traits(traits) {}

    bool go_further() const { return true; }

    void intersection(const Query& query, const Primitive& primitive)
    {
      boost::optional< typename KDOPTraits::template Intersection_and_primitive_id<Query>::Type >
      intersection = m_traits.intersection_object()(query, primitive);

      if (intersection) {
        *m_out_it++ = *intersection;
      }
    }

    bool do_intersect(const Query& query, const Kdop& kdop_query, const Node& node) const
    {
      return m_traits.do_intersect_object()(query, kdop_query, node.support_heights());
    }

  private:
    OutputIterator m_out_it;
    const KDOPTraits& m_traits;
  };

  /*!
   * \ingroup PkgKDOPTreeTraitsClasses
   * @class Listing_primitive_traits
   * Traits used in the k-dop tree traversal to get intersected primitives.
   */
  template<typename KDOPTraits, typename Query, typename OutputIterator>
  class Listing_primitive_traits
  {
    typedef typename KDOPTraits::Primitive Primitive;
    typedef typename KDOPTraits::Kdop Kdop;

    typedef CGAL::KDOP_tree::internal::KDOP_node<KDOPTraits> Node;

  public:
    Listing_primitive_traits(OutputIterator out_it, const KDOPTraits& traits)
      : m_out_it(out_it), m_traits(traits) {}

    bool go_further() const { return true; }

    void intersection(const Query& query, const Primitive& primitive)
    {
      if ( m_traits.do_intersect_object()(query, primitive) ) {
        *m_out_it++ = primitive.id();
      }
    }

    bool do_intersect(const Query& query, const Kdop& kdop_query, const Node& node) const
    {
      return m_traits.do_intersect_object()(query, kdop_query, node.support_heights());
    }

  private:
    OutputIterator m_out_it;
    const KDOPTraits& m_traits;
  };

  /*!
   * \ingroup PkgKDOPTreeTraitsClasses
   * @class Do_intersect_traits
   * Traits used in the k-dop tree traversal to check intersection.
   */
  template<typename KDOPTraits, typename Query>
  class Do_intersect_traits
  {
    typedef typename KDOPTraits::FT FT;
    typedef typename KDOPTraits::Point_3 Point;
    typedef typename KDOPTraits::Primitive Primitive;
    typedef typename KDOPTraits::Kdop Kdop;
    typedef typename Kdop::Vec_direction Vec_direction;

    typedef CGAL::KDOP_tree::internal::KDOP_node<KDOPTraits> Node;

  public:
    Do_intersect_traits(const KDOPTraits& traits)
      : m_is_found(false), m_traits(traits)
    {}

    bool go_further() const { return !m_is_found; }

    void intersection(const Query& query, const Primitive& primitive)
    {
      if ( m_traits.do_intersect_object()(query, primitive) ) {
        m_is_found = true;
      }
    }

    bool do_intersect(const Query& query, const Kdop& kdop_query, const Node& node) const
    {
      return m_traits.do_intersect_object()(query, kdop_query, node.support_heights());
    }

    bool is_intersection_found() const { return m_is_found; }

  private:
    bool m_is_found;
    const KDOPTraits& m_traits;
  };

  /*!
   * \ingroup PkgKDOPTreeTraitsClasses
   * @class Compute_kdop_traits
   * Traits used in the k-dop tree traversal to compute k-dops.
   */
  template<typename KDOPTraits>
  class Compute_kdop_traits
  {
    typedef typename KDOPTraits::FT FT;
    typedef typename KDOPTraits::Point_3 Point;
    typedef typename KDOPTraits::Primitive Primitive;
    typedef typename KDOPTraits::Kdop Kdop;

    typedef CGAL::KDOP_tree::internal::KDOP_node<KDOPTraits> Node;

    typedef typename Kdop::Vec_direction Vec_direction;

  public:
    Compute_kdop_traits(const KDOPTraits& traits)
      : m_kdop()
      , m_traits(traits) {}

    Kdop compute_kdop(const Primitive& primitive,
                      const Vec_direction& directions)
    {
      m_kdop = m_traits.compute_kdop_object()(primitive, directions);

      return m_kdop;
    }

  private:
    Kdop m_kdop;
    const KDOPTraits& m_traits;
  };

  /*!
   * \ingroup PkgKDOPTreeTraitsClasses
   * @class Projection_traits
   * Traits used for distance queries.
   */
  template<typename KDOPTraits>
  class Projection_traits
  {
    typedef typename KDOPTraits::Geom_traits Geom_traits;
    typedef typename KDOPTraits::FT FT;
    typedef typename KDOPTraits::Point_3 Point;
    typedef typename KDOPTraits::Primitive Primitive;
    typedef typename KDOPTraits::Kdop Kdop;
    typedef typename KDOPTraits::Point_and_primitive_id Point_and_primitive_id;

    typedef CGAL::KDOP_tree::internal::KDOP_node<KDOPTraits> Node;

  public:
    Projection_traits(const Point& query,
                      const Point& hint,
                      const typename Primitive::Id& hint_primitive,
                      const KDOPTraits& traits)
      : m_query(query),
        m_closest_point(hint),
        m_closest_primitive(hint_primitive),
        m_traits(traits)
    {
      m_sq_distance = Geom_traits().compute_squared_distance_3_object()(query, hint);
    }

    bool go_further() const { return true; }

    void intersection(const Point& query, const Primitive& primitive)
    {
      Point new_closest_point = m_traits.closest_point_object()(query, primitive, m_closest_point);
      if ( !m_traits.equal_3_object()(new_closest_point, m_closest_point) ) {
        m_closest_primitive = primitive.id();
        m_closest_point = new_closest_point;
        m_sq_distance = Geom_traits().compute_squared_distance_3_object()(m_query, m_closest_point);
      }
    }

    bool do_intersect(const Point& query, const Kdop& kdop_query, const Node& node) const
    {
      return m_traits.compare_distance_object()(query, kdop_query, node.support_heights(), m_sq_distance);
    }

    Point closest_point() const { return m_closest_point; }

    Point_and_primitive_id closeset_point_and_primitive() const
    {
      return Point_and_primitive_id(m_closest_point, m_closest_primitive);
    }

  private:
    Point m_query, m_closest_point;
    typename Primitive::Id m_closest_primitive;
    const KDOPTraits& m_traits;
    FT m_sq_distance;
  };

  /// @}

} // namespace KDOP_tree
} // namespace CGAL


#endif // CGAL_KDOP_TREE_INTERNAL_KDOP_TRAVERSAL_TRAITS_H_
