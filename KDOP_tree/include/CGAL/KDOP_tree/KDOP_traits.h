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

#include <CGAL/disable_warnings.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Default.h>
#include <CGAL/KDOP_tree/KDOP_kdop.h>

#include <CGAL/internal/AABB_tree/Has_nested_type_Shared_data.h>
#include <CGAL/internal/AABB_tree/Is_ray_intersection_geomtraits.h>
#include <CGAL/internal/AABB_tree/Primitive_helper.h>

#include <boost/optional.hpp>
#include <boost/bind.hpp>

namespace CGAL {
namespace KDOP_tree {

template<class T>
struct Remove_optional { typedef T type; };

template<class T>
struct Remove_optional< ::boost::optional<T> > { typedef T type; };

// KDOP_traits_base brings in the Intersection_distance predicate for ray queries
template <unsigned int N, typename GeomTraits, bool ray_intersection_geom_traits = CGAL::internal::AABB_tree::Is_ray_intersection_geomtraits<GeomTraits>::value>
struct KDOP_traits_base;

template <unsigned int N, typename GeomTraits>
struct KDOP_traits_base<N, GeomTraits, false> {};

template <unsigned int N, typename GeomTraits>
struct KDOP_traits_base<N, GeomTraits, true>
{
  typedef typename GeomTraits::FT    FT;
  typedef typename GeomTraits::Ray_3 Ray_3;

  typedef typename CGAL::KDOP_tree::KDOP_kdop<GeomTraits, N> Kdop;

  typedef typename Kdop::Array_height Array_height;
  typedef typename Kdop::Array_height_ray Array_height_ray;

  // similar to ray/kdop overlap check
  struct Intersection_distance {
    boost::optional<FT> operator () (const Kdop& kdop_ray, const Array_height& support_heights) const {
      FT t_min = -DBL_MAX, t_max = DBL_MAX;

      const Array_height_ray& support_heights_ray = kdop_ray.support_heights_ray();

      for (int i = 0; i < N/2; ++i) { // consider half the number of directions
        const FT height_source = support_heights_ray[i].first;
        const FT height_second_point = support_heights_ray[i].second;

        if (height_second_point >= height_source &&
            height_source > support_heights[i]) {
          return boost::none;
        }
        if (-height_second_point >= -height_source &&
            -height_source > support_heights[i + N/2]) {
          return boost::none;
        }

        // compute intersection parameter t
        if (height_second_point != height_source) {
          FT t1 = (support_heights[i] - height_source)/(height_second_point - height_source);
          FT t2 = (-support_heights[i + N/2] - height_source)/(height_second_point - height_source);

          FT t_min_dir = std::min(t1, t2);
          FT t_max_dir = std::max(t1, t2);

          if (t_min > t_max_dir || t_max < t_min_dir) return boost::none;

          t_min = std::max( t_min, t_min_dir );
          t_max = std::min( t_max, t_max_dir );

          if (t_min > t_max || t_max < FT(0.)) return boost::none;
        }
        else { // ray parallel to the direction
          // do nothing
        }
      }

      if (t_min < FT(0.)) {
        return FT(0.);
      }
      else {
        return t_min;
      }

      return t_min;
    }
  };

  Intersection_distance intersection_distance_object() const { return Intersection_distance(); }

};

/// \addtogroup PkgKDOPTree
/// @{

/// \tparam GeomTraits must  be a model of the concept \ref KDOPGeomTraits,
/// and provide the geometric types as well as the intersection tests and computations.
/// \tparam KdopPrimitive provide the type of primitives stored in the KDOP_tree.
///   It is a model of the concept `KDOPPrimitive`.
///
/// \tparam KdopMap must be a model of `ReadablePropertyMap` that has as key type a primitive id,
///                 and as value type a `Kdop`.
///                 If the type is `Default` the `Datum` must have the
///                 member function `kdop()` that returns the k-dop of the primitive.
///
/// If the argument `GeomTraits` is a model of the concept \ref
/// KdopRayIntersectionGeomTraits, this class is also a model of \ref
/// KdopRayIntersectionTraits.
///
/// \sa `KDOPTraits`
/// \sa `KDOP_tree`
/// \sa `KDOPPrimitive`

  template<unsigned int N, typename GeomTraits, typename KDOPPrimitive, typename BboxMap = Default, typename KDOPMap = Default>
  class KDOP_traits:
      public CGAL::internal::AABB_tree::AABB_traits_base<KDOPPrimitive>,
      public KDOP_traits_base<N, GeomTraits>
  {
  public:

    /// \name Types
    /// @{

    /// Type of geometry traits (kernel)
    typedef GeomTraits Geom_traits;

    /// Type of fild number of the kernel
    typedef typename GeomTraits::FT FT;

    /// Type of k-dop traits
    typedef KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap> KT;

    /// Type of primitives
    typedef KDOPPrimitive Primitive;

    /// 3D point and Primitive Id type
    typedef typename std::pair<typename GeomTraits::Point_3, typename Primitive::Id> Point_and_primitive_id;

    /// Type of intersection result
    template<typename Query>
    struct Intersection_and_primitive_id {
      typedef typename cpp11::result_of<typename GeomTraits::Intersect_3(Query, typename Primitive::Datum)>::type Intersection_type;

      typedef std::pair< typename Remove_optional<Intersection_type>::type, typename Primitive::Id > Type;
    };

    /// Type of 3D point
    typedef typename GeomTraits::Point_3 Point_3;

    /// Type of bounding box
    typedef typename CGAL::Bbox_3 Bounding_box;

    /// Type of k-dop
    typedef typename CGAL::KDOP_tree::KDOP_kdop<GeomTraits, N> Kdop;

    typedef typename Kdop::Vec_direction Vec_direction;

    typedef typename Kdop::Array_height Array_height;

    typedef typename CGAL::KDOP_tree::Construct_kdop<GeomTraits, N> Construct_kdop;

    /// Type of a sphere
    typedef typename GeomTraits::Sphere_3 Sphere_3;
    typedef typename GeomTraits::Compute_squared_radius_3 Compute_squared_radius_3;

    KDOPMap kdm;

    BboxMap bbm;

    /// @}

    /// \name Constructor
    /// @{

    /// Default constructor
    KDOP_traits() { }

    /// Constructor with given k-dop map
    KDOP_traits(BboxMap bbm, KDOPMap kdm)
      : bbm(bbm), kdm(kdm)
    {}

    /// @}

    typedef typename GeomTraits::Compute_squared_distance_3 Squared_distance;
    Squared_distance squared_distance_object() const { return GeomTraits().compute_squared_distance_3_object(); }

    typedef typename GeomTraits::Equal_3 Equal_3;
    Equal_3 equal_3_object() const { return GeomTraits().equal_3_object(); }

    /**
     * Split a range of primitives defined by [first, beyond).
     *
     * @param first iterator on the first element
     * @param beyond iterator on the past-the-end element
     *
     * \todo Split the primitives with an octree or a binary tree.
     *
     */
    class Split_primitives
    {
      typedef KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap> Traits;

      const Traits& m_traits;

    public:
      Split_primitives(const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& traits)
        : m_traits(traits) {}

      typedef void result_type;

      template<typename PrimitiveIterator>
      void operator()(PrimitiveIterator first,
                      PrimitiveIterator beyond,
                      const typename KT::Bounding_box& bbox) const
      {
#ifdef DEBUG_
        std::cout << "split primitives:" << std::endl;
        for (PrimitiveIterator pIter = first; pIter != beyond; ++pIter) {
          std::cout << (*pIter).id() << ", ";
        }
        std::cout << std::endl;
#endif

        PrimitiveIterator middle = first + (beyond - first)/2;
        switch(Traits::longest_axis(bbox))
        {
        case KT::CGAL_AXIS_X: // sort along x
          std::nth_element(first, middle, beyond, boost::bind(Traits::less_x,_1,_2,m_traits));
          break;
        case KT::CGAL_AXIS_Y: // sort along y
          std::nth_element(first, middle, beyond, boost::bind(Traits::less_y,_1,_2,m_traits));
          break;
        case KT::CGAL_AXIS_Z: // sort along z
          std::nth_element(first, middle, beyond, boost::bind(Traits::less_z,_1,_2,m_traits));
          break;
        default:
          CGAL_error();
        }
      }

    };

    Split_primitives split_primitives_object() const {return Split_primitives(*this);}

    /*
     * Compute the bounding box of a set of primitives, for the sake of
     * splitting primitives with bounding boxes.
     * @param first an iterator on the first primitive
     * @param beyond an iterator on the past-the-end primitive
     * @return the bounding box of the primitives of the iterator range
     */
    class Compute_bbox {
      const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& m_traits;
    public:
      Compute_bbox(const KDOP_traits<N, GeomTraits,KDOPPrimitive, BboxMap, KDOPMap>& traits)
        :m_traits (traits) {}

      template<typename ConstPrimitiveIterator>
      typename KT::Bounding_box operator()(ConstPrimitiveIterator first,
          ConstPrimitiveIterator beyond) const
      {
        typename KT::Bounding_box bbox = m_traits.compute_bbox(*first,m_traits.bbm);
        for(++first; first != beyond; ++first)
        {
          bbox = bbox + m_traits.compute_bbox(*first,m_traits.bbm);
        }
        return bbox;
      }

    };

    Compute_bbox compute_bbox_object() const {return Compute_bbox(*this);}


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
      const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& m_traits;
    public:
      Compute_kdop(const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& traits)
        : m_traits (traits) {}

      typename KT::Kdop operator () (const Primitive& pr,
                                     const Vec_direction& directions) const
      {
        typename KT::Kdop kdop = m_traits.compute_kdop(pr, directions);
        return kdop;
      }
    };

    Compute_kdop compute_kdop_object() const {return Compute_kdop(*this);}

    /**
     * Check if the query intersects the primitive
     *
     * @param q query object
     * @param pr primitive
     *
     * @return a bool result
     */
    class Do_intersect
    {
      const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& m_traits;
    public:
      Do_intersect(const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& traits)
    : m_traits (traits) {}

      template<typename Query>
      bool operator () (const Query& query, const Kdop& kdop_query, const Array_height& support_heights) const
      {
        bool is_intersect = m_traits.do_intersect(query, kdop_query, support_heights);
        return is_intersect;
      }

      template<typename Query>
      bool operator () (const Query& q, const Primitive& pr) const
      {
        return GeomTraits().do_intersect_3_object()(q, CGAL::internal::Primitive_helper<KT>::get_datum(pr, m_traits));
      }
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
      const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& m_traits;
    public:
      Intersection(const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& traits)
        : m_traits(traits) {}

      template<typename Query>
      boost::optional< typename Intersection_and_primitive_id<Query>::Type >
      operator () (const Query& q, const Primitive& pr) const {
        typename cpp11::result_of<typename GeomTraits::Intersect_3(Query, typename Primitive::Datum) >::type
          inter_res = GeomTraits().intersect_3_object()(CGAL::internal::Primitive_helper<KT>::get_datum(pr, m_traits), q);

        if (!inter_res) return boost::none;

        return boost::make_optional( std::make_pair(*inter_res, pr.id()) );
      }
    };

    Intersection intersection_object() const {return Intersection(*this);}

    /**
     * Compute the closest point in the primitives to a point query
     */
    class Closest_point {
      typedef typename KT::Point_3 Point;
      typedef typename KT::Primitive Primitive;

      const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& m_traits;

    public:
      Closest_point(const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& traits)
        : m_traits(traits) {}

      Point operator () (const Point& p, const Primitive& pr, const Point& bound) const
      {
        GeomTraits geom_traits;
        Point closest_point = geom_traits.construct_projected_point_3_object()(CGAL::internal::Primitive_helper<KT>::get_datum(pr, m_traits), p);

        return geom_traits.compare_distance_3_object()(p, closest_point, bound) == CGAL::LARGER ? bound : closest_point;
      }

    };

    Closest_point closest_point_object() const { return Closest_point(*this); }

    /**
     * Compare distance of two points to primitives
     */
    class Compare_distance {
      typedef typename KT::FT FT;
      typedef typename KT::Point_3 Point;
      typedef typename KT::Primitive Primitive;
      typedef typename KT::Kdop Kdop;

      typedef typename Kdop::Array_height Array_height;

    public:
      template<typename Height>
      bool operator () (const Point& p,
                        const Kdop& kdop_p,
                        const Height& support_heights,
                        const Point& bound) const
      {
        FT sq_distance = GeomTraits().compute_squared_distance_3_object()(p, bound);
        return (&kdop_p)->do_overlap_object()(support_heights, sq_distance);
      }

      template<typename Height>
      bool operator () (const Point& p,
                        const Kdop& kdop_p,
                        const Height& support_heights,
                        const FT& sq_distance) const
      {
        return (&kdop_p)->do_overlap_object()(support_heights, sq_distance);
      }

    };

    Compare_distance compare_distance_object() const { return Compare_distance(); }

  private:
    /**
     * @brief Computes bounding box of one primitive
     * @param pr the primitive
     * @return the bounding box of the primitive \c pr
     */
    template <typename PM>
    Bounding_box compute_bbox(const Primitive& pr, const PM&)const
    {
      return get(bbm, pr.id());
    }

    Bounding_box compute_bbox(const Primitive& pr, const Default&)const
    {
      return CGAL::internal::Primitive_helper<KT>::get_datum(pr,*this).bbox();
    }

    /**
     * Compute the k-dop of a primitive
     *
     * @param pr primitive
     * @param directions directions of k-dop
     * @param direction_number number of directions
     *
     * @return the k-dop of the primitive \c pr
     *
     */
    Kdop compute_kdop(const Primitive& pr,
                      const Vec_direction& directions) const
    {
      Construct_kdop construct_kdop;

      return construct_kdop( CGAL::internal::Primitive_helper<KT>::get_datum(pr, *this) );
    }

    /**
     * Check the intersection with k-dop
     *
     * @param q query
     * @param kdop k-dop of the node
     * @param directions k-dop directions
     *
     * @return true if the query intersects the node; otherwise, false. \c pr
     *
     */
    template<typename Query>
    bool do_intersect(const Query& query,
                      const Kdop& kdop_query,
                      const Array_height& support_heights) const
    {
      bool is_intersect = (&kdop_query)->do_overlap_object()(support_heights, query);

      return is_intersect;
    }

    typedef enum { CGAL_AXIS_X = 0,
                   CGAL_AXIS_Y = 1,
                   CGAL_AXIS_Z = 2 } Axis;

    static Axis longest_axis(const Bounding_box& bbox);

#ifdef TEST_
    void compute_min_max(const Primitive& pr,
                         std::vector<double>& minCoord,
                         std::vector<double>& maxCoord) const;

    static Axis longest_axis(const std::vector<double>& minCoord,
                             const std::vector<double>& maxCoord);
#endif

    /// Comparison functions
    static bool less_x(const Primitive& pr1, const Primitive& pr2, const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& traits)
    {
      return GeomTraits().less_x_3_object()( CGAL::internal::Primitive_helper<KT>::get_reference_point(pr1, traits),
                                             CGAL::internal::Primitive_helper<KT>::get_reference_point(pr2, traits) );
    }

    static bool less_y(const Primitive& pr1, const Primitive& pr2, const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& traits)
    {
      return GeomTraits().less_y_3_object()( CGAL::internal::Primitive_helper<KT>::get_reference_point(pr1, traits),
                                             CGAL::internal::Primitive_helper<KT>::get_reference_point(pr2, traits) );
    }

    static bool less_z(const Primitive& pr1, const Primitive& pr2, const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& traits)
    {
      return GeomTraits().less_z_3_object()( CGAL::internal::Primitive_helper<KT>::get_reference_point(pr1, traits),
                                             CGAL::internal::Primitive_helper<KT>::get_reference_point(pr2, traits) );
    }

  }; // end class KDOP_traits

  //---------------------------------------------------------------------------
  // Private methods
  //---------------------------------------------------------------------------
  template<unsigned int N, typename GeomTraits, typename KDOPPrimitive, typename BboxMap, typename KDOPMap>
  typename KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>::Axis
  KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>::longest_axis (const Bounding_box& bbox)
  {
    const double dx = bbox.xmax() - bbox.xmin();
    const double dy = bbox.ymax() - bbox.ymin();
    const double dz = bbox.zmax() - bbox.zmin();

    if (dx >= dy) {
      if (dx >= dz) {
        return CGAL_AXIS_X;
      }
      else {
        return CGAL_AXIS_Z;
      }
    }
    else {
      if (dy >= dz) {
        return CGAL_AXIS_Y;
      }
      else {
        return CGAL_AXIS_Z;
      }
    }

  }

/// @}

} // end namespace KDOP_tree
} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_KDOP_TREE_KDOP_TRAITS_H_
