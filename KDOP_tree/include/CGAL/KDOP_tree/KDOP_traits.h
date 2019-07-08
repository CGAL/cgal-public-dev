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
#include <CGAL/KDOP_tree/internal/Primitive_helper.h>

#include <boost/optional.hpp>
#include <boost/bind.hpp>

namespace CGAL {
namespace KDOP_tree {

template<class T>
struct Remove_optional { typedef T type; };

template<class T>
struct Remove_optional< ::boost::optional<T> > { typedef T type; };

// helper controlling whether extra data should be stored in the KDOP_tree traits class
template <class Primitive, bool has_shared_data = CGAL::internal::Has_nested_type_Shared_data<Primitive>::value>
struct KDOP_traits_base;

template <class Primitive>
struct KDOP_traits_base<Primitive,false>{};

template <class Primitive>
struct KDOP_traits_base<Primitive, true> {
  typename Primitive::Shared_data m_primitive_data;

  template <typename ... T>
  void set_shared_data(T&& ... t){
    m_primitive_data=Primitive::construct_shared_data(std::forward<T>(t)...);
  }
  const typename Primitive::Shared_data& shared_data() const {return m_primitive_data;}
};

/// \addtogroup PkgKDOPTree
/// @{

  /**
   * Traits class
   * \todo Add KDOP_traits_base class.
   */

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
      public KDOP_traits_base<KDOPPrimitive>
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

#ifdef TEST_
      template<typename PrimitiveIterator>
      void operator () (PrimitiveIterator first,
                        PrimitiveIterator beyond) const
      {
        std::vector<double> minCoord, maxCoord;

        for (PrimitiveIterator pIter = first; pIter != beyond; ++pIter) {
          m_traits.compute_min_max(*pIter, minCoord, maxCoord);

          std::cout << (*pIter).id() << ", ";
        }
        std::cout << std::endl;

        PrimitiveIterator middle = first + (beyond - first)/2;
        switch(Traits::longest_axis(minCoord, maxCoord))
        {
        case KT::CGAL_AXIS_X: // split along x
          std::nth_element(first, middle, beyond, boost::bind(Traits::less_x, _1, _2, m_traits));
          break;
        case KT::CGAL_AXIS_Y: // split along y
          std::nth_element(first, middle, beyond, boost::bind(Traits::less_y, _1, _2, m_traits));
          break;
        case KT::CGAL_AXIS_Z: // split along z
          std::nth_element(first, middle, beyond, boost::bind(Traits::less_z, _1, _2, m_traits));
          break;
        default:
          CGAL_error();
        }
      }
#endif

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
     *
     * \todo Define operators to check intersection with k-dops.
     *
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
        return GeomTraits().do_intersect_3_object()(q, internal::Primitive_helper<KT>::get_datum(pr, m_traits));
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
     * \todo Define operators to compute intersection.
     *
     */
    class Intersection
    {
      //TODO define operator to compute intersection
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
        Point closest_point = geom_traits.construct_projected_point_3_object()(internal::Primitive_helper<KT>::get_datum(pr, m_traits), p);

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

    public:
      template<typename Solid>
      CGAL::Comparison_result operator () (const Point& p, const Kdop& p_kdop, const Solid& pr, const Point& bound) const
      {
        //\todo point-polytope distance (make full use of first three axis directions)
        //return GeomTraits().do_intersect_3_object()
        //                    (GeomTraits().construct_sphere_3_object()(p, GeomTraits().compute_squared_distance_3_object()(p, bound)), pr) ?
        //                    CGAL::SMALLER : CGAL::LARGER;
      }

      template<typename Solid>
      CGAL::Comparison_result operator () (const Point& p, const Kdop& p_kdop, const Solid& pr, const FT& sq_distance) const
      {
        //\todo point-polytope distance (make full use of first three axis directions)
        //return GeomTraits().do_intersect_3_object()
        //                    (GeomTraits().construct_sphere_3_object()(p, sq_distance), pr) ?
        //                    CGAL::SMALLER : CGAL::LARGER;
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
      return internal::Primitive_helper<KT>::get_datum(pr,*this).bbox();
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
#ifdef DEBUG_
      std::cout << "primitive: " << pr.id() << std::endl;
#endif

      Kdop kdop;

      (&kdop)->compute_support_heights_object()( directions, internal::Primitive_helper<KT>::get_datum(pr, *this) );

#ifdef DEBUG_
      std::cout << std::endl;
#endif

      return kdop;
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
      return GeomTraits().less_x_3_object()( internal::Primitive_helper<KT>::get_reference_point(pr1, traits),
                                             internal::Primitive_helper<KT>::get_reference_point(pr2, traits) );
    }

    static bool less_y(const Primitive& pr1, const Primitive& pr2, const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& traits)
    {
      return GeomTraits().less_y_3_object()( internal::Primitive_helper<KT>::get_reference_point(pr1, traits),
                                             internal::Primitive_helper<KT>::get_reference_point(pr2, traits) );
    }

    static bool less_z(const Primitive& pr1, const Primitive& pr2, const KDOP_traits<N, GeomTraits, KDOPPrimitive, BboxMap, KDOPMap>& traits)
    {
      return GeomTraits().less_z_3_object()( internal::Primitive_helper<KT>::get_reference_point(pr1, traits),
                                             internal::Primitive_helper<KT>::get_reference_point(pr2, traits) );
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

#ifdef TEST_
  //todo need to generalise it to other primitives, currently only triangles.
  template<typename GT, typename KP, typename KM>
  void KDOP_traits<GT, KP, KM>::compute_min_max(const Primitive& pr,
                                                std::vector<double>& minCoord,
                                                std::vector<double>& maxCoord) const
  {
    Point_3 p1 = internal::Primitive_helper<KT>::get_datum(pr, *this).vertex(0);
    Point_3 p2 = internal::Primitive_helper<KT>::get_datum(pr, *this).vertex(1);
    Point_3 p3 = internal::Primitive_helper<KT>::get_datum(pr, *this).vertex(2);

    std::vector<double> xCoord, yCoord, zCoord;
    xCoord.push_back(p1.x()), xCoord.push_back(p2.x()), xCoord.push_back(p3.x());
    yCoord.push_back(p1.y()), yCoord.push_back(p2.y()), yCoord.push_back(p3.y());
    zCoord.push_back(p1.z()), zCoord.push_back(p2.z()), zCoord.push_back(p3.z());

    double xmin = *std::min_element(xCoord.begin(), xCoord.end());
    double xmax = *std::max_element(xCoord.begin(), xCoord.end());

    double ymin = *std::min_element(yCoord.begin(), yCoord.end());
    double ymax = *std::max_element(yCoord.begin(), yCoord.end());

    double zmin = *std::min_element(zCoord.begin(), zCoord.end());
    double zmax = *std::max_element(zCoord.begin(), zCoord.end());

    if (minCoord.empty()) {
      minCoord.push_back(xmin), minCoord.push_back(ymin), minCoord.push_back(zmin);
    }
    else {
      if (xmin < minCoord[0]) minCoord[0] = xmin;
      if (ymin < minCoord[1]) minCoord[1] = ymin;
      if (zmin < minCoord[2]) minCoord[2] = zmin;
    }

    if (maxCoord.empty()) {
      maxCoord.push_back(xmax), maxCoord.push_back(ymax), maxCoord.push_back(zmax);
    }
    else {
      if (xmax > maxCoord[0]) maxCoord[0] = xmax;
      if (ymax > maxCoord[1]) maxCoord[1] = ymax;
      if (zmax > maxCoord[2]) maxCoord[2] = zmax;
    }

  }

  template<typename GeomTraits, typename KDOPPrimitive, typename KDOPMap>
  typename KDOP_traits<GeomTraits, KDOPPrimitive, KDOPMap>::Axis
  KDOP_traits<GeomTraits, KDOPPrimitive, KDOPMap>::longest_axis
  (const std::vector<double>& minCoord, const std::vector<double>& maxCoord)
  {
    const double dx = maxCoord[0] - minCoord[0];
    const double dy = maxCoord[1] - minCoord[1];
    const double dz = maxCoord[2] - minCoord[2];

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
#endif

/// @}

} // end namespace KDOP_tree
} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_KDOP_TREE_KDOP_TRAITS_H_
