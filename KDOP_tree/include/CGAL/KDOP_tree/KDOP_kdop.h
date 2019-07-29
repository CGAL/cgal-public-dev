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

#ifndef CGAL_KDOP_TREE_KDOP_KDOP_H_
#define CGAL_KDOP_TREE_KDOP_KDOP_H_

#include <vector>
#include <algorithm>
#include <array>

#include <CGAL/double.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel/Same_uncertainty.h>
#include <CGAL/Coercion_traits.h>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {
namespace KDOP_tree {

/// \addtogroup PkgKDOPTree
/// @{

  /**
   * Class KDOP_kdop is a data structure to compute and store k-dop of primitives.
   */

  template<typename GeomTraits, unsigned int N>
  class KDOP_kdop
  {
  public:
    const static unsigned int num_directions = N;

    typedef GeomTraits R;

    typedef typename R::FT FT;

    typedef typename R::Point_3 Point_3;
    typedef typename R::Segment_3 Segment_3;
    typedef typename R::Triangle_3 Triangle_3;
    typedef typename R::Ray_3 Ray_3;
    typedef typename R::Sphere_3 Sphere_3;

    /// \name Types
    /// @{

    /// Type of directions
    typedef Point_3 Direction_type;

    /// Type of vector of k directions
    typedef std::vector<Direction_type> Vec_direction;

    /// Type of support heights
    typedef std::array<FT, N> Array_height;

    typedef std::array< std::pair<FT, FT>, N > Array_height_ray;

    typedef KDOP_kdop Kdop;

    /// @}

    /// \name Constructors
    /// @{

    /// Default constructor with default directions
    KDOP_kdop() { }

    /// Constructor with support heights
    KDOP_kdop(const Array_height& support_heights)
      : array_heights_(support_heights) { }

    /// Constructor with support heights defined for rays
    KDOP_kdop(const Array_height_ray& support_heights)
      : array_heights_ray_(support_heights) { }

    /// @}

    /// Inline function to set support heights in all directions.
    void set_support_heights(const Array_height& support_heights) { array_heights_ = support_heights; }

    /// Inline function to return support heights in all directions.
    const Array_height& support_heights() const { return array_heights_; }

    const Array_height_ray& support_heights_ray() const { return array_heights_ray_; }

    /// Inline function to return the minimum support height.
    double min_height() const { return *std::min_element( array_heights_, array_heights_ + num_directions ); }

    /// Inline function to return the maximum support height.
    double max_height() const { return *std::max_element( array_heights_, array_heights_ + num_directions ); }

    //-------------------------------------------------------------------------
    class Do_overlap
    {
      const KDOP_kdop<GeomTraits, N> * m_kdop;
    public:
      Do_overlap(const KDOP_kdop<GeomTraits, N> * kdop) : m_kdop(kdop) { }

      bool operator () (const Array_height& support_heights, const Triangle_3& t) const {
        return m_kdop->do_overlap_kdop(support_heights);
      }

      bool operator () (const Array_height& support_heights, const Ray_3& r) const {
        return m_kdop->do_overlap_ray(support_heights);
      }

      bool operator () (const Array_height& support_heights, const FT& squared_radius) const {
        return m_kdop->do_overlap_sphere(support_heights, squared_radius);
      }

    };

    Do_overlap do_overlap_object() const { return Do_overlap(this); }

    //-------------------------------------------------------------------------
    // Primitive except ray and line (segment)
    bool do_overlap_kdop(const Array_height& support_heights) const;

    // Ray
    bool do_overlap_ray(const Array_height& support_heights) const;

    // Sphere
    bool do_overlap_sphere(const Array_height& support_heights, const FT& squared_radius) const;

  private:
    std::array<FT, num_directions> array_heights_;

    std::array< std::pair<FT, FT>, num_directions > array_heights_ray_; // <source, second_point> heights of a ray

  };

  template<typename GeomTraits, unsigned int N>
  bool KDOP_kdop<GeomTraits, N>::do_overlap_kdop(const Array_height& support_heights) const
  {
    bool is_overlap = true;

    const Array_height& support_heights_query = this->support_heights();

    for (int i = 0; i < num_directions/2; ++i) {
      if ( support_heights[i] + support_heights_query[i] < 0. ) return false;
      if ( support_heights[i + num_directions/2] + support_heights_query[i + num_directions/2] < 0. ) return false;
    }

    return is_overlap;
  }

  template<typename GeomTraits, unsigned int N>
  bool KDOP_kdop<GeomTraits, N>::do_overlap_ray(const Array_height& support_heights) const
  {
    bool is_overlap = true;

    const Array_height_ray& support_heights_ray = this->support_heights_ray();

    int num_inside_source = 0, num_inside_target = 0;

    typedef typename Coercion_traits<double, FT>::Type CFT;

    CFT dmin = 0., dmax = 0., rmin = 0., rmax = 0.;

    bool is_non_parallel_occur = true;
    bool is_non_parallel_first = false;
    int num_parallels = 0;

    for (int i = 0; i < num_directions/2; ++i) { // consider half the number of directions
      const FT height_source = support_heights_ray[i].first;
      const FT height_target = support_heights_ray[i].second;

      if (height_target >= height_source &&
          height_source > support_heights[i]) { // ray must outside the i-th direction
        return false;
      }

      // the opposite direction
      if (-height_target >= -height_source &&
          -height_source > support_heights[i + num_directions/2]) { // ray must outside the (i + N/2)-th direction
        return false;
      }

      CFT dmin_dir = 0., dmax_dir = 0.;
      CFT rmin_dir = 0., rmax_dir = 0.;

      if (height_target > height_source) {
        dmax_dir = support_heights[i] - height_source;
        dmin_dir = -support_heights[i + num_directions/2] - height_source;

        rmax_dir = height_target - height_source;
        rmin_dir = rmax_dir;

        if ( is_non_parallel_occur == false ) {
          is_non_parallel_first = true;
          is_non_parallel_occur = true;
        }
        else {
          is_non_parallel_first = false;
        }
      }
      else if (height_target < height_source) {
        dmin_dir = -support_heights[i] + height_source;
        dmax_dir = support_heights[i + num_directions/2] + height_source;

        rmin_dir = -height_target + height_source;
        rmax_dir = rmin_dir;

        if ( is_non_parallel_occur == false ) {
          is_non_parallel_first = true;
          is_non_parallel_occur = true;
        }
        else {
          is_non_parallel_first = false;
        }
      }
      else {
        // height_source = height_target, i.e. the ray is parallel to the i-th direction.
        // This case has already been checked previously, so just inherit tmin and tmax.
        dmin_dir = dmin;
        dmax_dir = dmax;
        rmin_dir = rmin;
        rmax_dir = rmax;

        if (num_parallels == i) {
          is_non_parallel_occur = false;
          num_parallels += 1;
        }
      }

      if (i == 0 || is_non_parallel_first == true) {
        dmin = dmin_dir;
        dmax = dmax_dir;
        rmin = rmin_dir;
        rmax = rmax_dir;
      }
      else {
        if ( is_non_parallel_occur == true && is_non_parallel_first == false ) {
          // no overlapping if tmin > tmax_dir or tmax < tmin_dir
          if ( dmin*rmax_dir > rmin*dmax_dir || dmax*rmin_dir < rmax*dmin_dir ) return false;

          // if tmin < tmin_dir, tmin = tmin_dir (narrowing the scope of t)
          if ( dmin*rmin_dir < rmin*dmin_dir ) dmin = dmin_dir, rmin = rmin_dir;

          // if tmax > tmax_dir, tmax = tmax_dir (narrowing the scope of t)
          if ( dmax*rmax_dir > rmax*dmax_dir ) dmax = dmax_dir, rmax = rmax_dir;
        }
      }

    }

    return is_overlap;
  }

  template<typename GeomTraits, unsigned int N>
  bool KDOP_kdop<GeomTraits, N>::do_overlap_sphere(const Array_height& support_heights,
                                                   const FT& squared_radius) const
  {
    bool is_overlap = true;

    const Array_height& support_heights_point = this->support_heights();

    typedef typename Coercion_traits<double, FT>::Type CFT;

    CFT distance = 0.;

    for (int i = 0; i < 3; ++i) {
      CFT d = 0.;
      if ( support_heights_point[i] > support_heights[i] ) {
        d = support_heights_point[i] - support_heights[i]; // distance between the point and the slab
        if (d * d > squared_radius) return false;
        distance += d * d; // compute Cartesian distance
        if (distance > squared_radius) return false;
      }
      else if ( -support_heights_point[i] > support_heights[i + num_directions/2] ) {
        d = -support_heights_point[i] - support_heights[i + num_directions/2];
        if (d * d > squared_radius) return false;
        distance += d * d;
        if (distance > squared_radius) return false;
      }
    }

    //if (distance > squared_radius) return false;

    for (int i = 3; i < num_directions/2; ++i) {
      CFT d = 0.;
      if ( support_heights_point[i] > support_heights[i] ) {
        d = support_heights_point[i] - support_heights[i]; // distance between the point and the slab
        if (d * d > squared_radius) return false; // compare point/slab distance and squared radius of the sphere
      }
      else if ( -support_heights_point[i] > support_heights[i + num_directions/2] ) {
        d = -support_heights_point[i] - support_heights[i + num_directions/2];
        if (d * d > squared_radius) return false; // compare point/slab distance and squared radius of the sphere
      }
    }

    return is_overlap;
  }

/// @}

  template<typename GeomTraits, unsigned int N>
  class Construct_kdop
  {
    const static unsigned int num_directions = N;

    typedef GeomTraits R;

    typedef typename R::FT FT;

    typedef typename R::Vector_3 Vector_3;

    typedef typename R::Point_3 Point_3;
    typedef typename R::Segment_3 Segment_3;
    typedef typename R::Triangle_3 Triangle_3;
    typedef typename R::Ray_3 Ray_3;
    typedef typename R::Sphere_3 Sphere_3;

  public:
    typedef std::vector< Vector_3 > Vec_direction;

    typedef std::array<FT, N> Array_height;
    typedef std::array< std::pair<FT, FT>, N > Array_height_ray;

    typedef CGAL::KDOP_tree::KDOP_kdop<GeomTraits, N> result_type;

    // a point
    result_type
    operator () (const Point_3& p) const {
      Array_height heights;
      switch(N)
      {
      case 6:
      {
        heights[0] = p.x(); heights[1] = p.y(); heights[2] = p.z();
        break;
      }
      case 14:
      {
        const double& norm = 1./std::sqrt(3.);
        heights[0] = p.x(); heights[1] = p.y(); heights[2] = p.z();
        heights[3] = (p.x() + p.y() + p.z())*norm;
        heights[4] = (-p.x() + p.y() + p.z())*norm;
        heights[5] = (-p.x() - p.y() + p.z())*norm;
        heights[6] = (p.x() - p.y() + p.z())*norm;
        break;
      }
      case 18:
      {
        const double& norm = 1./std::sqrt(2.);
        heights[0] = p.x(); heights[1] = p.y(); heights[2] = p.z();
        heights[3] = (p.x() + p.y())*norm; heights[4] = (p.x() + p.z())*norm; heights[5] = (p.y() + p.z())*norm;
        heights[6] = (p.x() - p.y())*norm; heights[7] = (p.x() - p.z())*norm; heights[8] = (p.y() - p.z())*norm;
        break;
      }
      case 26:
      {
        const double& norm1 = 1./std::sqrt(3.);
        const double& norm2 = 1./std::sqrt(2.);
        heights[0] = p.x(); heights[1] = p.y(); heights[2] = p.z();
        heights[3] = (p.x() + p.y() + p.z())*norm1;
        heights[4] = (-p.x() + p.y() + p.z())*norm1;
        heights[5] = (-p.x() - p.y() + p.z())*norm1;
        heights[6] = (p.x() - p.y() + p.z())*norm1;
        heights[7] = (p.x() + p.y())*norm2; heights[8] = (p.x() + p.z())*norm2; heights[9] = (p.y() + p.z())*norm2;
        heights[10] = (p.x() - p.y())*norm2; heights[11] = (p.x() - p.z())*norm2; heights[12] = (p.y() - p.z())*norm2;
        break;
      }
      default:
        std::cerr << "Number of directions should be 6, 14, 18 or 26!" << std::endl;
      }

      result_type kdop(heights);
      return kdop;
    }

    // a ray
    result_type
    operator () (const Ray_3& r) const {
      const Point_3& source = r.source();
      const Point_3& second_point = r.second_point();

      typename CGAL::KDOP_tree::Construct_kdop<R, N> construct_kdop;

      result_type kdop_source = construct_kdop(source);
      result_type kdop_second_point = construct_kdop(second_point);

      const Array_height& heights_source = kdop_source.support_heights();
      const Array_height& heights_second_point = kdop_second_point.support_heights();

      Array_height_ray heights_ray;
      for (int i = 0; i < num_directions/2; ++i) {
        const FT& height_source = heights_source[i];
        const FT& height_second_point = heights_second_point[i];
        heights_ray[i] = std::make_pair( height_source, height_second_point );
      }

      result_type kdop(heights_ray);
      return kdop;
    }

    // a triangle
    result_type
    operator () (const Triangle_3& t) const {
      Array_height heights;
      for (int i = 0; i < 3; ++i) {
        const Point_3& v = t.vertex(i);

        typename CGAL::KDOP_tree::Construct_kdop<R, N> construct_kdop;

        result_type kdop_vertex = construct_kdop(v);
        const Array_height& heights_vertex = kdop_vertex.support_heights();

        for (int j = 0; j < num_directions/2; ++j) {
          const FT& height = heights_vertex[j];
          if (i == 0) {
            heights[j] = height;
            heights[j + num_directions/2] = -height;
          }
          else {
            if (heights[j] < height) heights[j] = height;
            if (heights[j + num_directions/2] < -height) heights[j + num_directions/2] = -height;
          }
        }

      }

      result_type kdop(heights);
      return kdop;
    }

  };

}
}



#endif // CGAL_KDOP_TREE_KDOP_KDOP_H_
