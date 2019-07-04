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

    /// @}

    /// Inline function to set support heights in all directions.
    void set_support_heights(const Array_height& support_heights) { array_heights_ = support_heights; }

    /// Inline function to return support heights in all directions.
    const Array_height& support_heights() const { return array_heights_; }

    const Array_height_ray& support_heights_ray() const { return array_heights_ray_; }

    /// Function to compute support heights in all directions.
    void compute_support_heights(const Vec_direction& directions, const Triangle_3& t);

    void compute_support_heights_vertex(const Point_3& vertex, Array_height& heights);

    /// Function to compute support heights of a ray in all directions.
    void compute_support_heights_ray(const Vec_direction& directions, const Ray_3& r);

    /// Inline function to return the minimum support height.
    double min_height() const { return *std::min_element( array_heights_, array_heights_ + num_directions ); }

    /// Inline function to return the maximum support height.
    double max_height() const { return *std::max_element( array_heights_, array_heights_ + num_directions ); }

    /*!
     * @brief Check if a line segment overlaps by comparing support heights of
     * the two k-dops.
     * @param q query
     * @param directions k-dop directions
     * @return true if the query overlaps the kdop; otherwise, false.
     */
    template<typename Query>
    bool do_overlap_segment(const Query& q, const Vec_direction& directions) const;

    // Ray
    bool do_overlap_ray(const Array_height& support_heights) const;

  private:
    std::array<FT, num_directions> array_heights_;

    std::array< std::pair<FT, FT>, num_directions > array_heights_ray_; // store <source, second_point> heights of a ray

  };

  template<typename GeomTraits, unsigned int N>
  void KDOP_kdop<GeomTraits, N>::compute_support_heights(const Vec_direction& directions, const Triangle_3& t)
  {
/*
    for (int i = 0; i < num_directions/2; ++i) { // consider half the number of directions
      const Point_3& direction = directions[i];

      for (int j = 0; j < 3; ++j) { // number of vertices
        const Point_3& v = t.vertex(j);
        FT height = v.x()*direction.x() + v.y()*direction.y() + v.z()*direction.z();

        if ( (j == 0) || array_heights_[i] < height ) {
          array_heights_[i] = height;
        }
        if ( (j == 0) || array_heights_[i + num_directions/2] < -height ) {
          array_heights_[i + num_directions/2] = -height;
        }
      }
    }
*/

    for (int i = 0; i < 3; ++i) {
      const Point_3& v = t.vertex(i);

      Array_height vertex_heights;
      this->compute_support_heights_vertex(v, vertex_heights);

      for (int j = 0; j < num_directions/2; ++j) {
        const double& height = vertex_heights[j];
        if (i == 0) {
          array_heights_[j] = height;
          array_heights_[j + num_directions/2] = -height;
        }
        else {
          if (array_heights_[j] < height) array_heights_[j] = height;
          if (array_heights_[j + num_directions/2] < -height) array_heights_[j + num_directions/2] = -height;
        }
      }
    }

  }

  template<typename GeomTraits, unsigned int N>
  void KDOP_kdop<GeomTraits, N>::compute_support_heights_ray(const Vec_direction& directions, const Ray_3& r)
  {
    const Point_3& source = r.source();
    const Point_3& target = r.second_point();
/*
    for (int i = 0; i < num_directions/2; ++i) { // consider half the number of directions
      const Point_3& direction = directions[i];

      FT height_source = source.x()*direction.x() + source.y()*direction.y() + source.z()*direction.z();
      FT height_target = target.x()*direction.x() + target.y()*direction.y() + target.z()*direction.z();

      std::pair<FT, FT> height_ray = std::make_pair(height_source, height_target);
      array_heights_ray_[i] = height_ray;

      std::pair<FT, FT> height_ray_opposite = std::make_pair(-height_source, -height_target);
      array_heights_ray_[i + num_directions/2] = height_ray_opposite;
    }
*/

    Array_height heights_source, heights_target;

    this->compute_support_heights_vertex(source, heights_source);
    this->compute_support_heights_vertex(target, heights_target);

    for (int i = 0; i < num_directions/2; ++i) {
      const double& height_source = heights_source[i];
      const double& height_target = heights_target[i];
      array_heights_ray_[i] = std::make_pair( height_source, height_target );
      //array_heights_ray_[i + num_directions/2] = std::make_pair( -height_source, -height_target );
    }

  }

  template<typename GeomTraits, unsigned int N>
  template<typename Query>
  bool KDOP_kdop<GeomTraits, N>::do_overlap_segment(const Query& q, const Vec_direction& directions) const
  {
    bool is_overlap = true;

    const Point_3& source = q.source();
    const Point_3& target = q.target();

    const Array_height& array_heights = this->support_heights();

    Array_height array_heights_query;
    for (int i = 0; i < num_directions; ++i) {
      const Point_3& direction = directions[i];

      double height1 = -source.x()*direction.x() - source.y()*direction.y() - source.z()*direction.z();
      double height2 = -target.x()*direction.x() - target.y()*direction.y() - target.z()*direction.z();

      if (height1 <= height2) array_heights_query[i] = height2;
      else array_heights_query[i] = height1;

      if (array_heights[i] + array_heights_query[i] < 0) return false; // line segment outside the i-th direction

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

  //---------------------------------------------------------------------------
  // compute support heights with fixed directions
  template<typename GeomTraits, unsigned N>
  void KDOP_kdop<GeomTraits, N>::compute_support_heights_vertex(const Point_3& v, Array_height& heights)
  {
    switch(N)
    {
    case 6:
    {
      heights[0] = v.x(); heights[1] = v.y(); heights[2] = v.z();
      break;
    }
    case 14:
    {
      heights[0] = v.x(); heights[1] = v.y(); heights[2] = v.z();
      heights[3] = v.x() + v.y() + v.z();
      heights[4] = -v.x() + v.y() + v.z();
      heights[5] = -v.x() - v.y() + v.z();
      heights[6] = v.x() - v.y() + v.z();
      break;
    }
    case 18:
    {
      heights[0] = v.x(); heights[1] = v.y(); heights[2] = v.z();
      heights[3] = v.x() + v.y(); heights[4] = v.x() + v.z(); heights[5] = v.y() + v.z();
      heights[6] = v.x() - v.y(); heights[7] = v.x() - v.z(); heights[8] = v.y() - v.z();
      break;
    }
    case 26:
    {
      heights[0] = v.x(); heights[1] = v.y(); heights[2] = v.z();
      heights[3] = v.x() + v.y() + v.z();
      heights[4] = -v.x() + v.y() + v.z();
      heights[5] = -v.x() - v.y() + v.z();
      heights[6] = v.x() - v.y() + v.z();
      heights[7] = v.x() + v.y(); heights[8] = v.x() + v.z(); heights[9] = v.y() + v.z();
      heights[10] = v.x() - v.y(); heights[11] = v.x() - v.z(); heights[12] = v.y() - v.z();
      break;
    }
    default:
      std::cerr << "Number of directions shoule be 6, 14, 18 or 26!" << std::endl;
    }
  };

/// @}

}
}



#endif // CGAL_KDOP_TREE_KDOP_KDOP_H_
