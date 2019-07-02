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
    /// \todo Define default directions for some selected numbers N.
    KDOP_kdop() { }

    /// @}

    /// Inline function to set support heights in all directions.
    void set_support_heights(const Array_height& support_heights) { array_height_ = support_heights; }

    /// Inline function to return support heights in all directions.
    const Array_height& support_heights() const { return array_height_; }

    const Array_height_ray& support_heights_ray() const { return array_heights_ray_; }

    /// Function to compute support heights in all directions.
    void compute_support_heights(const Vec_direction& directions, const Triangle_3& t);

    /// Function to compute support heights of a ray in all directions.
    void compute_support_heights_ray(const Vec_direction& directions, const Ray_3& r);

    /// Inline function to return the minimum support height.
    double min_height() const { return *std::min_element( array_height_, array_height_ + num_directions ); }

    /// Inline function to return the maximum support height.
    double max_height() const { return *std::max_element( array_height_, array_height_ + num_directions ); }

#ifdef TEST_
    /// Inline function to set directions
    void set_directions(const Vec_direction& vec_direction) { vector_direction_.clear(); vector_direction_ = vec_direction; }

    /// Inline function to return directions.
    const Vec_direction& directions() const { return vector_direction_; }

    /*!
     * @brief Add a new direction to existing directions.
     * @param new_direction the new direction to be added
     * @return updated vector of (k + 1) directions
     */
    void add_direction(Direction_type new_direction) { vector_direction_.push_back(new_direction); }
#endif

    /*!
     * @brief Check if a line segment overlaps by comparing support heights of
     * the two k-dops.
     * @param q query
     * @param directions k-dop directions
     * @return true if the query overlaps the kdop; otherwise, false.
     * \todo support heights of the query should be located before traversal (need discussion)
     */
    template<typename Query>
    bool do_overlap_segment(const Query& q, const Vec_direction& directions) const;

    // Ray
    // \todo support heights of the query should be located before traversal (need discussion)
    bool do_overlap_ray(const Kdop& kdop_query, const Vec_direction& directions) const;


  private:
    std::array<FT, num_directions> array_height_;

    std::array< std::pair<FT, FT>, num_directions > array_heights_ray_; // store <source, second_point> heights of a ray

  };

  template<typename GeomTraits, unsigned int N>
  void KDOP_kdop<GeomTraits, N>::compute_support_heights(const Vec_direction& directions, const Triangle_3& t)
  {
    for (int i = 0; i < num_directions/2; ++i) { // consider half the number of directions
      const Point_3& direction = directions[i];

      for (int j = 0; j < 3; ++j) { // number of vertices
        //const Point_3& v = t.vertex(j);
        //double height = v.x()*direction.x() + v.y()*direction.y() + v.z()*direction.z();

        //---------------------------------------------------------------------
        FT height = 0.;
        switch(i) {
        case 0: height = v.x(); break;
        case 1: height = v.y(); break;
        case 2: height = v.z(); break;
        case 3: height = v.x() + v.y() + v.z(); break;
        case 4: height = -v.x() + v.y() + v.z(); break;
        case 5: height = -v.x() - v.y() + v.z(); break;
        case 6: height = v.x() - v.y() + v.z(); break;
        }
        //---------------------------------------------------------------------

        if ( (j == 0) || array_height_[i] < height ) {
          array_height_[i] = height;
        }
        if ( (j == 0) || array_height_[i + num_directions/2] < -height ) {
          array_height_[i + num_directions/2] = -height;
        }
      }

    }

  }

  template<typename GeomTraits, unsigned int N>
  void KDOP_kdop<GeomTraits, N>::compute_support_heights_ray(const Vec_direction& directions, const Ray_3& r)
  {
    const Point_3& source = r.source();
    const Point_3& target = r.second_point();

    for (int i = 0; i < num_directions/2; ++i) { // consider half the number of directions
      //const Point_3& direction = directions[i];

      //double height_source = source.x()*direction.x() + source.y()*direction.y() + source.z()*direction.z();
      //double height_target = target.x()*direction.x() + target.y()*direction.y() + target.z()*direction.z();

      //-----------------------------------------------------------------------
      FT height_source = 0., height_target = 0.;
      switch(i) {
      case 0: height_source = source.x(); height_target = target.x(); break;
      case 1: height_source = source.y(); height_target = target.y(); break;
      case 2: height_source = source.z(); height_target = target.z(); break;
      case 3:
      {
        height_source = source.x() + source.y() + source.z();
        height_target = target.x() + target.y() + target.z();
        break;
      }
      case 4:
      {
        height_source = -source.x() + source.y() + source.z();
        height_target = -target.x() + target.y() + target.z();
        break;
      }
      case 5:
      {
        height_source = -source.x() - source.y() + source.z();
        height_target = -target.x() - target.y() + target.z();
        break;
      }
      case 6:
      {
        height_source = source.x() - source.y() + source.z();
        height_target = target.x() - target.y() + target.z();
        break;
      }
      }
      //-----------------------------------------------------------------------

      std::pair<FT, FT> height_ray = std::make_pair(height_source, height_target);
      array_heights_ray_[i] = height_ray;

      std::pair<FT, FT> height_ray_opposite = std::make_pair(-height_source, -height_target);
      array_heights_ray_[i + num_directions/2] = height_ray_opposite;
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
  bool KDOP_kdop<GeomTraits, N>::do_overlap_ray(const Kdop& kdop_query, const Vec_direction& directions) const
  {
    bool is_overlap = true;

    const Array_height& array_heights = this->support_heights();

    const Array_height_ray& array_heights_ray = kdop_query.support_heights_ray();

    int num_inside_source = 0, num_inside_target = 0;

    typedef typename Coercion_traits<double, FT>::Type CFT;

    //---------------------------------------------------------------------------------------------
    //double dmin = 0., dmax = 0., rmin = 0., rmax = 0.;
    CFT dmin = 0., dmax = 0., rmin = 0., rmax = 0.;

    bool is_non_parallel_occur = true;
    bool is_non_parallel_first = false;
    int num_parallels = 0;
    //----------------------------------------------------------------------------------------------

    for (int i = 0; i < num_directions/2; ++i) { // consider half the number of directions
      const FT height_source = array_heights_ray[i].first;
      const FT height_target = array_heights_ray[i].second;

      if (height_target >= height_source &&
          height_source > array_heights[i]) { // ray must outside the i-th direction
        return false;
      }

      // the opposite direction
      if (-height_target >= -height_source &&
          -height_source > array_heights[i + num_directions/2]) { // ray must outside the (i + N/2)-th direction
        return false;
      }

      //-------------------------------------------------------------------------------------------
      //double dmin_dir = 0., dmax_dir = 0.;
      //double rmin_dir = 0., rmax_dir = 0.;
      CFT dmin_dir = 0., dmax_dir = 0.;
      CFT rmin_dir = 0., rmax_dir = 0.;

      if (height_target > height_source) {
        dmax_dir = array_heights[i] - height_source;
        dmin_dir = -array_heights[i + num_directions/2] - height_source;

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
        dmin_dir = -array_heights[i] + height_source;
        dmax_dir = array_heights[i + num_directions/2] + height_source;

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

      //----------------------------------------------------------------------------------------------
/*
      if (num_inside_source == i || num_inside_target == i) {
        if (height_source <= array_heights[i] &&
            -height_source <= array_heights[i + num_directions/2]) {
          num_inside_source += 1;
        }

        if (height_target <= array_heights[i] &&
            -height_target <= array_heights[i + num_directions/2]) {
          num_inside_target += 1;
        }
      }
*/
    }
/*
    // definitely inside
    if (num_inside_source == num_directions/2 ||
        num_inside_target == num_directions/2) { // ray must intersect the k-dop
      return true;
    }
*/
#ifdef DEBUG_
    // other more complex cases
    // the following code encodes: tmin = dmin/rmin; tmax = dmax/rmax;
    // where, tmin and tmax are the scope of ray parameters intersecting the slabs defined by the directions.

    //double tmin = 0., tmax = 0.; // the smallest intersection interval between the ray and the slabs
    double dmin = 0., dmax = 0., rmin = 0., rmax = 0.;

    bool is_non_parallel_occur = true;
    bool is_non_parallel_first = false;
    int num_parallels = 0;

    for (int i = 0; i < num_directions/2; ++i) { // only consider half of the directions
      const double height_source = array_heights_ray[i].first;
      const double height_target = array_heights_ray[i].second;

      double dmin_dir = 0., dmax_dir = 0.;
      double rmin_dir = 0., rmax_dir = 0.;

      if (height_target > height_source) {
        dmax_dir = array_heights[i] - height_source;
        dmin_dir = -array_heights[i + num_directions/2] - height_source;

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
        /*
         * should take the negatives of the following results, convenient for later usage.
        dmin_dir = array_heights[i] - height_source;
        dmax_dir = -array_heights[i + num_directions/2] - height_source;

        r_dir = height_target - height_source;
        */

        dmin_dir = -array_heights[i] + height_source;
        dmax_dir = array_heights[i + num_directions/2] + height_source;

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
#endif
    return is_overlap;
  }

/// @}

}
}



#endif // CGAL_KDOP_TREE_KDOP_KDOP_H_
