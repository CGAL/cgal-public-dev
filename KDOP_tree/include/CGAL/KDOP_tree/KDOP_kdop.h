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
    const static unsigned int direction_number = N;

    typedef GeomTraits R;

    typedef typename R::Point_3 Point_3;
    typedef typename R::Segment_3 Segment_3;
    typedef typename R::Triangle_3 Triangle_3;

    /// \name Types
    /// @{

    /// Type of directions
    typedef Point_3 Direction_type;

    /// Type of vector of k directions
    typedef std::vector<Direction_type> Vec_direction;

    /// Type of support heights
    typedef std::array<double, N> Array_height;

    typedef KDOP_kdop Kdop;

    /// @}

    /// \name Constructors
    /// @{

    /// Default constructor with default directions
    /// \todo Define default directions for some selected numbers N.
    KDOP_kdop() { }

    /// Constructor with directions given.
    KDOP_kdop(Vec_direction vector_direction)
      : vector_direction_(vector_direction)
    {}

    /// @}

    /// Inline function to set support heights in all directions.
    void set_support_heights(const Array_height& support_heights) { array_height_ = support_heights; }

    /// Inline function to return support heights in all directions.
    Array_height support_heights() const& { return array_height_; }

    /// Function to compute support heights in all directions.
    void compute_support_heights(const Vec_direction& directions, const Triangle_3& t);

    /// Inline function to return the minimum support height.
    double min_height() const { return *std::min_element( array_height_, array_height_ + direction_number ); }

    /// Inline function to return the maximum support height.
    double max_height() const { return *std::max_element( array_height_, array_height_ + direction_number ); }

#ifdef DEBUG_
    /// Inline function to set directions
    void set_directions(const Vec_direction& vec_direction) { vector_direction_.clear(); vector_direction_ = vec_direction; }

    /// Inline function to return directions.
    Vec_direction directions() const& { return vector_direction_; }

    /*!
     * @brief Add a new direction to existing directions.
     * @param new_direction the new direction to be added
     * @return updated vector of (k + 1) directions
     */
    void add_direction(Direction_type new_direction) { vector_direction_.push_back(new_direction); }
#endif

    /*!
     * @brief Check if two k-dops overlap by comparing support heights of
     * the two k-dops.
     * @param kdop1 the first k-dop
     * @param kdop2 the second k-dop
     * @return true if the two k-dops overlap; otherwise, false.
     * \todo Add the checking function.
     */
    bool do_overlap(const Kdop& kdop1, const Kdop& kdop2);

  private:
    Vec_direction vector_direction_; // Redundant at the moment

    std::array<double, direction_number> array_height_;

  };

  template<typename GeomTraits, unsigned N>
  void KDOP_kdop<GeomTraits, N>::compute_support_heights(const Vec_direction& directions, const Triangle_3& t)
  {
    int num_directions = directions.size();

    for (int i = 0; i < num_directions; ++i) { // number of directions
      Point_3 direction = directions[i];

      std::cout << "direction: " << direction.x() << ", " << direction.y() << ", " << direction.z() << std::endl;

      std::vector<double> heights;
      for (int j = 0; j < 3; ++j) { // number of vertices
        Point_3 v = t.vertex(j);
        double height = v.x()*direction.x() + v.y()*direction.y() + v.z()*direction.z();
        heights.push_back(height);

        std::cout << "vertex: " << v.x() << ", " << v.y() << ", " << v.z() << ": height = " << height << std::endl;
      }

      double height_max = *std::max_element(heights.begin(), heights.end());

      array_height_[i] = height_max;
    }
  }

  template<typename GeomTraits, unsigned int N>
  bool KDOP_kdop<GeomTraits, N>::do_overlap(const Kdop& kdop1, const Kdop& kdop2)
  {
    bool is_overlap = true;

    std::array<double, direction_number> support_heights1 = kdop1.support_heights();
    std::array<double, direction_number> support_heights2 = kdop2.support_heights();

    int num_support_heights = direction_number;

    int num_non_overlap = 0;
    for (int i = 0; i < num_support_heights; ++i) {
      double height1 = support_heights1[i];
      double height2 = support_heights2[i];

      if (height1 + height2 < 0) num_non_overlap += 1; // can "break" after this to reduce checks!
    }

    if (num_non_overlap != 0) is_overlap = false;

    return is_overlap;
  }

/// @}

}
}



#endif // CGAL_KDOP_TREE_KDOP_KDOP_H_
