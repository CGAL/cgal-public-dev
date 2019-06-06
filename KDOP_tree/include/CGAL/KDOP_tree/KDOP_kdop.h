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

namespace CGAL {
namespace KDOP_tree {

  template<typename T>
  struct Simple_cartesian;

/// \addtogroup PkgKDOPTree
/// @{

  /**
   * Class KDOP_kdop is a data structure to compute and store k-dop of primitives.
   */

  template<unsigned int N>
  class KDOP_kdop
  {
  public:
    typedef Simple_cartesian<double> R;

    /// \name Types
    /// @{

    /// Type of directions
    typedef std::vector<double> Direction_type;

    /// Type of vector of k directions
    typedef std::vector<Direction_type> Vec_direction;

    /// Type of support heights
    typedef std::vector<double> Vec_height;

    /// @}

    /// \name Constructors
    /// @{

    /// Default constructor with directions inferred from template parameter N.
    /// \todo Define default directions for some selected numbers N.
    KDOP_kdop() { }

    /// Constructor with directions given.
    KDOP_kdop(Vec_direction vector_direction)
      : vector_direction_(vector_direction)
    {}

    /// @}

    /// Inline function to return support heights in all directions.
    Vec_height give_support_heights() const { return vector_height_; }

    /// Inline function to return the minimum support height.
    double min_height() const { return *std::min_element( vector_height_.begin(), vector_height_.end() ); }

    /// Inline function to return the maximum support height.
    double max_height() const { return *std::max_element( vector_height_.begin(), vector_height_.end() ); }

    /*!
     * @brief Add a new direction to existing directions.
     * @param new_direction the new direction to be added
     * @return updated vector of (k + 1) directions
     */
    void add_direction(Direction_type new_direction) { vector_direction_.push_back(new_direction); }

    /*!
     * @brief Check if two k-dops overlap by comparing support heights of
     * the two k-dops.
     * @param kdop1 the first k-dop
     * @param kdop2 the second k-dop
     * @return true if the two k-dops overlap; otherwise, false.
     * \todo Add the checking function.
     */
    bool do_overlap(const KDOP_kdop& kdop1, const KDOP_kdop& kdop2);

  private:
    Vec_direction vector_direction_;

    std::vector<double> vector_height_;

  };

/// @}

}
}



#endif // CGAL_KDOP_TREE_KDOP_KDOP_H_
