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
    typedef std::vector< std::vector<double> > Vec_direction;

    /// Type of support heights
    typedef std::vector<double> Vec_height;

    /// @}

    /// \name Constructors
    /// @{

    /// Default constructor (directions inferred from template parameter N)
    KDOP_kdop() { }

    /// Constructor with directions given
    KDOP_kdop(Vec_direction vector_direction)
      : vector_direction_(vector_direction)
    {}

    /// @}

    //TODO some flexibility can be provided for users to choose whether to use user-defined directions or pre-computed directions.

    /// inline function to compute the minimum support height
    inline double min_height() const { return *std::min_element( vector_height_.begin(), vector_height_.end() ); }

    /// inline function to compute the maximum support height
    inline double max_height() const { return *std::max_element( vector_height_.begin(), vector_height_.end() ); }

  private:
    Vec_direction vector_direction_;

    std::vector<double> vector_height_;

  };

/// @}

}
}



#endif // CGAL_KDOP_TREE_KDOP_KDOP_H_
