// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_NUMBER_UTILS_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_NUMBER_UTILS_H

#include <CGAL/license/Levels_of_detail.h>

// Boost headers.
#include <boost/mpl/has_xxx.hpp>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel/global_functions.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits> 
  class Default_sqrt {
    
  private:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

  public:
    FT operator()(const FT value) const { 
      
      CGAL_precondition(value >= FT(0));
      return static_cast<FT>(CGAL::sqrt(CGAL::to_double(value)));
    }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<typename GeomTraits, 
  bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
  class Get_sqrt {
        
  public:
    using Traits = GeomTraits;
    using Sqrt = Default_sqrt<Traits>;

    static Sqrt sqrt_object(const Traits& ) { 
      return Sqrt();
    }
  };

  // Case: do_not_use_default = true.
  template<typename GeomTraits>
  class Get_sqrt<GeomTraits, true> {
        
  public:
    using Traits = GeomTraits;
    using Sqrt = typename Traits::Sqrt;

    static Sqrt sqrt_object(const Traits& traits) { 
      return traits.sqrt_object();
    }
  };

  template<typename FT>
  static FT max_value() {
    return FT(1000000000000);
  }

  template<typename Vector_3>
  typename Kernel_traits<Vector_3>::Kernel::FT
  vector_length(const Vector_3& v) {

    using Traits = typename Kernel_traits<Vector_3>::Kernel;
    using FT = typename Traits::FT;
    using Get_sqrt = Get_sqrt<Traits>;
    using Sqrt = typename Get_sqrt::Sqrt;
    const Sqrt sqrt;
    return static_cast<FT>(sqrt(v * v));
  }


  template<typename Point>
  typename Kernel_traits<Point>::Kernel::FT
  distance(
    const Point& p, 
    const Point& q) {
      
    using Traits = typename Kernel_traits<Point>::Kernel;
    using FT = typename Traits::FT;
    return static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(p, q))));
  }

  template<
  typename Item_range,
  typename Point_map_2,
  typename Line_2,
  typename Point_2>
  void boundary_points_on_line_2(
    const Item_range& item_range,
    const Point_map_2& point_map_2,
    const std::vector<std::size_t>& indices,
    const Line_2& line,
    Point_2 &p,
    Point_2 &q) {

    using Traits = typename Kernel_traits<Line_2>::Kernel;
    using FT = typename Traits::FT;
    using Vector_2 = typename Traits::Vector_2;

    FT min_proj_value = max_value<FT>();
    FT max_proj_value = -max_value<FT>();

    const Vector_2 ref_vector = line.to_vector();
    const Point_2& ref_point = get(point_map_2, item_range[indices[0]]);
    
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const Point_2& point = get(point_map_2, item_range[indices[i]]);
      
      const Vector_2 curr_vector(ref_point, point);
      const FT value = CGAL::scalar_product(curr_vector, ref_vector);
      
      if (value < min_proj_value) {
        min_proj_value = value;
        p = point; }
      if (value > max_proj_value) {
        max_proj_value = value;
        q = point; }
    }
  }

  template<
  typename Item_range,
  typename Point_map_2,
  typename Line_2>
  typename Kernel_traits<Line_2>::Kernel::FT
  points_squared_length_2(
    const Item_range& item_range,
    const Point_map_2& point_map_2,
    const std::vector<std::size_t>& indices,
    const Line_2& line) {

    using Traits = typename Kernel_traits<Line_2>::Kernel;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    Point_2 p, q;
    boundary_points_on_line_2(item_range, point_map_2, indices, line, p, q);
    const FT squared_length = 
    CGAL::squared_distance(line.projection(p), line.projection(q));
    return squared_length;
  }

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_NUMBER_UTILS_H
