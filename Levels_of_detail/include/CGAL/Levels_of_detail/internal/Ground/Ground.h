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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_GROUND_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_GROUND_H

#include <CGAL/license/Levels_of_detail.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/property_map.h>
#include <CGAL/Levels_of_detail/internal/Spacial_search/Knn_search_2.h>
#include <CGAL/Levels_of_detail/internal/Ground/Planar_ground_builder.h>
#include <CGAL/Levels_of_detail/internal/Ground/Smooth_ground_builder.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Ground {

  public:
    using Data_structure = DataStructure;
    
    using Traits = typename Data_structure::Traits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    using Ground_base = internal::Ground_base<Traits>;

    using Filtered_range = typename Data_structure::Filtered_range;
    using Filter_iterator = typename Data_structure::Filter_iterator;
    using Point_map = typename Data_structure::Point_map;

    using Point_2_from_iterator_map =
    Point_2_from_iterator_property_map<Filter_iterator, Point_2, Point_map>;

    using Knn_search_2 = 
    Knn_search_2<Traits, Filter_iterator, Point_2_from_iterator_map>;

    /*
    using Planar_ground_builder =
    internal::Planar_ground_builder<Traits>;
    using Smooth_ground_builder = 
    internal::Smooth_ground_builder<Traits, Filtered_range, Point_map, Knn_search>; */

    Ground(const Data_structure& data) : 
    m_data(data)
    { }

    void make_planar() {

      if (m_data.verbose) 
        std::cout << std::endl << "- Computing planar ground" << std::endl;
      
      if (m_data.verbose) 
        std::cout << "* fitting plane" << std::endl;

      /*
      const auto& ground_points = m_data.ground_points();
      const auto& point_map = m_data.point_map;

      internal::plane_from_points_3(
        ground_points, 
        point_map, 
        m_planar_ground.plane);

      internal::bounding_box_2(
        ground_points, 
        point_map, 
        m_planar_ground.bbox);

      if (m_data.verbose) 
        std::cout << "* creating triangulation" << std::endl;

      Planar_ground_builder builder(m_planar_ground.plane);
      builder.initialize(m_planar_ground.bbox);
      builder.finilize(); */
    }

    void make_smooth(const FT ground_precision) {
      
      if (m_data.verbose) 
        std::cout << std::endl << "- Computing smooth ground" << std::endl;

      if (m_data.verbose) 
        std::cout << "* fitting plane" << std::endl;

      m_smooth_ground.plane = m_planar_ground.plane;

      if (m_data.verbose) 
        std::cout << "* creating triangulation" << std::endl;

      const auto& ground_points = m_data.ground_points();
      const auto& point_map = m_data.point_map;

      Knn_search_2 knn_search(
        ground_points,
        Point_2_from_iterator_map(point_map), 
        6);

      /*
      Smooth_ground_builder builder(
        m_planar_ground.plane,
        ground_points,
        point_map,
        knn_search,
        ground_precision
      );

      builder.initialize(m_planar_ground.bbox);
      builder.finilize(); */
    }

  private:
    const Data_structure& m_data;
    
    Ground_base m_planar_ground;
    Ground_base m_smooth_ground;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_GROUND_H
