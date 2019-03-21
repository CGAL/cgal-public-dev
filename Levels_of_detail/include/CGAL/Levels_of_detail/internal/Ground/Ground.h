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

// STL includes.
#include <vector>
#include <utility>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/property_map.h>
#include <CGAL/Levels_of_detail/internal/Spacial_search/K_neighbor_query.h>
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
    using Input_range = typename Data_structure::Input_range;
    using Point_map = typename Data_structure::Point_map;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    using Ground_base = internal::Ground_base<Traits>;
    
    using Point_map_3 = 
    internal::Item_property_map<Input_range, Point_map>;
    using Point_map_3_to_2 = 
    internal::Point_2_from_point_3_property_map<Point_map, Point_2>;
    using Point_map_2 =
    internal::Item_property_map<Input_range, Point_map_3_to_2>;

    using Ground_points = std::vector<std::size_t>;

    using Neighbor_query = 
    K_neighbor_query<Traits, Ground_points, Point_map_2>;
    using Planar_ground_builder =
    internal::Planar_ground_builder<Ground_base>;
    using Smooth_ground_builder = 
    internal::Smooth_ground_builder<Ground_base, 
    Ground_points, Point_map_3, Neighbor_query>;

    Ground(const Data_structure& data) : 
    m_data(data),
    m_point_map_3(m_data.input_range, m_data.point_map),
    m_point_map_3_to_2(m_data.point_map),
    m_point_map_2(m_data.input_range, m_point_map_3_to_2) { 
      m_data.ground_points(m_ground_points);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    std::pair<VerticesOutputIterator, FacesOutputIterator> 
    output(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      Reconstruction_type ground_type) const {

      switch (ground_type) {
        case Reconstruction_type::PLANAR_GROUND: {
          const auto& tri = m_planar_ground.triangulation;
          return tri.output(vertices, faces);
        }
        case Reconstruction_type::SMOOTH_GROUND: {
          const auto& tri = m_smooth_ground.triangulation;
          return tri.output(vertices, faces);
        }
        default:
          return std::make_pair(vertices, faces);
      }
    }

    void make_planar() {

      if (m_data.verbose) 
        std::cout << std::endl << "- Computing planar ground" << std::endl;
      
      if (m_data.verbose) 
        std::cout << "* initializing" << std::endl;
      initialize(m_planar_ground);

      if (m_data.verbose) 
        std::cout << "* creating triangulation" << std::endl;
      Planar_ground_builder builder(m_planar_ground);
      builder.initialize();
      builder.finilize();
    }

    void make_smooth(const FT ground_precision) {
      
      if (m_data.verbose) 
        std::cout << std::endl << "- Computing smooth ground" << std::endl;

      if (m_data.verbose) 
        std::cout << "* initializing" << std::endl;
      initialize(m_smooth_ground);

      if (m_data.verbose) 
        std::cout << "* creating triangulation" << std::endl;
      Neighbor_query neighbor_query(
        m_ground_points, 6, m_point_map_2);
      Smooth_ground_builder builder(
        m_smooth_ground, 
        m_ground_points, m_point_map_3,
        neighbor_query,
        ground_precision);
      builder.initialize();
      builder.finilize();
    }

  private:
    const Data_structure& m_data;
    
    Point_map_3 m_point_map_3;
    Point_map_3_to_2 m_point_map_3_to_2;
    Point_map_2 m_point_map_2;

    Ground_points m_ground_points;

    Ground_base m_planar_ground;
    Ground_base m_smooth_ground;

    void initialize(Ground_base& ground_base) const {
      internal::plane_from_points_3(
        m_ground_points, 
        m_point_map_3, 
        ground_base.plane);
      internal::bounding_box_2(
        m_ground_points, 
        m_point_map_2, 
        ground_base.bbox);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_GROUND_H
