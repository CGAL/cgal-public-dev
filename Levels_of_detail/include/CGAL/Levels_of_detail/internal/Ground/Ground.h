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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_GROUND_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_GROUND_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <memory>
#include <vector>
#include <utility>

// Boost includes.
#include <boost/optional/optional.hpp>
#include <boost/function_output_iterator.hpp>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>
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
    using Point_map_2 = typename Data_structure::Point_map_2;
    using Point_map_3 = typename Data_structure::Point_map_3;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Ground_base = internal::Ground_base<Traits>;
    using Ground_points = std::vector<std::size_t>;

    using Neighbor_query =
    internal::K_neighbor_query<Traits, Ground_points, Point_map_2>;
    using Planar_builder =
    internal::Planar_ground_builder<Ground_base>;
    using Smooth_builder =
    internal::Smooth_ground_builder<Ground_base,
    Ground_points, Point_map_3, Neighbor_query>;

    using Neighbor_query_ptr = std::shared_ptr<Neighbor_query>;
    using Planar_builder_ptr = std::shared_ptr<Planar_builder>;
    using Smooth_builder_ptr = std::shared_ptr<Smooth_builder>;

    Ground(const Data_structure& data) :
    m_data(data) {
      m_data.points(Semantic_label::GROUND, m_ground_points);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_ground(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type ground_type,
      const FT ground_precision) const {
      if (empty())
        return boost::none;

      switch (ground_type) {
        case Reconstruction_type::PLANAR_GROUND: {
          Ground_base planar_ground;
          make_planar_ground(planar_ground);
          return planar_ground.output_for_lod0(vertices, faces); }
        case Reconstruction_type::SMOOTH_GROUND: {
          Ground_base smooth_ground;
          make_smooth_ground(smooth_ground, ground_precision);
          return smooth_ground.output_for_lod0(vertices, faces);
        }
        default: {
          return boost::none; }
      }
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire0(
      OutputIterator output) const {

      Ground_base planar_ground;
      make_planar_ground(planar_ground);
      return planar_ground.output_all_edges(output);
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire12(
      OutputIterator output,
      const FT ground_precision) const {

      Ground_base smooth_ground;
      make_smooth_ground(smooth_ground, ground_precision);
      return smooth_ground.output_all_edges(output);
    }

    void make_planar_ground(
      Ground_base& planar_ground) const {
      if (empty())
        return;

      if (m_data.verbose)
        std::cout << std::endl << "- Computing planar ground" << std::endl;
      initialize(planar_ground);

      if (m_data.verbose)
        std::cout << "* creating triangulation" << std::endl;
      auto planar_builder_ptr = make_planar_builder(planar_ground);
      planar_builder_ptr->initialize();
      planar_builder_ptr->finilize();
    }

    void make_smooth_ground(
      Ground_base& smooth_ground,
      const FT ground_precision) const {
      if (empty())
        return;

      if (m_data.verbose)
        std::cout << std::endl << "- Computing smooth ground" << std::endl;
      initialize(smooth_ground);

      if (m_data.verbose)
        std::cout << "* creating triangulation" << std::endl;
      auto neighbor_query_ptr = make_neighbor_query();
      auto smooth_builder_ptr = make_smooth_builder(
        smooth_ground, *neighbor_query_ptr, ground_precision);
      smooth_builder_ptr->initialize();
      smooth_builder_ptr->finilize();
    }

    void initialize(Ground_base& ground_base) const {
      if (m_data.verbose)
        std::cout << "* initializing ground base" << std::endl;
      internal::plane_from_points_3(
        m_ground_points,
        m_data.point_map_3,
        ground_base.plane);
      internal::bounding_box_2(
        m_ground_points,
        m_data.point_map_2,
        ground_base.bbox);
    }

    Neighbor_query_ptr make_neighbor_query() const {
      return std::make_shared<Neighbor_query>(
        m_ground_points, 6, m_data.point_map_2);
    }

    Planar_builder_ptr make_planar_builder(
      Ground_base& ground_base) const {
      return std::make_shared<Planar_builder>(
        ground_base);
    }

    Smooth_builder_ptr make_smooth_builder(
      Ground_base& ground_base,
      Neighbor_query& neighbor_query,
      const FT ground_precision) const {
      return std::make_shared<Smooth_builder>(
        ground_base,
        m_ground_points, m_data.point_map_3,
        neighbor_query,
        ground_precision);
    }

    bool empty() const {
      return m_ground_points.empty();
    }

  private:
    const Data_structure& m_data;
    Ground_points m_ground_points;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_GROUND_H
