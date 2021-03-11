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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDINGS_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDINGS_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <memory>
#include <vector>
#include <utility>

// Boost includes.
#include <boost/optional/optional.hpp>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/Clustering/Connected_components.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Buildings_site.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Buildings {

  public:
    using Data_structure = DataStructure;

    using Traits = typename Data_structure::Traits;
    using Point_map_3 = typename Data_structure::Point_map_3;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;

    using Building = internal::Building<Traits>;
    using Building_ptr = std::shared_ptr<Building>;
    using Building_points = std::vector<std::size_t>;

    using Data_range = std::vector< std::pair<std::size_t, bool> >;
    using Data_map = internal::From_pair_property_map<Point_map_3>;
    using Clustering =
    internal::Connected_components<Traits, Data_range, Data_map>;
    using Construction_site =
    internal::Buildings_site<Data_structure>;

    using Indexer = internal::Indexer<Point_3>;

    Buildings(const Data_structure& data) :
    m_data(data) {
      m_data.points(Semantic_label::BUILDING_INTERIOR, m_interior_points);
      m_data.points(Semantic_label::BUILDING_BOUNDARY, m_boundary_points);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_buildings(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type lod_type) const {
      if (empty())
        return boost::none;
      switch (lod_type) {
        case Reconstruction_type::BUILDINGS0: {
          return lod0(vertices, faces); }
        case Reconstruction_type::BUILDINGS1: {
          return lod1(vertices, faces); }
        case Reconstruction_type::BUILDINGS2: {
          return lod2(vertices, faces); }
        default: {
          return boost::none; }
      }
    }

    void make_buildings() {
      initialize();
      detect_boundaries();
      compute_footprints();
      extrude_footprints();
      detect_roofs();
      compute_roofs();
    }

    void initialize() {
      if (empty())
        return;
      if (m_data.verbose)
        std::cout << std::endl << "- Initializing buildings" << std::endl;
      create_clusters();
      create_construction_sites();
    }

    void detect_boundaries() {
      if (empty())
        return;
      if (m_data.verbose)
        std::cout << std::endl << "- Detecting building boundaries" << std::endl;

      std::size_t num_sites = m_sites.size();
      std::cout << "* sites processed: " << std::endl;
      std::size_t idx = 1;
      for (auto& site : m_sites) {
        site.detect_boundaries();
        std::cout << idx << " / " << num_sites << std::endl;
        ++idx;
      }
    }

    void compute_footprints() {
      if (empty())
        return;
      if (m_data.verbose)
        std::cout << std::endl << "- Computing building footprints" << std::endl;

      std::size_t num_sites = m_sites.size();
      std::cout << "* sites processed: " << std::endl;
      std::size_t idx = 1;
      for (auto& site : m_sites) {
        site.compute_footprints();
        std::cout << idx << " / " << num_sites << std::endl;
        ++idx;
      }
    }

    void extrude_footprints() {
      if (empty())
        return;
      if (m_data.verbose)
        std::cout << std::endl << "- Extruding building footprints" << std::endl;
      for (auto& site : m_sites)
        site.extrude_footprints();
    }

    void detect_roofs() {
      if (empty())
        return;
      if (m_data.verbose)
        std::cout << std::endl << "- Detecting building roofs" << std::endl;

      std::size_t num_sites = m_sites.size();
      std::cout << "* sites processed: " << std::endl;
      std::size_t idx = 1;
      for (auto& site : m_sites) {
        site.detect_roofs();
        std::cout << idx << " / " << num_sites << std::endl;
        ++idx;
      }
    }

    void compute_roofs() {
      if (empty())
        return;
      if (m_data.verbose)
        std::cout << std::endl << "- Computing building roofs" << std::endl;

      std::size_t num_sites = m_sites.size();
      std::cout << "* sites processed: " << std::endl;
      std::size_t idx = 1;
      for (auto& site : m_sites) {
        site.compute_roofs();
        std::cout << idx << " / " << num_sites << std::endl;
        ++idx;
      }
    }

    void get_buildings(std::vector<Building_ptr>& buildings) const {
      buildings.clear();
      std::vector<Building_ptr> ptrs;
      for (const auto& site : m_sites) {
        site.get_buildings(ptrs);
        for (const auto& ptr : ptrs)
          buildings.push_back(ptr);
      }
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_building_clusters(OutputIterator output) const {

      if (m_sites.empty())
        return boost::none;

      CGAL_assertion(
        m_interior_clusters.size() == m_boundary_clusters.size());

      const std::size_t n = m_interior_clusters.size();
      for (std::size_t i = 0; i < n; ++i) {
        const auto& cluster1 = m_interior_clusters[i];
        const auto& cluster2 = m_boundary_clusters[i];

        for (const std::size_t idx : cluster1)
          *(output++) = std::make_pair(get(m_data.point_map_3, idx), i);
        for (const std::size_t idx : cluster2)
          *(output++) = std::make_pair(get(m_data.point_map_3, idx), i);
      }
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_boundary_points(OutputIterator output) const {

      if (m_sites.empty())
        return boost::none;

      for (const auto& site : m_sites)
        site.get_boundary_points(output);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_wall_points(OutputIterator output) const {

      if (m_sites.empty())
        return boost::none;

      for (const auto& site : m_sites)
        site.get_wall_points(output);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_approximate_boundaries(
      OutputIterator output) const {

      if (m_sites.empty())
        return boost::none;

      for (const auto& site : m_sites)
        site.get_approximate_boundaries(output);
      return output;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_partitioning_2(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (m_sites.empty())
        return boost::none;

      Indexer indexer; FT offset = FT(0);
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites) {
        const auto& plane = site.ground_plane();
        const FT z = plane.point().z() + offset;
        site.get_partitioning_2(indexer, num_vertices, vertices, faces, z);
        offset += FT(1);
      }
      return std::make_pair(vertices, faces);
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_building_points(
      OutputIterator output) const {

      if (m_sites.empty())
        return boost::none;

      std::size_t building_index = 0;
      for (const auto& site : m_sites)
        site.get_building_points(output, building_index);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_building_boundaries(
      OutputIterator output) const {

      if (m_sites.empty())
        return boost::none;

      std::size_t building_index = 0;
      for (const auto& site : m_sites)
        site.get_building_boundaries(output, building_index);
      return output;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_building_footprints(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (m_sites.empty())
        return boost::none;

      Indexer indexer; std::size_t building_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_building_footprints(
          indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_extruded_building_boundaries(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (m_sites.empty())
        return boost::none;

      Indexer indexer; std::size_t building_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_extruded_building_boundaries(
          indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_extruded_building_footprints(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (m_sites.empty())
        return boost::none;

      Indexer indexer; std::size_t building_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_extruded_building_footprints(
          indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_roof_points(OutputIterator output) const {

      if (m_sites.empty())
        return boost::none;

      for (const auto& site : m_sites)
        site.get_roof_points(output);
      return output;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_building_approximate_bounds(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (m_sites.empty())
        return boost::none;

      Indexer indexer; std::size_t building_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_building_approximate_bounds(
          indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_building_partitioning_3(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (m_sites.empty())
        return boost::none;

      Indexer indexer; std::size_t building_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_building_partitioning_3(
          indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_building_walls(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (m_sites.empty())
        return boost::none;

      Indexer indexer; std::size_t building_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_building_walls(
          indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_building_roofs(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (m_sites.empty())
        return boost::none;

      Indexer indexer; std::size_t building_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_building_roofs(
          indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire0(
      OutputIterator output) const {

      std::vector<Building_ptr> buildings;
      get_buildings(buildings);
      if (buildings.empty())
        return boost::none;

      for (const auto& building : buildings)
        building->output_lod0_wire(output);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire1(
      OutputIterator output) const {

      std::vector<Building_ptr> buildings;
      get_buildings(buildings);
      if (buildings.empty())
        return boost::none;

      for (const auto& building : buildings)
        building->output_lod1_wire(building->base1.triangulation, output);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire2(
      OutputIterator output) const {

      std::vector<Building_ptr> buildings;
      get_buildings(buildings);
      if (buildings.empty())
        return boost::none;

      for (const auto& building : buildings)
        building->output_lod2_wire(building->base2.triangulation, output);
      return output;
    }

    bool empty() const {
      return m_interior_points.empty() && m_boundary_points.empty();
    }

  private:
    const Data_structure& m_data;
    Building_points m_interior_points;
    Building_points m_boundary_points;

    std::vector<Building_points> m_interior_clusters;
    std::vector<Building_points> m_boundary_clusters;
    std::vector<Building_points> m_exterior_clusters;
    std::vector<Construction_site> m_sites;

    void create_clusters() {
      if (m_data.verbose)
        std::cout << "* clustering (buildings)" << std::endl;

      // Merge interior and boundary points.
      Data_range data;
      data.reserve(m_interior_points.size() + m_boundary_points.size());
      for (const std::size_t idx : m_interior_points)
        data.push_back(std::make_pair(idx, true));
      for (const std::size_t idx : m_boundary_points)
        data.push_back(std::make_pair(idx, false));
      Data_map dmap(m_data.point_map_3);

      // Create clusters.
      std::vector<Data_range> clusters;
      Clustering clustering(
        data, dmap,
        m_data.parameters.buildings.cluster_scale,
        m_data.parameters.buildings.min_cluster_size);

      clustering.create_clusters(clusters);
      CGAL_assertion(!clusters.empty());

      // Split boundary and interior points.
      m_interior_clusters.clear();
      m_boundary_clusters.clear();
      m_exterior_clusters.clear();

      /* // Exterior points. Remove, if necessary.
      std::vector<Point_2> bbox; */

      Building_points interior, boundary, exterior;
      for (const auto& cluster : clusters) {

        interior.clear(); boundary.clear();
        for (const auto& item : cluster) {
          if (item.second) interior.push_back(item.first);
          else boundary.push_back(item.first);
        }

        if (interior.size() >= m_data.parameters.buildings.min_cluster_size &&
            boundary.size() >= 0) {

          m_interior_clusters.push_back(interior);
          m_boundary_clusters.push_back(boundary);

          /* // Exterior points. Remove, if necessary.
          internal::bounding_box_2(cluster, dmap, bbox);
          internal::scale_polygon_2(FT(6) / FT(5), bbox);
          const auto& minp = bbox[0];
          const auto& maxp = bbox[2];

          exterior.clear();
          for (std::size_t i = 0; i < m_data.input_range.size(); ++i) {
            const auto& p = get(m_data.point_map_2, *(m_data.input_range.begin() + i));
            if (p.x() > minp.x() && p.x() < maxp.x() &&
                p.y() > minp.y() && p.y() < maxp.y()) {

              if (
                (std::find(interior.begin(), interior.end(), i) == interior.end()) &&
                (std::find(boundary.begin(), boundary.end(), i) == boundary.end()) )
              exterior.push_back(i);
            }
          } */
          m_exterior_clusters.push_back(exterior);
        }
      }
    }

    void create_construction_sites() {
      if (m_data.verbose)
        std::cout << "* creating construction sites (buildings)" << std::endl;
      m_sites.clear();

      CGAL_assertion(
        m_interior_clusters.size() == m_boundary_clusters.size());
      CGAL_assertion(
        m_exterior_clusters.size() == m_interior_clusters.size());

      m_sites.reserve(m_interior_clusters.size());
      for (std::size_t i = 0; i < m_interior_clusters.size(); ++i)
        m_sites.push_back(Construction_site(
          m_data,
          m_interior_clusters[i],
          m_boundary_clusters[i],
          m_exterior_clusters[i],
          i));
      CGAL_assertion(m_sites.size() == m_interior_clusters.size());
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod0(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      Indexer indexer;
      std::size_t num_vertices = 0;
      std::vector<Building_ptr> buildings;
      get_buildings(buildings);

      if (buildings.empty())
        return boost::none;

      for (const auto& building : buildings)
        building->output_lod0(indexer, num_vertices, vertices, faces);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod1(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      Indexer indexer;
      std::size_t num_vertices = 0;
      std::vector<Building_ptr> buildings;
      get_buildings(buildings);

      if (buildings.empty())
        return boost::none;

      for (const auto& building : buildings)
        building->output_lod1(building->base1.triangulation,
        indexer, num_vertices, vertices, faces, true);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod2(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      Indexer indexer;
      std::size_t num_vertices = 0;
      std::vector<Building_ptr> buildings;
      get_buildings(buildings);

      if (buildings.empty())
        return boost::none;

      for (const auto& building : buildings)
        building->output_lod2(building->base2.triangulation,
        indexer, num_vertices, vertices, faces, true);
      return std::make_pair(vertices, faces);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDINGS_H
