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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_SITE_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_SITE_H

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
#include <CGAL/Levels_of_detail/internal/Clustering/Vegetation_clustering.h>
#include <CGAL/Levels_of_detail/internal/Trees/Tree_model_estimator.h>
#include <CGAL/Levels_of_detail/internal/Trees/Tree_builder.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Trees_site {

  public:
    using Data_structure = DataStructure;

    using Traits = typename Data_structure::Traits;
    using Point_map = typename Data_structure::Point_map;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_3 = typename Traits::Segment_3;

    using Points = std::vector<std::size_t>;
    using Point_map_3 = typename Data_structure::Point_map_3;
    using Vegetation_clustering =
    internal::Vegetation_clustering<Traits, Points, Point_map_3>;
    using Iterator = Points::const_iterator;

    using Tree_model_estimator =
    internal::Tree_model_estimator<Traits, Points, Point_map_3>;
    using Tree_model = internal::Tree_model<Traits>;
    using Tree = internal::Tree<Traits>;
    using Tree_ptr = std::shared_ptr<Tree>;
    using Tree_builder = internal::Tree_builder<Traits, Points, Point_map_3>;
    using Indexer = internal::Indexer<Point_3>;

    Trees_site(
      const Data_structure& data,
      const Points& points,
      const std::size_t site_index) :
    m_data(data),
    m_points(points),
    m_site_index(site_index),
    m_footprints_computed(false),
    m_footprints_extruded(false),
    m_crowns_computed(false) {
      CGAL_precondition(m_points.size() > 0);
    }

    void compute_footprints() {
      cluster_points(
        m_data.parameters.trees.grid_cell_width_2,
        m_data.parameters.trees.min_height);
      estimate_tree_models(
        m_data.parameters.trees.min_radius_2);
      initialize_trees();
      compute_tree_footprints(
        m_data.parameters.trees.min_faces_per_footprint);
    }

    void extrude_footprints() {
      extrude_tree_footprints(
        m_data.parameters.trees.extrusion_type);
    }

    void compute_crowns() {
      estimate_crown_icons();
      compute_tree_crowns();
    }

    void get_trees(
      std::vector<Tree_ptr>& trees) const {

      trees.clear();
      if (m_trees.empty()) return;
      for (const auto& tree : m_trees)
        trees.push_back(std::make_shared<Tree>(tree));
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_tree_points(
      OutputIterator output,
      std::size_t& tree_index) const {

      for (const auto& model : m_tree_models) {
        for (const auto& it : m_clusters[model.cluster_index])
          *(output++) = std::make_pair(get(m_data.point_map_3, *it), tree_index);
        ++tree_index;
      }
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_tree_boundaries(
      OutputIterator output,
      std::size_t& tree_index) const {

      if (m_trees.empty())
        return boost::none;

      for (const auto& tree : m_trees) {
        for (const auto& edge : tree.edges0) {
          const Point_2& s = edge.segment.source();
          const Point_2& t = edge.segment.target();
          const FT z = edge.z;
          *(output++) = std::make_pair(
            Segment_3(Point_3(s.x(), s.y(), z), Point_3(t.x(), t.y(), z)),
            tree_index);
        }
        ++tree_index;
      }
      return output;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_tree_footprints(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      std::size_t& tree_index) const {

      for (const auto& tree : m_trees) {
        tree.base0.output_for_object(
          indexer, num_vertices, vertices, faces, tree_index);
        ++tree_index;
      }
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_extruded_tree_boundaries(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      std::size_t& tree_index) const {

      for (const auto& tree : m_trees) {
        tree.trunk1.output_for_object(
          indexer, num_vertices, vertices, faces, tree_index);
        ++tree_index;
      }
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_extruded_tree_footprints(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      std::size_t& tree_index) const {

      for (const auto& tree : m_trees) {
        tree.crown1.output_for_object(
          indexer, num_vertices, vertices, faces, tree_index);
        ++tree_index;
      }
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_tree_trunks(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      std::size_t& tree_index) const {

      for (const auto& tree : m_trees) {
        tree.trunk2.output_for_object(
          indexer, num_vertices, vertices, faces, tree_index);
        ++tree_index;
      }
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_tree_crowns(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      std::size_t& tree_index) const {

      for (const auto& tree : m_trees) {
        tree.crown2.output_for_object(
          indexer, num_vertices, vertices, faces, tree_index);
        ++tree_index;
      }
      return std::make_pair(vertices, faces);
    }

  private:
    const Data_structure& m_data;
    const Points& m_points;
    const std::size_t m_site_index;

    std::vector< std::vector<Iterator> > m_clusters;
    std::vector<Tree_model> m_tree_models;
    std::vector<Tree> m_trees;

    bool m_footprints_computed;
    bool m_footprints_extruded;
    bool m_crowns_computed;

    void cluster_points(
      const FT grid_cell_width_2,
      const FT min_height) {

      if (m_points.empty()) return;
      Vegetation_clustering clustering(m_points, m_data.point_map_3);
      clustering.detect(grid_cell_width_2, min_height, m_clusters);
    }

    void estimate_tree_models(
      const FT min_radius_2) {

      if (m_clusters.empty()) return;
      const Tree_model_estimator estimator(m_clusters, m_data.point_map_3);
      estimator.estimate_model_parameters(min_radius_2, m_tree_models);
      for (std::size_t i = 0; i < m_tree_models.size(); ++i)
        m_tree_models[i].index = i;
    }

    void initialize_trees() {

      if (m_tree_models.empty()) return;
      m_trees.clear();
      m_trees.resize(m_tree_models.size());
      for (std::size_t i = 0; i < m_trees.size(); ++i)
        m_trees[i].index = i;
    }

    void compute_tree_footprints(
      const std::size_t min_faces_per_footprint) {

      if (m_trees.empty()) return;
      CGAL_assertion(m_trees.size() == m_tree_models.size());
      const Tree_builder builder;
      for (std::size_t i = 0; i < m_tree_models.size(); ++i)
        builder.add_lod0(
          m_tree_models[i],
          m_clusters[m_tree_models[i].cluster_index],
          m_data.point_map_3,
          min_faces_per_footprint,
          m_trees[i]);
      m_footprints_computed = true;
    }

    void extrude_tree_footprints(
      const Extrusion_type extrusion_type) {

      if (!m_footprints_computed) return;
      CGAL_assertion(m_trees.size() == m_tree_models.size());
      const Tree_builder builder;
      for (std::size_t i = 0; i < m_tree_models.size(); ++i)
        builder.add_lod1(
          extrusion_type,
          m_clusters[m_tree_models[i].cluster_index],
          m_data.point_map_3,
          m_trees[i]);
      m_footprints_extruded = true;
    }

    void estimate_crown_icons() {

      if (!m_footprints_extruded) return;
      const Tree_model_estimator estimator(m_clusters, m_data.point_map_3);
      estimator.estimate_crown_parameters(m_tree_models);
    }

    void compute_tree_crowns() {

      if (!m_footprints_extruded) return;
      CGAL_assertion(m_trees.size() == m_tree_models.size());
      const Tree_builder builder;
      for (std::size_t i = 0; i < m_tree_models.size(); ++i)
        builder.add_lod2(m_tree_models[i], m_trees[i]);
      m_crowns_computed = true;
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_SITE_H
