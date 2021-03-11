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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_H

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

// Clustering.
#include <CGAL/Levels_of_detail/internal/Clustering/Connected_components.h>

// Trees.
#include <CGAL/Levels_of_detail/internal/Trees/Trees_site.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Trees {

  public:
    using Data_structure = DataStructure;

    using Traits = typename Data_structure::Traits;
    using Point_map_3 = typename Data_structure::Point_map_3;

    using Point_3 = typename Traits::Point_3;

    using Tree = internal::Tree<Traits>;
    using Tree_ptr = std::shared_ptr<Tree>;
    using Vegetation_points = std::vector<std::size_t>;

    using Clustering =
    internal::Connected_components<Traits, Vegetation_points, Point_map_3>;
    using Construction_site =
    internal::Trees_site<Data_structure>;

    using Indexer = internal::Indexer<Point_3>;

    Trees(const Data_structure& data) :
    m_data(data) {
      m_data.points(Semantic_label::VEGETATION, m_vegetation_points);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_trees(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type lod_type) const {
      if (empty())
        return boost::none;
      switch (lod_type) {
        case Reconstruction_type::TREES0: {
          return lod0(vertices, faces); }
        case Reconstruction_type::TREES1: {
          return lod1(vertices, faces); }
        case Reconstruction_type::TREES2: {
          return lod2(vertices, faces); }
        default: {
          return boost::none; }
      }
    }

    void make_trees() {
      initialize();
      compute_footprints();
      extrude_footprints();
      compute_crowns();
    }

    void initialize() {
      if (empty())
        return;
      if (m_data.verbose)
        std::cout << std::endl << "- Initializing trees" << std::endl;
      create_clusters();
      create_construction_sites();
    }

    void compute_footprints() {
      if (empty())
        return;
      if (m_data.verbose)
        std::cout << std::endl << "- Computing tree footprints" << std::endl;
      for (auto& site : m_sites)
        site.compute_footprints();
    }

    void extrude_footprints() {
      if (empty())
        return;
      if (m_data.verbose)
        std::cout << std::endl << "- Extruding tree footprints" << std::endl;
      for (auto& site : m_sites)
        site.extrude_footprints();
    }

    void compute_crowns() {
      if (empty())
        return;
      if (m_data.verbose)
        std::cout << std::endl << "- Computing tree crowns" << std::endl;
      for (auto& site : m_sites)
        site.compute_crowns();
    }

    void get_trees(std::vector<Tree_ptr>& trees) const {
      trees.clear();
      std::vector<Tree_ptr> ptrs;
      for (const auto& site : m_sites) {
        site.get_trees(ptrs);
        for (const auto& ptr : ptrs)
          trees.push_back(ptr);
      }
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_tree_clusters(OutputIterator output) const {

      for (std::size_t i = 0; i < m_clusters.size(); ++i)
        for (const std::size_t idx : m_clusters[i])
          *(output++) = std::make_pair(get(m_data.point_map_3, idx), i);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_tree_points(
      OutputIterator output) const {

      std::size_t tree_index = 0;
      for (const auto& site : m_sites)
        site.get_tree_points(output, tree_index);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    get_tree_boundaries(
      OutputIterator output) const {

      std::size_t tree_index = 0;
      for (const auto& site : m_sites)
        site.get_tree_boundaries(output, tree_index);
      return output;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_tree_footprints(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      Indexer indexer; std::size_t tree_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_tree_footprints(
          indexer, num_vertices, vertices, faces, tree_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_extruded_tree_boundaries(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      Indexer indexer; std::size_t tree_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_extruded_tree_boundaries(
          indexer, num_vertices, vertices, faces, tree_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_extruded_tree_footprints(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      Indexer indexer; std::size_t tree_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_extruded_tree_footprints(
          indexer, num_vertices, vertices, faces, tree_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_tree_trunks(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      Indexer indexer; std::size_t tree_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_tree_trunks(
          indexer, num_vertices, vertices, faces, tree_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    get_tree_crowns(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      Indexer indexer; std::size_t tree_index = 0;
      std::size_t num_vertices = 0;
      for (const auto& site : m_sites)
        site.get_tree_crowns(
          indexer, num_vertices, vertices, faces, tree_index);
      return std::make_pair(vertices, faces);
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire0(
      OutputIterator output) const {

      std::vector<Tree_ptr> trees;
      get_trees(trees);
      if (trees.empty())
        return boost::none;

      for (const auto& tree : trees)
        tree->output_lod0_wire(output);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire1(
      OutputIterator output) const {

      std::vector<Tree_ptr> trees;
      get_trees(trees);
      if (trees.empty())
        return boost::none;

      for (const auto& tree : trees)
        tree->output_lod1_wire(tree->base1.triangulation, output);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire2(
      OutputIterator output) const {

      std::vector<Tree_ptr> trees;
      get_trees(trees);
      if (trees.empty())
        return boost::none;

      for (const auto& tree : trees)
        tree->output_lod2_wire(tree->base2.triangulation, output);
      return output;
    }

    bool empty() const {
      return m_vegetation_points.empty();
    }

  private:
    const Data_structure& m_data;
    Vegetation_points m_vegetation_points;

    std::vector<Vegetation_points> m_clusters;
    std::vector<Construction_site> m_sites;

    void create_clusters() {
      if (m_data.verbose)
        std::cout << "* clustering (trees)" << std::endl;

      Clustering clustering(
        m_vegetation_points, m_data.point_map_3,
        m_data.parameters.trees.cluster_scale,
        m_data.parameters.trees.min_cluster_size);

      clustering.create_clusters(m_clusters);
      CGAL_assertion(!m_clusters.empty());
    }

    void create_construction_sites() {
      if (m_data.verbose)
        std::cout << "* creating construction sites (trees)" << std::endl;
      m_sites.clear();

      m_sites.reserve(m_clusters.size());
      for (std::size_t i = 0; i < m_clusters.size(); ++i)
        m_sites.push_back(Construction_site(m_data, m_clusters[i], i));
      CGAL_assertion(m_sites.size() == m_clusters.size());
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
      std::vector<Tree_ptr> trees;
      get_trees(trees);

      if (trees.empty())
        return boost::none;

      for (const auto& tree : trees)
        tree->output_lod0(indexer, num_vertices, vertices, faces);
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
      std::vector<Tree_ptr> trees;
      get_trees(trees);

      if (trees.empty())
        return boost::none;

      for (const auto& tree : trees)
        tree->output_lod1(tree->base1.triangulation,
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
      std::vector<Tree_ptr> trees;
      get_trees(trees);

      if (trees.empty())
        return boost::none;

      for (const auto& tree : trees)
        tree->output_lod2(tree->base2.triangulation,
        indexer, num_vertices, vertices, faces, true);
      return std::make_pair(vertices, faces);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_H
