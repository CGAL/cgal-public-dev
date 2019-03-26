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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_CONSTRUCTION_SITE_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_CONSTRUCTION_SITE_H

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
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/Clustering/Vegetation_clustering.h>
#include <CGAL/Levels_of_detail/internal/Trees/Estimate_tree_models.h>
#include <CGAL/Levels_of_detail/internal/Trees/Tree_builder.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Construction_site {

  public:
    using Data_structure = DataStructure;

    using Traits = typename Data_structure::Traits;
    using Point_map = typename Data_structure::Point_map;

    using FT = typename Traits::FT;

    using Points = std::vector<std::size_t>;
    using Point_map_3 = typename Data_structure::Point_map_3;
    using Vegetation_clustering = 
    internal::Vegetation_clustering<Traits, Points, Point_map_3>;
    using Iterator = Points::const_iterator;
    
    using Estimator = 
    internal::Estimate_tree_models<Traits, Points, Point_map_3>;
    using Tree_model = internal::Tree_model<Traits>;
    using Tree = internal::Tree<Traits>;
    using Tree_ptr = std::shared_ptr<Tree>;
    using Tree_builder = internal::Tree_builder<Traits>;

    Construction_site(
      const Data_structure& data,
      const Points& points,
      const std::size_t site_index) : 
    m_data(data),
    m_points(points),
    m_site_index(site_index) { 
      CGAL_precondition(m_points.size() > 0);
    }

    void compute_footprints() {
      cluster_points(
        m_data.parameters.trees.grid_cell_width_2, 
        m_data.parameters.trees.min_height);
      estimate_tree_models(
        m_data.parameters.trees.min_radius_2);
      initialize_trees();
      create_tree_footprints(
        m_data.parameters.trees.min_faces_per_footprint);
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> get_tree_points(
      OutputIterator output,
      std::size_t& tree_idx) const {
      for (const auto& model : m_tree_models) {
        for (const auto& it : m_clusters[model.cluster_index])
          *(output++) = std::make_pair(get(m_data.point_map_3, *it), tree_idx);
        ++tree_idx;
      }
      return output;
    }

    void get_trees(std::vector<Tree_ptr>& trees) const {
      if (m_trees.empty()) return;
      trees.clear();
      for (const auto& tree : m_trees)
        trees.push_back(std::make_shared<Tree>(tree));
    }

  private:
    const Data_structure& m_data;
    const Points& m_points;
    const std::size_t m_site_index;

    std::vector< std::vector<Iterator> > m_clusters;
    std::vector<Tree_model> m_tree_models;
    std::vector<Tree> m_trees;

    void cluster_points(const FT grid_cell_width_2, const FT min_height) {
      if (m_points.empty()) return;
      Vegetation_clustering clustering(m_points, m_data.point_map_3);
      clustering.detect(grid_cell_width_2, min_height, m_clusters);
    }

    void estimate_tree_models(const FT min_radius_2) {
      if (m_clusters.empty()) return;
      const Estimator estimator(m_clusters, m_data.point_map_3);
      estimator.estimate(min_radius_2, m_tree_models);
    }

    void initialize_trees() {
      if (m_tree_models.empty()) return;
      m_trees.clear();
      m_trees.resize(m_tree_models.size());
    }

    void create_tree_footprints(const std::size_t min_faces_per_footprint) {
      if (m_trees.empty()) return;
      CGAL_assertion(m_trees.size() == m_tree_models.size());
      Tree_builder builder;
      for (std::size_t i = 0; i < m_tree_models.size(); ++i)
        builder.add_lod0(m_tree_models[i], min_faces_per_footprint, m_trees[i]);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_CONSTRUCTION_SITE_H
