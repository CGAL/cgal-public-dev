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
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/Clustering/Connected_components.h>
#include <CGAL/Levels_of_detail/internal/Trees/Construction_site.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Trees {

  public:
    using Data_structure = DataStructure;
    
    using Traits = typename Data_structure::Traits;
    using Point_map_3 = typename Data_structure::Point_map_3;

    using Tree = internal::Tree<Traits>;
    using Tree_ptr = std::shared_ptr<Tree>;
    using Vegetation_points = std::vector<std::size_t>;

    using Clustering = 
    internal::Connected_components<Traits, Vegetation_points, Point_map_3>;
    using Construction_site=
    internal::Construction_site<Data_structure>;

    Trees(const Data_structure& data) : 
    m_data(data) { 
      m_data.points(Semantic_label::VEGETATION, m_vegetation_points);
      
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

    void get_trees(std::vector<Tree_ptr>& trees) const {
      trees.clear();
      std::vector<Tree_ptr> ptrs;
      for (const auto& site : m_sites) {
        site.get_trees(ptrs);
        for (const auto& ptr : ptrs)
          trees.push_back(ptr);
      }
    }

    bool empty() const {
      return m_vegetation_points.empty();
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> get_tree_clusters(OutputIterator output) const {
      for (std::size_t i = 0; i < m_clusters.size(); ++i)
        for (const std::size_t idx : m_clusters[i])
          *(output++) = std::make_pair(get(m_data.point_map_3, idx), i);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> get_tree_points(
      OutputIterator output) const {
      
      std::size_t tree_idx = 0;
      for (const auto& site : m_sites)
        site.get_tree_points(output, tree_idx);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> get_tree_boundaries(
      OutputIterator output) const {
      return boost::none;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> > 
    get_tree_footprints(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      return boost::none;
    }

  private:
    const Data_structure& m_data;
    Vegetation_points m_vegetation_points;
    
    std::vector<Tree_ptr> m_trees;
    std::vector<Vegetation_points> m_clusters;
    std::vector<Construction_site> m_sites;

    void create_clusters() {
      if (m_data.verbose) 
        std::cout << "* clustering (trees)" << std::endl;
      
      Clustering clustering(
        m_vegetation_points, m_data.point_map_3, 
        m_data.parameters.trees.cluster_scale, 
        m_data.parameters.min_cluster_size);
        
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
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_H
