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
#include <vector>
#include <utility>

// Boost includes.
#include <boost/optional/optional.hpp>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/Clustering/Connected_components.h>

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
    using Vegetation_points = std::vector<std::size_t>;

    using Clustering = 
    internal::Connected_components<Traits, Vegetation_points, Point_map_3>;

    Trees(const Data_structure& data) : 
    m_data(data) { 
      m_data.points(Semantic_label::VEGETATION, m_vegetation_points);
      create_clusters();
    }

    void compute_footprints() {
      if (empty())
        return;

      if (m_data.verbose) 
        std::cout << std::endl << "- Computing tree footprints" << std::endl;
    }

    boost::optional<const std::vector<Tree>&> trees() const {
      if (m_trees.empty())
        return boost::none;
      return m_trees;
    }

    bool empty() const {
      return m_vegetation_points.empty();
    }

  private:
    const Data_structure& m_data;
    Vegetation_points m_vegetation_points;
    
    std::vector<Tree> m_trees;
    std::vector< std::vector<std::size_t> > m_clusters;

    void create_clusters() {
      if (m_data.verbose) 
        std::cout << "* clustering" << std::endl;
      m_clusters.clear();
      
      Clustering clustering(
        m_vegetation_points, m_data.point_map_3, m_data.parameters.scale);
        
      clustering.create_clusters(std::back_inserter(m_clusters));
      CGAL_assertion(!m_clusters.empty());
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_TREES_H
