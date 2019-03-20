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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_SPACIAL_KNN_SEARCH_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_SPACIAL_KNN_SEARCH_2_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/property_map.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename ValueType,
  typename PointMap>
  class Knn_search_2 {
    
  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;

    using Neighbors = std::vector<ValueType>;
    using Search_traits_2 = CGAL::Search_traits_2<Traits>;
              
    using Search_traits = CGAL::Search_traits_adapter<ValueType, PointMap, Search_traits_2>;
    using Search_tree   = CGAL::Kd_tree<Search_traits>;

    using Knn = CGAL::Orthogonal_k_neighbor_search<Search_traits>;
    using Distance = typename Knn::Distance;
    using Splitter = typename Knn::Splitter;

    template<typename InputRange>
    Knn_search_2(
      const InputRange& input, 
      PointMap point_map, 
      const std::size_t nb_neighbors) :
    m_tree(NULL),
    m_tree_point_map(point_map),
    m_nb_neighbors(nb_neighbors) { 

      create_tree_2(input);
    }

    ~Knn_search_2() {
      if (m_tree != NULL)
        delete m_tree;
    }

    void get_neighbors(
      const Point_2& query, 
      Neighbors& neighbors) const {
      
      neighbors.clear();
      Distance distance(m_tree_point_map);
      Knn search(*m_tree, query, m_nb_neighbors, 0, true, distance);
      for (typename Knn::iterator it = search.begin(); it != search.end(); ++ it)
        neighbors.push_back (it->first);
    }

    inline const PointMap& point_map() const {
      return m_tree_point_map;
    }

  private:
    Search_tree *m_tree;
    PointMap m_tree_point_map;
    const std::size_t m_nb_neighbors;

    template<typename InputRange>
    void create_tree_2(const InputRange& input) {

      std::vector<typename InputRange::const_iterator> input_iterators;
      input_iterators.reserve(input.size());
      for (typename InputRange::const_iterator it = input.begin();
      it != input.end(); ++it)
        input_iterators.push_back(it);
        
      m_tree = new Search_tree(
        input_iterators.begin(),
        input_iterators.end(),
        Splitter(),
        Search_traits(m_tree_point_map));
    }

  }; // Knn_search_2
  
} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_SPACIAL_KNN_SEARCH_2_H
