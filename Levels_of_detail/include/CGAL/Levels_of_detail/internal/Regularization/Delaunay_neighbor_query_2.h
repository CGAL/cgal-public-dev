// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_DELAUNEY_NEIGHBOR_QUERY_2
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_DELAUNEY_NEIGHBOR_QUERY_2

#include <CGAL/license/Levels_of_detail.h>

#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <vector>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
    typename GeomTraits,
    typename InputRange,
    typename SegmentMap>
  class Delaunay_neighbor_query_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    using Point = typename GeomTraits::Point_2;
    using Segment = typename GeomTraits::Segment_2;

    using VB = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, GeomTraits>;
    using DS = CGAL::Triangulation_data_structure_2<VB>;
    using DT = CGAL::Delaunay_triangulation_2<GeomTraits, DS>;

    using Vertex_circulator = typename DT::Vertex_circulator;
    using Indices_map = std::vector <std::vector <std::size_t>>;

    Delaunay_neighbor_query_2(
      InputRange& input_range,
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map) { }

    void operator()(const std::size_t query_index, std::vector<std::size_t> & neighbors) {
      neighbors.clear();
      if(m_map_of_neighbors.size() == 0)
        return;
      CGAL_precondition(query_index >= 0 && query_index < m_map_of_neighbors.size());
      neighbors = m_map_of_neighbors[query_index];
    }

    template<typename Range, typename IndexMap = CGAL::Identity_property_map<std::size_t>>
  	void add_group(const Range& group, const IndexMap index_map = IndexMap()) {
      std::vector<std::size_t> gr;
      for (const auto & item : group) {
        const std::size_t seg_index = get(index_map, item);
        gr.push_back(seg_index);
      }

      if (gr.size() > 1) {
        build_delaunay_triangulation(gr);
        find_neighbors();
      }
    }

    void clear() {
      m_map_of_neighbors.clear();
    }

  private:
    Input_range& m_input_range;
    const Segment_map m_segment_map;
    DT m_dt;
    Indices_map m_map_of_neighbors;

    void build_delaunay_triangulation(const std::vector<std::size_t> & v) {
      m_dt.clear();
      for (std::size_t i = 0; i < v.size(); ++i) {
        const Segment& seg = get(m_segment_map, *(m_input_range.begin() + v[i]));
        const Point& source = seg.source();
        const Point& target = seg.target();
        const Point middle_point = internal::middle_point_2(source, target);

        auto vh = m_dt.insert(middle_point);
        vh->info() = v[i];
      }
    }

    void find_neighbors() {
      if(m_map_of_neighbors.size() < m_input_range.size()) {
        m_map_of_neighbors.clear();
        m_map_of_neighbors.resize(m_input_range.size());
      }
      for (auto vit = m_dt.finite_vertices_begin(); vit != m_dt.finite_vertices_end(); ++vit) {
        Vertex_circulator vc(vit);
        do {
          if(!m_dt.is_infinite(vc)) {
            CGAL_precondition(vit->info() >= 0 && vit->info() < m_input_range.size()
              && vc->info() >= 0 && vc->info() < m_input_range.size());
            m_map_of_neighbors[vit->info()].push_back(vc->info());
          }
          --vc;
        } while (vc != m_dt.incident_vertices(vit));
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_DELAUNEY_NEIGHBOR_QUERY_2
