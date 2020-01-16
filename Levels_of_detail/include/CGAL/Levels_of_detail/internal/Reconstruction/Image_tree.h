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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_TREE_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_TREE_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <queue>
#include <vector>
#include <utility>
#include <memory>

// CGAL includes.
#include <CGAL/enum.h>
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename Vertex,
typename Edge,
typename Halfedge,
typename Face,
typename Edge_type>
class Image_tree {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Segment_2 = typename Traits::Segment_2;
  using Segment_3 = typename Traits::Segment_3;
  
  using Image = internal::Image<Traits>;
  using Point_type = typename Image::Point_type;

  using FT_pair = std::pair<FT, FT>;
  using Indices = std::vector<std::size_t>;
  using Size_pair = std::pair<std::size_t, std::size_t>;

  using Saver = Saver<Traits>;
  using Inserter = Polygon_inserter<Traits>;
  using Indexer = internal::Indexer<Point_3>;
  using Color = CGAL::Color;

  Image_tree(
    std::vector<Vertex>& vertices,
    std::vector<Edge>& edges,
    std::vector<Halfedge>& halfedges, 
    std::vector<Face>& faces) :
  m_vertices(vertices),
  m_edges(edges),
  m_halfedges(halfedges),
  m_faces(faces),
  m_pi(static_cast<FT>(CGAL_PI)) 
  { }

  void build() {  
    clear(); build_tree();
  }

  void cut(const std::size_t level) {
    const std::size_t ref_label = m_faces[level].label;
    for (auto& face : m_faces)
      face.label = ref_label;
  }

  void clear() {
    
  }

  void check_vertex_information() {

    std::cout.precision(30);
    std::cout << "num vertices: " << m_vertices.size() << std::endl;
    for (const auto& vertex : m_vertices) {
      if (
        vertex.type != Point_type::FREE && 
        vertex.type != Point_type::LINEAR) {
        
        std::cout << 
        int(vertex.type) << " : " << 
        vertex.bd_idx << " , " <<
        vertex.labels.size() << " , " << 
        vertex.hedges.size() << " , " <<
        vertex.neighbors.size() << " , " <<
        vertex.faces.size() << std::endl;
      }
    }
  }

  void check_edge_information() {
    
    std::cout << "num edges: " << m_edges.size() << std::endl;
    std::cout << "num halfedges: " << m_halfedges.size() << std::endl;
    for (const auto& edge : m_edges) {
      std::cout << 
        int(edge.type) << " : " << 
        edge.labels.first << " " << 
        edge.labels.second << std::endl;
    }
  }

  void check_halfedge_information() {
    
    std::cout << "num edges: " << m_edges.size() << std::endl;
    std::cout << "num halfedges: " << m_halfedges.size() << std::endl;
    for (const auto& he : m_halfedges) {
      std::cout <<
      he.from_vertex << " " << he.to_vertex << " : " <<  
      int(m_vertices[he.from_vertex].type) << " " <<
      int(m_vertices[he.to_vertex].type) << " : " <<
      he.index << " " << 
      he.edg_idx << " " << 
      he.next << " " << 
      he.opposite << std::endl;
    }
  }

  void check_face_information() {
    
    std::cout << "num faces: " << m_faces.size() << std::endl;
    for (const auto& face : m_faces) {
      std::cout <<
      face.index << " : " << 
      int(face.type) << " , " <<
      face.label << " , " <<
      face.probs.size() << " , " <<
      face.hedges.size() << " , " <<
      face.neighbors.size() << " , " <<
      face.area << " , " <<
      face.skip << std::endl;
    }
  }

private:
  std::vector<Vertex>& m_vertices;
  std::vector<Edge>& m_edges;
  std::vector<Halfedge>& m_halfedges;
  std::vector<Face>& m_faces;
  const FT m_pi;

  void build_tree() {

  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_TREE_H
