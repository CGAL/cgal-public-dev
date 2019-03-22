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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_STRUCT_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_STRUCT_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <utility>

// Boost includes.
#include <boost/optional/optional.hpp>

// CGAL includes.
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/number_utils.h>
#include <CGAL/Levels_of_detail/internal/parameters.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  struct Partition_face_2 {

    using Traits = GeomTraits;
    
    bool valid = true;
    Visibility_label visibility = Visibility_label::OUTSIDE;
  };

  template<typename GeomTraits>
  struct Partition_2 {

    using Traits = GeomTraits;
    using Face = Partition_face_2<Traits>;

    std::vector<Face> faces;
  };

  template<typename GeomTraits>
  struct Partition_facet_3 {

    using Traits = GeomTraits;
    
    bool valid = true;
    Visibility_label visibility = Visibility_label::OUTSIDE;
  };

  template<typename GeomTraits>
  struct Partition_3 {

    using Traits = GeomTraits;
    using Facet = Partition_facet_3<Traits>;

    std::vector<Facet> facets;
  };

  template<typename GeomTraits>
  struct Vertex_info {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;

    const FT default_z = internal::max_value<FT>();
    FT z = default_z;
  };

  template<typename GeomTraits>
  struct Face_info {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;

    bool tagged = false;
    Urban_object_type urban_tag = Urban_object_type::GROUND;
    std::size_t object_index = std::size_t(-1);

    const FT default_z = internal::max_value<FT>();
    std::vector<FT> z{default_z, default_z, default_z};
  };

  template<typename GeomTraits>
  struct Triangulation {

    using Traits = GeomTraits;
    
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;

    using VI = Vertex_info<Traits>;
    using FI = Face_info<Traits>;

    using VBI = CGAL::Triangulation_vertex_base_with_info_2<VI, Traits>;
    using FBI = CGAL::Triangulation_face_base_with_info_2<FI, Traits>;
    using CFB = CGAL::Constrained_triangulation_face_base_2<Traits, FBI>;
    using TAG = CGAL::Exact_predicates_tag;
    using TDS = CGAL::Triangulation_data_structure_2<VBI, CFB>;

    using Delaunay = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS, TAG>;
    using Vertex_handle = typename Delaunay::Vertex_handle;
    using Face_handle = typename Delaunay::Finite_faces_iterator;
    using Indexer = internal::Indexer<Point_3>;

    Delaunay delaunay;

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (auto fh = delaunay.finite_faces_begin(); 
      fh != delaunay.finite_faces_end(); ++fh) {
        
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3 p = get_point_3(fh->vertex(k));
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, fh->info().urban_tag);
      }
      return std::make_pair(vertices, faces);
    }

    void locate(const Point_2& query) const { }
    void locate(const Point_3& query) const { }

    bool empty() const {
      return ( delaunay.number_of_faces() == 0 );
    }

    Point_3 get_point_3(const Vertex_handle& vh) const {
      const Point_2& p = vh->point(); const FT z = get_z(vh);
      return Point_3(p.x(), p.y(), z);
    }

    FT get_z(const Vertex_handle& vh) const {
      auto fh = delaunay.incident_faces(vh);
      const auto start = fh;

      FT sum = FT(0); FT num_faces = FT(0);
      do { if (delaunay.is_infinite(fh)) {
          ++fh; continue; }
        sum += fh->info().z[fh->index(vh)];
        num_faces += FT(1); ++fh;
      } while (fh != start);

      CGAL_assertion(num_faces > FT(0));
      const FT z = sum / num_faces;
      return z;
    }
  };

  template<typename GeomTraits>
  struct Ground_base {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Plane_3 = typename Traits::Plane_3;
    using Triangulation = Triangulation<Traits>;
    using Indexer = internal::Indexer<Point_3>;

    Triangulation triangulation;
    Plane_3 plane = Plane_3(FT(0), FT(0), FT(1), FT(0));
    std::vector<Point_2> bbox{
      Point_2(FT(-1), FT(-1)),
      Point_2(FT(1), FT(-1)),
      Point_2(FT(1), FT(1)),
      Point_2(FT(-1), FT(1))
    };

    bool empty() const {
      return triangulation.empty();
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      Indexer indexer;
      std::size_t num_vertices = 0;
      return triangulation.output(indexer, num_vertices, vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      return triangulation.output(indexer, num_vertices, vertices, faces);
    }
  };

  template<typename GeomTraits>
  struct Building_wall {

    using Traits = GeomTraits;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    
    std::vector<Triangle_3> triangles;
  };

  template<typename GeomTraits>
  struct Building_roof {

    using Traits = GeomTraits;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    
    std::vector<Triangle_3> triangles;
  };

  template<typename GeomTraits>
  struct Boundary {

    using Traits = GeomTraits;
    using Segment_2 = typename Traits::Segment_2;

    Segment_2 segment;
  };

  template<typename GeomTraits>
  struct Building {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Base = Ground_base<Traits>;
    using Wall = Building_wall<Traits>;
    using Roof = Building_roof<Traits>;
    using Edge = Boundary<Traits>;

    using Indexer = internal::Indexer<Point_3>;

    Base base;
    std::vector<Wall> walls;
    std::vector<Roof> roofs;
    std::vector<Edge> boundaries;

    const FT default_z = internal::max_value<FT>();
    FT bottom_z = default_z;
    FT top_z = default_z;

    std::size_t index = std::size_t(-1);
    const Urban_object_type urban_tag = Urban_object_type::BUILDING;

    const std::vector<Edge>& boundaries0() const { 
      return boundaries;
    }
    const std::vector<Edge>& boundaries1() const { 
      return boundaries;
    }
    const std::vector<Edge>& boundaries2() const { 
      return boundaries;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod1(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      
      return boost::none;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod2(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      
      return boost::none;
    }
  };

  template<typename GeomTraits>
  struct Tree_trunk {

    using Traits = GeomTraits;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    
    std::vector<Triangle_3> triangles;
  };

  template<typename GeomTraits>
  struct Tree_crown {

    using Traits = GeomTraits;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    
    std::vector<Triangle_3> triangles;
  };

  template<typename GeomTraits>
  struct Tree {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Base = Ground_base<Traits>;
    using Trunk = Tree_trunk<Traits>;
    using Crown = Tree_crown<Traits>;
    using Edge = Boundary<Traits>;

    using Indexer = internal::Indexer<Point_3>;

    Base base;
    Trunk trunk;
    Crown crown;
    std::vector<Edge> boundaries;

    const FT default_z = internal::max_value<FT>();
    FT bottom_z = default_z;
    FT top_z = default_z;

    std::size_t index = std::size_t(-1);
    const Urban_object_type urban_tag = Urban_object_type::TREE;

    const std::vector<Edge>& boundaries0() const { 
      return boundaries;
    }
    const std::vector<Edge>& boundaries1() const { 
      return boundaries;
    }
    const std::vector<Edge>& boundaries2() const { 
      return boundaries;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod1(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      
      return boost::none;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod2(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      
      return boost::none;
    }
  };

  template<
  typename GeomTraits, 
  typename InputRange, 
  typename PointMap,
  typename SemanticMap, 
  typename VisibilityMap>
  struct Data_structure {

    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Semantic_map = SemanticMap;
    using Visibility_map = VisibilityMap;

    using FT = typename Traits::FT;
    internal::Parameters<FT> parameters;

    const Input_range& input_range;
    const Point_map& point_map;
    const Semantic_map& semantic_map;
    const Visibility_map& visibility_map;
    const bool verbose;

    Data_structure(
      const Input_range& input_range_, 
      const Point_map& point_map_,
      const Semantic_map& semantic_map_, 
      const Visibility_map& visibility_map_,
      const bool verbose_ = false) : 
    input_range(input_range_),
    point_map(point_map_),
    semantic_map(semantic_map_),
    visibility_map(visibility_map_),
    verbose(verbose_) 
    { }

    ~Data_structure() 
    { }

    void ground_points(std::vector<std::size_t>& indices) const {

      for (std::size_t i = 0; i < input_range.size(); ++i) {
        const Semantic_label label = 
        get(semantic_map, *(input_range.begin() + i));
        if (label == Semantic_label::GROUND)
          indices.push_back(i);
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_STRUCT_H
