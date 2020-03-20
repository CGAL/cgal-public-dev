// Copyright (c) 2020 SARL GeometryFactory (France).
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

#ifndef CGAL_URBAN_AREA_PROCESSING_STRUCT_H
#define CGAL_URBAN_AREA_PROCESSING_STRUCT_H

// #include <CGAL/license/Urban_area_processing.h>

// STL includes.
#include <vector>
#include <utility>
#include <unordered_map>

// Boost includes.
#include <boost/optional/optional.hpp>

// CGAL includes.
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Internal includes.
#include <CGAL/Urban_area_processing/enum.h>
#include <CGAL/Urban_area_processing/parameters.h>

#include <CGAL/Urban_area_processing/internal/utils.h>
#include <CGAL/Urban_area_processing/internal/property_map.h>

namespace CGAL {
namespace Urban_area_processing {

  template<typename GeomTraits>
  struct Edge {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;

    Segment_2 segment;
    
    const FT default_z = internal::max_value<FT>();
    FT z = default_z;
  };

  template<typename GeomTraits>
  struct Vertex_info {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;

    const FT default_z = internal::max_value<FT>();
    FT z = default_z;
    std::size_t index = std::size_t(-1);
    std::size_t object_index = std::size_t(-1);
  };

  template<typename GeomTraits>
  struct Face_info {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;

    const FT default_z = internal::max_value<FT>();
    std::vector<FT> z{default_z, default_z, default_z};
    std::size_t index = std::size_t(-1);
    std::size_t object_index = std::size_t(-1);

    bool tagged = false;
    Urban_object_type urban_tag = Urban_object_type::GROUND;

    std::size_t label = std::size_t(-1);
    bool used = false;
  };

  template<typename GeomTraits>
  struct Triangulation {

    using Traits = GeomTraits;
    
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Triangle_2 = typename Traits::Triangle_2;
    using Triangle_3 = typename Traits::Triangle_3;

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
    using Location_type = typename Delaunay::Locate_type;

    Delaunay delaunay;

    bool empty() const {
      return delaunay.number_of_faces() == 0;
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
      Point_2(FT( 1), FT(-1)),
      Point_2(FT( 1), FT( 1)),
      Point_2(FT(-1), FT( 1))
    };

    bool empty() const {
      return triangulation.empty();
    }
  };

  template<typename GeomTraits>
  struct Building_wall {

    using Traits = GeomTraits;
    using Point_2 = typename GeomTraits::Point_2;
    using Point_3 = typename GeomTraits::Point_3;
    using Segment_3 = typename GeomTraits::Segment_3;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    using Indexer = internal::Indexer<Point_3>;

    std::vector<Segment_3> segments;
    std::vector<Triangle_3> triangles;

    bool empty() const {
      return triangles.empty();
    }
  };

  template<typename GeomTraits>
  struct Building_roof {

    using Traits = GeomTraits;
    using Point_3 = typename GeomTraits::Point_3;
    using Segment_3 = typename GeomTraits::Segment_3;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    using Indexer = internal::Indexer<Point_3>;
    
    std::vector<Segment_3> segments;
    std::vector<Triangle_3> triangles;
    
    bool empty() const {
      return triangles.empty();
    }
  };

  template<typename GeomTraits>
  struct Building {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;

    using Base = Ground_base<Traits>;
    using Wall = Building_wall<Traits>;
    using Roof = Building_roof<Traits>;
    using Edge = Edge<Traits>;

    using Indexer = internal::Indexer<Point_3>;

    Base base0, base1;
    std::vector<Wall> walls1;
    std::vector<Roof> roofs1;
    std::vector<Edge> edges0, edges1;

    const FT default_z = internal::max_value<FT>();
    FT bottom_z = default_z;
    FT top_z = default_z;

    std::size_t index = std::size_t(-1);
  };

  template<typename GeomTraits>
  struct Tree_trunk {

    using Traits = GeomTraits;
    using Point_2 = typename GeomTraits::Point_2;
    using Point_3 = typename GeomTraits::Point_3;
    using Segment_3 = typename GeomTraits::Segment_3;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    using Indexer = internal::Indexer<Point_3>;

    std::vector<Segment_3> segments;
    std::vector<Triangle_3> triangles;

    bool empty() const {
      return triangles.empty();
    }
  };

  template<typename GeomTraits>
  struct Tree_crown {

    using Traits = GeomTraits;
    using Point_3 = typename GeomTraits::Point_3;
    using Segment_3 = typename GeomTraits::Segment_3;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    using Indexer = internal::Indexer<Point_3>;

    std::vector<Segment_3> segments;
    std::vector<Triangle_3> triangles;
    
    bool empty() const {
      return triangles.empty();
    }
  };

  template<typename GeomTraits>
  struct Tree {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Base = Ground_base<Traits>;
    using Trunk = Tree_trunk<Traits>;
    using Crown = Tree_crown<Traits>;
    using Edge = Edge<Traits>;

    using Indexer = internal::Indexer<Point_3>;

    Base base0, base1;
    Trunk trunk1;
    Crown crown1;
    std::vector<Edge> edges0, edges1;

    const FT default_z = internal::max_value<FT>();
    FT bottom_z = default_z;
    FT top_z = default_z;

    std::size_t index = std::size_t(-1);
  };

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_STRUCT_H
