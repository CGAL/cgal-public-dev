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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_STRUCT_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_STRUCT_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <utility>
#include <unordered_map>

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
#include <CGAL/Levels_of_detail/internal/property_map.h>
#include <CGAL/Levels_of_detail/internal/parameters.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  struct Vertex_info {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;

    FT default_z = internal::max_value<FT>();
    FT z = default_z;
  };

  template<typename GeomTraits>
  struct Face_info {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;

    bool tagged = false;
    Urban_object_type urban_tag = Urban_object_type::GROUND;
    std::size_t object_index = std::size_t(-1);

    FT default_z = internal::max_value<FT>();
    std::vector<FT> z{default_z, default_z, default_z};
    bool interior = true;
  };

  template<typename GeomTraits>
  struct Triangulation {

    using Traits = GeomTraits;
    
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Plane_3 = typename Traits::Plane_3;
    using Triangle_2 = typename Traits::Triangle_2;

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
    output_for_lod0(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (auto fh = delaunay.finite_faces_begin(); 
      fh != delaunay.finite_faces_end(); ++fh) {
        if (!fh->info().interior) continue;
        
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

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_visibility(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Visibility_label visibility,
      const FT z) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (auto fh = delaunay.finite_faces_begin(); 
      fh != delaunay.finite_faces_end(); ++fh) {
        
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_2& q = fh->vertex(k)->point();
          const Point_3 p = Point_3(q.x(), q.y(), z);
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, visibility);
      }
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_lod12(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (auto fh = delaunay.finite_faces_begin(); 
      fh != delaunay.finite_faces_end(); ++fh) {
        if (fh->info().tagged) continue;
        
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

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_object(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t object_index) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (auto fh = delaunay.finite_faces_begin(); 
      fh != delaunay.finite_faces_end(); ++fh) {
        if (!fh->info().tagged) continue;
        
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3 p = get_point_3(fh->vertex(k));
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, object_index);
      }
      return std::make_pair(vertices, faces);
    }

    FT area() const {
      
      FT total_area = FT(0); Triangle_2 triangle;
      for (auto fh = delaunay.finite_faces_begin(); 
      fh != delaunay.finite_faces_end(); ++fh) {

        triangle = Triangle_2(
          fh->vertex(0)->point(),
          fh->vertex(1)->point(),
          fh->vertex(2)->point());
        total_area += triangle.area();
      }
      return total_area;
    }

    Point_3 locate(const Point_3& q) const {
      const Point_2 p = Point_2(q.x(), q.y());
      const auto fh = delaunay.locate(p);
      Point_3 p1, p2, p3;
      internal::point_3(fh, 0, p1);
      internal::point_3(fh, 1, p2);
      internal::point_3(fh, 2, p3);

      const Plane_3 plane = Plane_3(p1, p2, p3);
      const Point_3 res = internal::position_on_plane_3(p, plane);
      return res;
    }

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
    output_for_lod0(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      Indexer indexer;
      std::size_t num_vertices = 0;
      return triangulation.output_for_lod0(
        indexer, num_vertices, vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_lod0(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      
      return triangulation.output_for_lod0(
        indexer, num_vertices, vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_lod12(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      
      return triangulation.output_for_lod12(
        indexer, num_vertices, vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_object(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t object_index) const {
      
      return triangulation.output_for_object(
        indexer, num_vertices, vertices, faces, object_index);
    }
  };

  template<typename GeomTraits>
  struct Building_wall {

    using Traits = GeomTraits;
    using Point_3 = typename GeomTraits::Point_3;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    using Indexer = internal::Indexer<Point_3>;

    std::vector<Triangle_3> triangles;

    bool empty() const {
      return triangles.empty();
    }

    template<
    typename InputTriangulation,
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_lod(
      const InputTriangulation& tri,
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3); std::size_t tri_idx = 0;
      for (const auto& triangle : triangles) {
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3& q = triangle[k];

          Point_3 p;
          if (tri_idx % 2 == 0 && (k == 0 || k == 1))
            p = tri.locate(q);
          else if (tri_idx % 2 != 0 && k == 2)
            p = tri.locate(q);
          else
            p = q;

          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, Urban_object_type::BUILDING);
        ++tri_idx;
      }
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_object(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t object_index) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (const auto& triangle : triangles) {
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3& p = triangle[k];
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, object_index);
      }
      return std::make_pair(vertices, faces);
    }
  };

  template<typename GeomTraits>
  struct Building_roof {

    using Traits = GeomTraits;
    using Point_3 = typename GeomTraits::Point_3;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    using Indexer = internal::Indexer<Point_3>;
    
    std::vector<Triangle_3> triangles;

    bool empty() const {
      return triangles.empty();
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_lod(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (const auto& triangle : triangles) {
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3& p = triangle[k];
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, Urban_object_type::BUILDING);
      }
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_object(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t object_index) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (const auto& triangle : triangles) {
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3& p = triangle[k];
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, object_index);
      }
      return std::make_pair(vertices, faces);
    }
  };

  template<typename GeomTraits>
  struct Boundary {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;

    Segment_2 segment;
    
    FT default_z = internal::max_value<FT>();
    FT z = default_z;
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

    Base base0, base1, base2;
    std::vector<Wall> walls1, walls2;
    std::vector<Roof> roofs1, roofs2;
    std::vector<Edge> edges0, edges1, edges2;

    FT default_z = internal::max_value<FT>();
    FT bottom_z = default_z;
    FT top_z = default_z;

    std::size_t index = std::size_t(-1);
    Urban_object_type urban_tag = Urban_object_type::BUILDING;

    std::size_t cluster_index = std::size_t(-1);

    bool empty0() const {
      return base0.empty();
    }

    bool empty1() const {
      return base1.empty() && walls1.empty() && roofs1.empty();
    }

    bool empty2() const {
      return base2.empty() && walls2.empty() && roofs2.empty();
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod0(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      
      if (empty0())
        return boost::none;
        
      return base0.output_for_lod0(indexer, num_vertices, vertices, faces);
    }

    template<
    typename InputTriangulation,
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod1(
      const InputTriangulation& tri,
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const bool out_base = false) const {
      
      if (empty1())
        return boost::none;

      if (out_base) 
        base1.output_for_lod0(indexer, num_vertices, vertices, faces);
      for (const auto& wall : walls1)
        wall.output_for_lod(tri, indexer, num_vertices, vertices, faces);
      for (const auto& roof : roofs1)
        roof.output_for_lod(indexer, num_vertices, vertices, faces);
      return std::make_pair(vertices, faces);
    }

    template<
    typename InputTriangulation,
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod2(
      const InputTriangulation& tri,
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const bool out_base = false) const {
      
      if (empty2())
        return boost::none;

      if (out_base) 
        base2.output_for_lod0(indexer, num_vertices, vertices, faces);
      for (const auto& wall : walls2)
        wall.output_for_lod(tri, indexer, num_vertices, vertices, faces);
      for (const auto& roof : roofs2)
        roof.output_for_lod(indexer, num_vertices, vertices, faces);
      return std::make_pair(vertices, faces);
    }
  };

  template<typename GeomTraits>
  struct Tree_trunk {

    using Traits = GeomTraits;
    using Point_3 = typename GeomTraits::Point_3;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    using Indexer = internal::Indexer<Point_3>;

    std::vector<Triangle_3> triangles;

    bool empty() const {
      return triangles.empty();
    }

    template<
    typename InputTriangulation,
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_lod(
      const InputTriangulation& tri,
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3); std::size_t tri_idx = 0;
      for (const auto& triangle : triangles) {
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3& q = triangle[k];

          Point_3 p;
          if (tri_idx % 2 == 0 && (k == 0 || k == 1))
            p = tri.locate(q);
          else if (tri_idx % 2 != 0 && k == 2)
            p = tri.locate(q);
          else
            p = q;

          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, Urban_object_type::TREE);
        ++tri_idx;
      }
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_object(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t object_index) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (const auto& triangle : triangles) {
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3& p = triangle[k];
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, object_index);
      }
      return std::make_pair(vertices, faces);
    }
  };

  template<typename GeomTraits>
  struct Tree_crown {

    using Traits = GeomTraits;
    using Point_3 = typename GeomTraits::Point_3;
    using Triangle_3 = typename GeomTraits::Triangle_3;
    using Indexer = internal::Indexer<Point_3>;

    std::vector<Triangle_3> triangles;

    bool empty() const {
      return triangles.empty();
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_lod(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (const auto& triangle : triangles) {
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3& p = triangle[k];
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, Urban_object_type::TREE);
      }
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_object(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t object_index) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      for (const auto& triangle : triangles) {
        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3& p = triangle[k];
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, object_index);
      }
      return std::make_pair(vertices, faces);
    }
  };

  template<typename GeomTraits>
  struct Tree_model {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    Tree_model(
      const Point_2 center_, 
      const FT radius_, 
      const std::size_t cluster_index_) :
    center(center_), radius(radius_), 
    cluster_index(cluster_index_)
    { }

    Point_2 center = Point_2(FT(0), FT(0));
    FT radius = FT(0);
    std::size_t cluster_index = std::size_t(-1);
    std::size_t index = std::size_t(-1);

    FT default_z = internal::max_value<FT>();
    std::vector<FT> crown_z{default_z, default_z, default_z, default_z, default_z};
    std::vector<FT> crown_r{-FT(1), -FT(1), -FT(1), -FT(1)};

    FT trunk2_radius() const {
      return radius / FT(2);
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
    using Edge = Boundary<Traits>;

    using Indexer = internal::Indexer<Point_3>;

    Base base0, base1, base2;
    Trunk trunk1, trunk2;
    Crown crown1, crown2;
    std::vector<Edge> edges0, edges1, edges2;

    FT default_z = internal::max_value<FT>();
    FT bottom_z = default_z;
    FT top_z = default_z;

    std::size_t index = std::size_t(-1);
    Urban_object_type urban_tag = Urban_object_type::TREE;

    bool empty0() const {
      return base0.empty();
    }

    bool empty1() const {
      return base1.empty() && trunk1.empty() && crown1.empty();
    }

    bool empty2() const {
      return base2.empty() && trunk2.empty() && crown2.empty();
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod0(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      
      if (empty0())
        return boost::none;

      return base0.output_for_lod0(indexer, num_vertices, vertices, faces);
    }

    template<
    typename InputTriangulation,
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod1(
      const InputTriangulation& tri,
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const bool out_base = false) const {
      
      if (empty1())
        return boost::none;

      if (out_base) 
        base1.output_for_lod0(indexer, num_vertices, vertices, faces);
      trunk1.output_for_lod(tri, indexer, num_vertices, vertices, faces);
      crown1.output_for_lod(indexer, num_vertices, vertices, faces);
      return std::make_pair(vertices, faces);
    }

    template<
    typename InputTriangulation,
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod2(
      const InputTriangulation& tri,
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const bool out_base = false) const {
      
      if (empty2())
        return boost::none;

      if (out_base) 
        base2.output_for_lod0(indexer, num_vertices, vertices, faces);
      trunk2.output_for_lod(tri, indexer, num_vertices, vertices, faces);
      crown2.output_for_lod(indexer, num_vertices, vertices, faces);
      return std::make_pair(vertices, faces);
    }
  };

  template<typename GeomTraits>
  struct Partition_face_2 {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Triangulation = Triangulation<Traits>;
    using Vertex_handle = typename Triangulation::Delaunay::Vertex_handle;
    using Indexer = internal::Indexer<Point_3>;

    Triangulation base;
    std::vector<int> neighbors;
    Visibility_label visibility = Visibility_label::OUTSIDE;
    bool exterior = false;
    FT inside = FT(0);
    FT outside = FT(1);
    FT weight = FT(0);
    std::vector<Segment_2> edges;
    std::unordered_map<int, bool> constraints;

    void get_neighbors(std::vector<std::size_t>& out) const {
      out.clear();
      for (const int idx : neighbors)
        out.push_back(std::size_t(idx));
    }

    void compute_weight() {
      weight = base.area();
    }

    Partition_face_2() { }
    Partition_face_2(
      const std::vector<Point_2>& polygon) {

      neighbors.clear();
      base.delaunay.clear();
      std::vector<Vertex_handle> vhs;
      vhs.reserve(polygon.size());
      for (const Point_2& p : polygon)
        vhs.push_back(base.delaunay.insert(p));

      for (std::size_t i = 0; i < vhs.size(); ++i) {
        const std::size_t ip = (i + 1) % vhs.size();
        if (vhs[i] != vhs[ip])
          base.delaunay.insert_constraint(vhs[i], vhs[ip]);
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_visibility(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const FT z) const {

      return base.output_for_visibility(
        indexer, num_vertices, vertices, faces, visibility, z);
    }
  };

  template<typename GeomTraits>
  struct Partition_edge_2 {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    
    std::pair<int, int> neighbors;
    Segment_2 segment;
    FT weight = FT(0);

    void compute_weight() {
      weight = internal::distance(segment.source(), segment.target());
    }

    Partition_edge_2(
      const Point_2 p, const Point_2 q, 
      const int fidxi, const int fidxj) {

      segment = Segment_2(p, q);
      neighbors.first = fidxi;
      neighbors.second = fidxj;
    }
  };

  template<typename GeomTraits>
  struct Partition_2 {

    using Traits = GeomTraits;
    using Face = Partition_face_2<Traits>;
    using Edge = Partition_edge_2<Traits>;

    std::vector<Face> faces;
    std::vector<Edge> edges;

    void clear() {
      faces.clear();
      edges.clear();
    }

    const bool empty() const {
      return faces.empty() || edges.empty();
    }
  };

  template<typename GeomTraits>
  struct Partition_face_3 {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    
    Visibility_label visibility = Visibility_label::OUTSIDE;
    bool exterior = false;
    FT inside = FT(0);
    FT outside = FT(1);
    FT weight = FT(0);

    void compute_weight() {
      
    }
  };

  template<typename GeomTraits>
  struct Partition_edge_3 {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Triangle_3 = typename Traits::Triangle_3;
    using Indexer = internal::Indexer<Point_3>;

    std::pair<int, int> neighbors;
    FT weight = FT(0);

    std::vector<Point_3> polygon;

    const bool empty() const {
      return polygon.empty();
    }

    void compute_weight() {
      
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_for_object(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t object_index) const {

      if (empty())
        return boost::none;

      std::vector<std::size_t> face(3);
      const Point_3& ref = polygon[0];
      for (std::size_t i = 1; i < polygon.size(); ++i) {
        const std::size_t ip = (i + 1) % polygon.size();

        const Point_3& p1 = ref;
        const Point_3& p2 = polygon[i];
        const Point_3& p3 = polygon[ip];
        const Triangle_3 triangle = Triangle_3(p1, p2, p3);

        for (std::size_t k = 0; k < 3; ++k) {
          const Point_3& p = triangle[k];
          const std::size_t idx = indexer(p);
          if (idx == num_vertices) {
            *(vertices++) = p; 
            ++num_vertices;
          }
          face[k] = idx;
        }
        *(faces++) = std::make_pair(face, object_index);
      }
      return std::make_pair(vertices, faces);
    }
  };

  template<typename GeomTraits>
  struct Partition_3 {

    using Traits = GeomTraits;
    using Face = Partition_face_3<Traits>;
    using Edge = Partition_edge_3<Traits>;

    std::vector<Face> faces;
    std::vector<Edge> edges;

    void clear() {
      faces.clear();
      edges.clear();
    }

    const bool empty() const {
      return faces.empty() || edges.empty();
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
    using Point_2 = typename Traits::Point_2;

    using Point_map_3 = 
    internal::Item_property_map<Input_range, Point_map>;
    using Point_map_3_to_2 = 
    internal::Point_2_from_point_3_property_map<Point_map, Point_2>;
    using Point_map_2 =
    internal::Item_property_map<Input_range, Point_map_3_to_2>;
    using Visibility_map_d =
    internal::Item_property_map<Input_range, Visibility_map>;

    const Input_range& input_range;
    const Point_map& point_map;
    const Semantic_map& semantic_map;
    const Visibility_map& visibility_map;
    const bool verbose;

    Point_map_3 point_map_3;
    Point_map_3_to_2 point_map_3_to_2;
    Point_map_2 point_map_2;
    Visibility_map_d visibility_map_d;

    internal::Parameters<FT> parameters;

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
    verbose(verbose_),
    point_map_3(input_range, point_map),
    point_map_3_to_2(point_map),
    point_map_2(input_range, point_map_3_to_2),
    visibility_map_d(input_range, visibility_map, false) 
    { }

    ~Data_structure() 
    { }

    void points(
      const Semantic_label output_label,
      std::vector<std::size_t>& indices) const {

      indices.clear();
      for (std::size_t i = 0; i < input_range.size(); ++i) {
        const Semantic_label label = 
        get(semantic_map, *(input_range.begin() + i));
        if (label == output_label)
          indices.push_back(i);
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_STRUCT_H
