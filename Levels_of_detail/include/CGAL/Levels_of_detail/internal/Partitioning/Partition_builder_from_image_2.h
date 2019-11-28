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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_H

// STL includes.
#include <map>
#include <set>
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Planar_image_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Linear_image_region.h>
#include <CGAL/Levels_of_detail/internal/Spatial_search/Image_neighbor_query.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename ImagePointer>
class Partition_builder_from_image_2 {

public:
  using Traits = GeomTraits;

  using FT          = typename Traits::FT;
  using Point_2     = typename Traits::Point_2;
  using Point_3     = typename Traits::Point_3;
  using Vector_2    = typename Traits::Vector_2;
  using Vector_3    = typename Traits::Vector_3;
  using Line_2      = typename Traits::Line_2;
  using Line_3      = typename Traits::Line_3;
  using Segment_2   = typename Traits::Segment_2;
  using Segment_3   = typename Traits::Segment_3;
  using Plane_3     = typename Traits::Plane_3;
  using Intersect_3 = typename Traits::Intersect_3;

  using Triangulation    = internal::Triangulation<Traits>;
  using Face_handle      = typename Triangulation::Delaunay::Face_handle;
  using Vertex_handle    = typename Triangulation::Delaunay::Vertex_handle;
  using Partition_edge_2 = internal::Partition_edge_2<Traits>;
  using Partition_face_2 = internal::Partition_face_2<Traits>;

  using Partition_2   = internal::Partition_2<Traits>;
  using LF_circulator = typename Triangulation::Delaunay::Line_face_circulator;
  using Size_pair     = std::pair<std::size_t, std::size_t>;
  using Indices       = std::vector<std::size_t>;

  using Saver = Saver<Traits>;

  struct Pixel {

  };

  Partition_builder_from_image_2(
    const std::vector<Segment_2>& boundary,
    ImagePointer& image_ptr,
    Partition_2& partition_2) :
  m_boundary(boundary),
  m_image_ptr(image_ptr),
  m_partition_2(partition_2),
  m_pi(static_cast<FT>(CGAL_PI)) { 

    m_partition_2.clear();
  }

  void build() {


    save_partition_2(
      "/Users/monet/Documents/lod/logs/buildings/tmp/partition_step_0", false);
  }

  void add_constraints() {

  }

  void compute_visibility() {

  }

  void label_faces() {

  }

  void optimize() {

  }

private:
  const std::vector<Segment_2>& m_boundary;
  ImagePointer& m_image_ptr;
  Partition_2& m_partition_2;
  const FT m_pi;

  std::vector<Pixel> m_image;
  Triangulation m_base;

  Saver m_saver;

  void save_partition_2(
    const std::string path,
    const bool with_roof_colors) {

    const FT z = FT(0);
    std::size_t num_vertices = 0;
    internal::Indexer<Point_3> indexer;

    std::vector<Point_3> vertices; 
    std::vector<Indices> faces; 
    std::vector<CGAL::Color> fcolors;

    Polygon_inserter<Traits> inserter(faces, fcolors);
    auto output_vertices = std::back_inserter(vertices);
    auto output_faces = boost::make_function_output_iterator(inserter);

    for (const auto& face : m_partition_2.faces) {
      if (!with_roof_colors) {
        face.output_for_visibility(
          indexer, num_vertices, output_vertices, output_faces, z);
      } else {
        face.output_with_label_color(
          indexer, num_vertices, output_vertices, output_faces, z);
      }
    }
    m_saver.export_polygon_soup(vertices, faces, fcolors, path);
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARTITION_BUILDER_FROM_IMAGE_2_H
