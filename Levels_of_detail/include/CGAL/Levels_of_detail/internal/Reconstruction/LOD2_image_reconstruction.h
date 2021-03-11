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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_LOD2_IMAGE_RECONSTRUCTION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_LOD2_IMAGE_RECONSTRUCTION_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <queue>
#include <vector>
#include <utility>
#include <memory>

// CGAL includes.
#include <CGAL/property_map.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image_creator.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image_data_structure.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image_tree.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/Sphere_neighbor_query.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename ImagePointer,
typename Point_map>
class LOD2_image_reconstruction {

public:
  using Traits = GeomTraits;
  using Image_ptr = ImagePointer;

  using FT = typename Traits::FT;
  using Segment_2 = typename Traits::Segment_2;
  using Segment_3 = typename Traits::Segment_3;
  using Point_3 = typename Traits::Point_3;
  using Plane_3 = typename Traits::Plane_3;
  using Triangle_3 = typename Traits::Triangle_3;

  using Triangulation = internal::Triangulation<Traits>;
  using Partition_2 = internal::Partition_2<Traits>;
  using Building = internal::Building<Traits>;
  using Wall = typename Building::Wall;
  using Roof = typename Building::Roof;
  using Triangle_set_3 = std::vector<Triangle_3>;

  using Image_creator = internal::Image_creator<Traits, Image_ptr>;
  using Image_data_structure = internal::Image_data_structure<Traits>;
  using Image_tree = internal::Image_tree<
    Traits,
    typename Image_data_structure::Vertex,
    typename Image_data_structure::Edge,
    typename Image_data_structure::Halfedge,
    typename Image_data_structure::Face,
    typename Image_data_structure::Edge_type,
    typename Image_data_structure::Face_type>;

  using Random = CGAL::Random;
  using Indices = std::vector<std::size_t>;
  using Sphere_neighbor_query =
    internal::Sphere_neighbor_query<Traits, Indices, Point_map>;

  LOD2_image_reconstruction(
    const Indices& all_points,
    const Point_map& point_map,
    std::vector<Segment_2>& boundary,
    const std::vector<Segment_2>& directions,
    const Triangulation& lod0,
    ImagePointer& image_ptr,
    Partition_2& partition_2,
    const FT noise_level_2,
    const FT min_length_2,
    const FT angle_bound_2,
    const FT ordinate_bound_2,
    const FT max_height_difference,
    const FT beta,
    const FT top_z) :
  m_all_points(all_points),
  m_point_map(point_map),
  m_boundary(boundary),
  m_directions(directions),
  m_lod0(lod0),
  m_image_ptr(image_ptr),
  m_partition_2(partition_2),
  m_noise_level_2(noise_level_2),
  m_min_length_2(min_length_2),
  m_angle_bound_2(angle_bound_2),
  m_ordinate_bound_2(ordinate_bound_2),
  m_max_height_difference(max_height_difference),
  m_beta(beta),
  m_top_z(top_z),
  m_pi(static_cast<FT>(CGAL_PI)),
  m_random(0),
  m_image_creator(
    m_image_ptr, m_boundary, m_noise_level_2),
    m_snq(m_all_points, FT(1), m_point_map) {

    m_partition_2.clear();
  }

  void build() {
    m_image_creator.create_image();
    m_image_creator.clean_image();
    m_image_creator.create_label_pairs();

    /* m_image_creator.create_ridges_with_contours_v1(); */
    /* m_image_creator.create_ridges_with_contours_v2(); */

    m_image_creator.create_ridges_with_contours_v3();
    m_data_structure_ptr = std::make_shared<Image_data_structure>(
      m_boundary,
      m_image_creator.get_ridges(),
      m_image_creator.get_image(),
      m_image_ptr->get_plane_map(),
      m_noise_level_2,
      m_min_length_2,
      m_angle_bound_2,
      m_ordinate_bound_2,
      m_max_height_difference,
      m_top_z);

    m_data_structure_ptr->clear();
    m_data_structure_ptr->build();
  }

  void simplify() {

    if (m_data_structure_ptr->faces().size() == 0)
      return;

    m_data_structure_ptr->simplify();
    std::cout << "data structure simplified" << std::endl;
  }

  void create_tree() {

    if (m_data_structure_ptr->faces().size() == 0)
      return;

    /* m_data_structure_ptr->save_faces_ply("faces"); */
    m_tree_ptr = std::make_shared<Image_tree>(
      m_boundary,
      m_directions,
      m_image_ptr->get_plane_map(),
      m_min_length_2,
      m_max_height_difference,
      m_beta,
      m_data_structure_ptr->vertices(),
      m_data_structure_ptr->edges(),
      m_data_structure_ptr->halfedges(),
      m_data_structure_ptr->faces());

    m_tree_ptr->build();
    m_tree_ptr->check_tree_information(false, true);

    /* m_tree_ptr->check_vertex_information(); */
    /* m_tree_ptr->check_edge_information(); */
    /* m_tree_ptr->check_halfedge_information(); */
    /* m_tree_ptr->check_face_information(); */
    /* m_tree_ptr->apply_test(); */

    for (std::size_t i = 0; i < m_tree_ptr->num_levels(); ++i) {
      m_tree_ptr->cut(i);
      m_data_structure_ptr->save_all_faces_ply(i, "tree");
    }

    if (m_tree_ptr->num_levels() == 0)
      return;

    if (m_data_structure_ptr->faces().size() > 1) {

      std::cout << "cut 1" << std::endl;
      m_tree_ptr->cut(1); // 1 - base level
      std::cout << "cut finished" << std::endl;
      m_tree_ptr->merge_faces();
      /* m_data_structure_ptr->save_faces_ply("faces1"); */
    }

    if (m_data_structure_ptr->faces().size() > 1) {

      std::cout << "cut 2" << std::endl;
      m_tree_ptr->remove_one_neighbor_faces();
      std::cout << "cut finished" << std::endl;
      m_tree_ptr->merge_faces();
      /* m_data_structure_ptr->save_faces_ply("faces2"); */
    }

    std::cout << "data structure hierarchy built" << std::endl;
  }

  void regularize() {

    if (m_data_structure_ptr->faces().size() <= 1)
      return;

    std::cout << "regularize" << std::endl;
    m_data_structure_ptr->regularize(0);
    m_data_structure_ptr->project_linear(1);
    m_data_structure_ptr->regularize(2);
    m_data_structure_ptr->snap(3);
    m_data_structure_ptr->regularize(4);
    m_data_structure_ptr->merge_corners(5);
    m_data_structure_ptr->regularize(6);
    m_data_structure_ptr->merge_free_parts(7);
    m_data_structure_ptr->regularize(8);

    std::cout << "data structure regularized" << std::endl;
  }

  void get_roof_planes(
    std::vector<Plane_3>& roof_planes) {

    const auto& plane_map = m_image_ptr->get_plane_map();
    roof_planes.clear();
    roof_planes.reserve(plane_map.size());

    for (const auto& pair : plane_map) {
      const auto& plane = pair.second;
      roof_planes.push_back(plane);
    }
  }

  void get_lod2(Building& building) {

    const auto& edges0 = building.edges0;
    CGAL_assertion(!edges0.empty());
    auto& edges2 = building.edges2;
    edges2 = edges0;

    const auto& base0 = building.base0;
    CGAL_assertion(!base0.empty());
    auto& base2 = building.base2;
    base2 = base0;

    m_data_structure_ptr->create_all_face_edges();
    create_building_walls(
      building.bottom_z, building.walls2);
    create_building_roofs(
      building.top_z, building.roofs2);
  }

private:
  const Indices& m_all_points;
  const Point_map& m_point_map;
  std::vector<Segment_2>& m_boundary;
  const std::vector<Segment_2>& m_directions;
  const Triangulation& m_lod0;
  Image_ptr& m_image_ptr;
  Partition_2& m_partition_2;
  const FT m_noise_level_2;
  const FT m_min_length_2;
  const FT m_angle_bound_2;
  const FT m_ordinate_bound_2;
  const FT m_max_height_difference;
  const FT m_beta;
  const FT m_top_z;
  const FT m_pi;

  Random m_random;
  Image_creator m_image_creator;
  std::shared_ptr<Image_data_structure> m_data_structure_ptr;
  std::shared_ptr<Image_tree> m_tree_ptr;
  Sphere_neighbor_query m_snq;

  void create_building_walls(
    const FT bottom_z,
    std::vector<Wall>& walls) {

    std::vector<Segment_3> outer_segments_3;
    m_data_structure_ptr->get_wall_outer_segments(outer_segments_3);

    std::vector<Segment_3> inner_segments_3;
    m_data_structure_ptr->get_wall_inner_segments(inner_segments_3);

    walls.clear();
    add_walls(bottom_z, outer_segments_3, walls);
    for (const auto& wall : walls)
      for (const auto& triangle : wall.triangles)
        add_points(triangle);

    add_walls(bottom_z, inner_segments_3, walls);
    /*
    for (const auto& wall : walls)
      for (const auto& triangle : wall.triangles)
        add_points(triangle); */
  }

  void add_walls(
    const FT bottom_z,
    const std::vector<Segment_3>& segments_3,
    std::vector<Wall>& walls) {

    if (segments_3.empty()) return;
    for (const auto& segment_3 : segments_3) {
      const Point_3& s = segment_3.source();
      const Point_3& t = segment_3.target();

      const Point_3 p1 = Point_3(s.x(), s.y(), bottom_z);
      const Point_3 p2 = Point_3(t.x(), t.y(), bottom_z);
      const Point_3 p3 = Point_3(t.x(), t.y(), t.z());
      const Point_3 p4 = Point_3(s.x(), s.y(), s.z());

      Wall wall;
      wall.triangles.push_back(Triangle_3(p1, p2, p3));
      wall.triangles.push_back(Triangle_3(p3, p4, p1));
      walls.push_back(wall);
    }
  }

  void add_points(
    const Triangle_3& triangle) {

    return;
    std::ofstream outfile;
    outfile.precision(20);
    outfile.open(
      "/Users/monet/Documents/lod/logs/buildings/02_building_interior_points.ply",
      std::ios_base::app);

    std::vector<Point_3> samples;
    using Point_generator = CGAL::Random_points_in_triangle_3<Point_3>;

    Point_generator generator(triangle, m_random);
    std::copy_n(
      generator, 200, std::back_inserter(samples));
    for (const auto& sample : samples)
      outfile << sample << " 255 102 51 " << std::endl;
    outfile.close();
  }

  void create_building_roofs(
    const FT top_z,
    std::vector<Roof>& roofs) {

    std::vector<Triangle_set_3> triangle_sets_3;
    m_data_structure_ptr->get_roof_triangles(triangle_sets_3);

    roofs.clear();
    if (triangle_sets_3.empty()) return;
    for (const auto& triangle_set_3 : triangle_sets_3) {

      Roof roof;
      for (const auto& triangle_3 : triangle_set_3)
        roof.triangles.push_back(triangle_3);
      roofs.push_back(roof);
    }
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_LOD2_IMAGE_RECONSTRUCTION_H
