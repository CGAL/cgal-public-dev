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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_KINETIC_PARTITIONING_3_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_KINETIC_PARTITIONING_3_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <vector>
#include <algorithm>
#include <unordered_map>

// CGAL includes.
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

// Kinetic includes.
// #include "kinetic3/defs_cgal.h"
// #include "kinetic3/universe.h"
// #include "kinetic3/propagation_simple.h"
// #include "kinetic3/propagation_multiple.h"
// #include "kinetic3/support_plane.h"

#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic/include/defs_cgal.h>
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic/include/universe.h>
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic/include/propagation_simple.h>
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic/include/propagation_multiple.h>
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic/include/support_plane.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<class GeomTraits>
  class Kinetic_partitioning_3 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;

    using Partition_3 = internal::Partition_3<Traits>;
    using Edge = typename Partition_3::Edge;
    using Face = typename Partition_3::Face;

    using Polygon = std::vector<Point_3>;

    // Kinetic.
    using JP_FT = Skippy::FT;
    using JP_point_3 = Skippy::CGAL_Point_3;
    using JP_polygon = std::vector<JP_point_3>;
    using JP_polygons = std::vector<JP_polygon>;
    /* using JP_kinetic_propagation = Skippy::Kinetic_Propagation_Simple; */
    using JP_kinetic_propagation = Skippy::Kinetic_Propagation_Multiple;
    using JP_polyhedron = Skippy::Partition_Polyhedron;
    using JP_vertex = Skippy::Partition_Vertex;
    using JP_sequence  = std::list<JP_vertex*>;
    using JP_sequences = std::list<JP_sequence>;
    using JP_conversions = std::map<const JP_vertex*, int>;
    using JP_sequence_set = std::set<JP_vertex*>;
    using JP_edge = Skippy::Partition_Edge;
    using JP_facet_vertices = std::vector<JP_vertex*>;
    using JP_facet = Skippy::Partition_Facet;
    using Random = CGAL::Random;
    using Timer = CGAL::Timer;

    Kinetic_partitioning_3(
      std::vector<Edge>& outer_walls,
      std::vector<Edge>& inner_walls,
      std::vector<Edge>& roofs,
      Edge& ground,
      const std::size_t kinetic_max_intersections_3) :
    m_outer_walls(outer_walls),
    m_inner_walls(inner_walls),
    m_roofs(roofs),
    m_ground(ground),
    m_max_intersections(kinetic_max_intersections_3),
    m_up_scale(FT(3)),
    m_down_scale(FT(99) / FT(100)),
    m_z_scale(FT(10)),
    m_fixed_disc_radius(FT(1) / FT(1000)),
    m_num_points_in_disc(25),
    m_random(0)
    { }

    void compute(Partition_3& partition) {

      /*
      if (m_outer_walls.empty() || m_roofs.empty())
        return; */

      JP_polygons jp_polygons;
      const std::size_t input_size =
      1 + m_outer_walls.size() + m_inner_walls.size() + m_roofs.size();
      jp_polygons.reserve(input_size);
      set_ground(jp_polygons);
      set_outer_walls(jp_polygons);
      set_inner_walls(jp_polygons);
      set_roofs(jp_polygons);
      CGAL_assertion(jp_polygons.size() == input_size);

      create_partitioning(jp_polygons, partition);
    }

  private:
    std::vector<Edge>& m_outer_walls;
    std::vector<Edge>& m_inner_walls;
    std::vector<Edge>& m_roofs;
    Edge& m_ground;
    Random m_random;

    // External parameters.
    const std::size_t m_max_intersections;

    // Internal parameters.
    const FT m_up_scale;
    const FT m_down_scale;
    const FT m_z_scale;
    const FT m_fixed_disc_radius;
    const std::size_t m_num_points_in_disc;

    void set_ground(JP_polygons& jp_polygons) {
      process_polygon(m_ground.polygon, jp_polygons, m_up_scale, FT(1));
    }

    void set_outer_walls(JP_polygons& jp_polygons) {
      for (auto& wall : m_outer_walls)
        process_polygon(wall.polygon, jp_polygons, m_down_scale, m_z_scale);
    }

    void set_inner_walls(JP_polygons& jp_polygons) {
      for (auto& wall : m_inner_walls)
        process_polygon(wall.polygon, jp_polygons, m_up_scale, m_z_scale);
    }

    void set_roofs(JP_polygons& jp_polygons) {
      for (auto& roof : m_roofs)
        process_polygon(roof.polygon, jp_polygons, m_up_scale, FT(1));
    }

    void process_polygon(
      Polygon &polygon,
      JP_polygons& jp_polygons,
      const FT scale,
      const FT z_extender) {

      if (polygon.size() == 0) return;
      internal::scale_polygon_3(scale, z_extender, polygon);
      internal::perturb_polygon_vertices_3(
        m_fixed_disc_radius, m_num_points_in_disc, m_random, polygon);

      JP_polygon jp_polygon;
      jp_polygon.reserve(polygon.size());

      for (const auto& p : polygon)
        jp_polygon.push_back(
          JP_point_3(JP_FT(p.x()), JP_FT(p.y()), JP_FT(p.z())));
      CGAL_assertion(jp_polygon.size() == polygon.size());
      jp_polygons.push_back(jp_polygon);
    }

    void create_partitioning(
      const JP_polygons& jp_polygons,
      Partition_3& partition) const {

      JP_kinetic_propagation kinetic(jp_polygons);
      Skippy::Universe::params->K = m_max_intersections;

      /*
      Skippy::Universe::params->output_polyhedrons = true;
      Skippy::Universe::params->basename = "/Users/monet/Documents/lod/logs/polyhedrons/kinetic"; */

      // Subdivision options.
      Skippy::Universe::params->use_grid = true;
      Skippy::Universe::params->grid_x = 2;
      Skippy::Universe::params->grid_y = 2;
      Skippy::Universe::params->grid_z = 2;

	    if (!kinetic.data()) return;

      Timer timer;
      timer.start();
      kinetic.run();
      timer.stop();
      std::cout << "Finished in " << timer.time() << " seconds!" << std::endl;

      set_output(kinetic, partition);

      /* kinetic.delete_unique_kinetic_data_structure(); */ // for simple propagation
    }

    void set_output(
      const JP_kinetic_propagation& kinetic,
      Partition_3& partition) const {

      std::unordered_map<int, int> fmap;
      add_faces(kinetic, partition, fmap);
      add_edges(kinetic, fmap, partition);
    }

    void add_faces(
      const JP_kinetic_propagation& kinetic,
      Partition_3& partition,
      std::unordered_map<int, int>& fmap) const {

      partition.faces.clear(); int face_id = 0;
      for (auto it = kinetic.partition->polyhedrons_begin();
      it != kinetic.partition->polyhedrons_end(); ++it) {

        add_face(kinetic.partition->planes, *it, partition.faces);
        CGAL_assertion((*it)->id >= 0);
        fmap[(*it)->id] = face_id;
        ++face_id;
      }

      // Neighbors.
      std::size_t i = 0;
      for (auto it = kinetic.partition->polyhedrons_begin();
      it != kinetic.partition->polyhedrons_end(); ++it, ++i) {
        const JP_polyhedron* polyhedron = *it;

        const int poly_id = polyhedron->id;
        partition.faces[i].neighbors.clear();

        for (auto fit = polyhedron->facets_begin();
        fit != polyhedron->facets_end(); ++fit) {
          const JP_facet* facet = fit->first;

          int id = -1;
          const JP_polyhedron* poly1 = facet->get_polyhedron_1();
          const JP_polyhedron* poly2 = facet->get_polyhedron_2();

          // Neighbors.
          if (facet->p < 6) {
            CGAL_assertion(poly1 == nullptr || poly2 == nullptr);
            id = -1;
          } else {
            CGAL_assertion(poly1 != nullptr && poly2 != nullptr);
            if (poly1->id == poly_id) {
              CGAL_assertion(fmap.find(poly2->id) != fmap.end());
              id = fmap.at(poly2->id);
            } else {
              CGAL_assertion(fmap.find(poly1->id) != fmap.end());
              id = fmap.at(poly1->id);
            }
          }
          partition.faces[i].neighbors.push_back(id);
        }
      }
    }

    void add_face(
      const std::vector<Skippy::CGAL_Plane>& planes,
      const JP_polyhedron* polyhedron,
      std::vector<Face>& faces) const {

      JP_conversions conversions;
      JP_sequences sequences_per_side;

      Face face;
      get_polyhedron_vertices(
        planes,
        polyhedron,
        sequences_per_side, conversions,
        face.vertices);
      get_polyhedron_faces(
        sequences_per_side, conversions,
        face.faces);
      faces.push_back(face);
    }

    void get_polyhedron_vertices(
      const std::vector<Skippy::CGAL_Plane>& planes,
      const JP_polyhedron* polyhedron,
      JP_sequences& sequences_per_side,
      JP_conversions& conversions,
      std::vector<Point_3>& vertices) const {

      /*
      std::vector<Skippy::CGAL_Plane> planes;
      planes.reserve(Skippy::Universe::map_of_planes.size());
      for (auto it = Skippy::Universe::map_of_planes.begin();
      it != Skippy::Universe::map_of_planes.end(); ++it)
        planes.push_back((*it)->plane); */

      // Sequences.
      sequences_per_side.clear();
      JP_sequence_set vertices_used;
      JP_sequence facet_vertices;

      for (auto fit = polyhedron->facets_begin();
      fit != polyhedron->facets_end(); ++fit) {
        const JP_facet* facet = fit->first;

        facet_vertices.clear();
        facet->get_circular_sequence_of_vertices(planes, facet_vertices, !fit->second);
        for (auto vit = facet_vertices.begin();
        vit != facet_vertices.end(); ++vit)
          vertices_used.insert(*vit);

        sequences_per_side.push_back(facet_vertices);
      }

      // Vertices.
      vertices.clear();
      conversions.clear();
      for (auto vit = vertices_used.begin();
      vit != vertices_used.end(); ++vit) {

        const auto* v = *vit;
        const auto& p = v->M;

        const FT x = static_cast<FT>(CGAL::to_double(p.x()));
        const FT y = static_cast<FT>(CGAL::to_double(p.y()));
        const FT z = static_cast<FT>(CGAL::to_double(p.z()));

        const Point_3 vertex = Point_3(x, y, z);
        vertices.push_back(vertex);
        conversions[v] = static_cast<int>(vertices.size()) - 1;
      }
    }

    void get_polyhedron_faces(
      const JP_sequences& sequences_per_side,
      const JP_conversions& conversions,
      std::vector< std::vector<std::size_t> >& faces) const {

      // Faces.
      faces.clear(); std::vector<std::size_t> face;
      for (auto sit = sequences_per_side.begin();
      sit != sequences_per_side.end(); ++sit) {
        const auto& sequence = *sit;

        face.clear();
        for (auto vit = sequence.begin(); vit != sequence.end(); ++vit) {
          const int idx = conversions.at(*vit);
          CGAL_assertion(idx >= 0);
          face.push_back(static_cast<std::size_t>(idx));
        }
        faces.push_back(face);
      }
    }

    void add_edges(
      const JP_kinetic_propagation& kinetic,
      const std::unordered_map<int, int>& fmap,
      Partition_3& partition) const {

      /*
      std::vector<Skippy::CGAL_Plane> planes;
      planes.reserve(Skippy::Universe::map_of_planes.size());
      for (auto it = Skippy::Universe::map_of_planes.begin();
      it != Skippy::Universe::map_of_planes.end(); ++it)
        planes.push_back((*it)->plane); */

      const auto& planes = kinetic.partition->planes;
      auto& edges = partition.edges;
      edges.clear();

      JP_facet_vertices v;
      kinetic.partition->get_all_vertices_sorted_by_identifier(v);
      std::vector<Point_3> vertices;

      CGAL_assertion(v.size() > 0);
      for (const auto& vi : v) {
        const auto& p = vi->M;

        const FT x = static_cast<FT>(CGAL::to_double(p.x()));
        const FT y = static_cast<FT>(CGAL::to_double(p.y()));
        const FT z = static_cast<FT>(CGAL::to_double(p.z()));
        vertices.push_back(Point_3(x, y, z));
      }

      JP_sequence facet_vertices;
      const auto& facets = kinetic.partition->facets;
      std::vector<int> indices;

      Edge edge; int id1, id2;
      for (const auto& facet : facets) {
        for (auto fit = facet.begin(); fit != facet.end(); ++fit) {
          const JP_facet* f = *fit;

          const JP_polyhedron* poly1 = f->get_polyhedron_1();
          const JP_polyhedron* poly2 = f->get_polyhedron_2();

          // Neighbors.
          if (f->p < 6) {
            CGAL_assertion(poly1 == nullptr || poly2 == nullptr);
            if (poly1 != nullptr) {
              CGAL_assertion(fmap.find(poly1->id) != fmap.end());
              id1 = fmap.at(poly1->id);
              id2 = -1;
            } else {
              CGAL_assertion(fmap.find(poly2->id) != fmap.end());
              id1 = -1;
              id2 = fmap.at(poly2->id);
            }
          } else {
            CGAL_assertion(poly1 != nullptr && poly2 != nullptr);
            CGAL_assertion(fmap.find(poly1->id) != fmap.end());
            CGAL_assertion(fmap.find(poly2->id) != fmap.end());
            id1 = fmap.at(poly1->id);
            id2 = fmap.at(poly2->id);
          }

          edge.neighbors.first = id1;
          edge.neighbors.second = id2;

          // Edges.
          facet_vertices.clear();
          f->get_circular_sequence_of_vertices(planes, facet_vertices, true);

          indices.clear();
          for (auto vit = facet_vertices.begin();
          vit != facet_vertices.end(); ++vit)
            indices.push_back((*vit)->id);

          edge.polygon.clear();
          edge.polygon.reserve(indices.size());
          for (const int idx : indices) {
            CGAL_assertion(idx >= 0);
            edge.polygon.push_back(vertices[idx]);
          }

          // 6 faces of the bbox + 1 ground face + num walls
          const std::size_t rem = 7 + m_outer_walls.size() + m_inner_walls.size();
          if (f->p >= rem)
            edge.plane_index = f->p - rem;
          edges.push_back(edge);
        }
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_KINETIC_PARTITIONING_3_H
