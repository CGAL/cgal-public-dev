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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_RECONSTRUCTION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_RECONSTRUCTION_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <memory>
#include <vector>
#include <utility>

// Boost includes.
#include <boost/optional/optional.hpp>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/Trees/Trees.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Buildings.h>
#include <CGAL/Levels_of_detail/internal/Ground/Ground.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Reconstruction {

  public:
    using Data_structure = DataStructure;
    using Traits = typename Data_structure::Traits;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    using Ground = internal::Ground<Data_structure>;
    using Trees = internal::Trees<Data_structure>;
    using Buildings = internal::Buildings<Data_structure>;

    using Ground_base = typename Ground::Ground_base;
    using Triangulation = typename Ground_base::Triangulation;
    using Tree = internal::Tree<Traits>;
    using Tree_ptr = std::shared_ptr<Tree>;
    using Building = internal::Building<Traits>;
    using Building_ptr = std::shared_ptr<Building>;

    using Indexer = typename internal::Indexer<typename Traits::Point_3>;

    Reconstruction(
      const Data_structure& data,
      const Ground& ground,
      const Trees& trees,
      const Buildings& buildings) :
    m_data(data),
    m_ground(ground),
    m_trees(trees),
    m_buildings(buildings)
    { }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_lod(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type lod_type,
      const FT ground_precision) const {

      if (empty())
        return boost::none;

      switch (lod_type) {
        case Reconstruction_type::LOD0: {
          return lod0(vertices, faces); }
        case Reconstruction_type::LOD1: {
          return lod1(vertices, faces, ground_precision); }
        case Reconstruction_type::LOD2: {
          return lod2(vertices, faces, ground_precision); }
        default: {
          return boost::none; }
      }
    }

    bool empty() const {
      return ( m_ground.empty() && m_trees.empty() && m_buildings.empty() );
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire0(
      OutputIterator output) const {

      std::vector<Tree_ptr> trees;
      m_trees.get_trees(trees);
      std::vector<Building_ptr> buildings;
      m_buildings.get_buildings(buildings);

      // Create ground.
      Ground_base ground_base;
      m_ground.initialize(ground_base);
      auto builder_ptr = m_ground.make_planar_builder(ground_base);
      add_footprints(*builder_ptr, trees, buildings, Reconstruction_type::LOD0);

      // Output ground wire.
      ground_base.output_boundary_edges(output);
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire1(
      OutputIterator output,
      const FT ground_precision) const {
      return output_wire12(Reconstruction_type::LOD1, ground_precision, output);
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire2(
      OutputIterator output,
      const FT ground_precision) const {
      return output_wire12(Reconstruction_type::LOD2, ground_precision, output);
    }

  private:
    const Data_structure& m_data;
    const Ground& m_ground;
    const Trees& m_trees;
    const Buildings& m_buildings;

    template<typename OutputIterator>
    boost::optional<OutputIterator>
    output_wire12(
      const Reconstruction_type type,
      const FT ground_precision,
      OutputIterator output) const {

      std::vector<Tree_ptr> trees;
      m_trees.get_trees(trees);
      std::vector<Building_ptr> buildings;
      m_buildings.get_buildings(buildings);

      // Create ground.
      Ground_base ground_base;
      m_ground.initialize(ground_base);
      auto neighbor_query_ptr = m_ground.make_neighbor_query();
      auto builder_ptr = m_ground.make_smooth_builder(
        ground_base, *neighbor_query_ptr, ground_precision);
      add_footprints(*builder_ptr, trees, buildings, type);

      // Output objects wire.
      if (!trees.empty()) output_objects_wire(
        ground_base.triangulation,
        trees, type, output);
      if (!buildings.empty()) output_objects_wire(
        ground_base.triangulation,
        buildings, type, output);

      // Output ground wire.
      ground_base.output_all_edges(output);
      return output;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod0(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (m_data.verbose)
        std::cout << std::endl << "- Computing LOD0" << std::endl;

      std::vector<Tree_ptr> trees;
      m_trees.get_trees(trees);
      std::vector<Building_ptr> buildings;
      m_buildings.get_buildings(buildings);

      // Create ground.
      Ground_base ground_base;
      m_ground.initialize(ground_base);
      auto builder_ptr = m_ground.make_planar_builder(ground_base);
      add_footprints(*builder_ptr, trees, buildings, Reconstruction_type::LOD0);

      // Output.
      return ground_base.output_for_lod0(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod1(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const FT ground_precision) const {

      if (m_data.verbose)
        std::cout << std::endl << "- Computing LOD1" << std::endl;
      return lod12(vertices, faces, Reconstruction_type::LOD1, ground_precision);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod2(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const FT ground_precision) const {

      if (m_data.verbose)
        std::cout << std::endl << "- Computing LOD2" << std::endl;
      return lod12(vertices, faces, Reconstruction_type::LOD2, ground_precision);
    }

    template<typename Builder>
    void add_footprints(
      Builder& builder,
      const std::vector<Tree_ptr>& trees,
      const std::vector<Building_ptr>& buildings,
      const Reconstruction_type type) const {

      builder.initialize();
      if (!trees.empty()) add_constraints(builder, trees, type);
      if (!buildings.empty()) add_constraints(builder, buildings, type);
      if (!buildings.empty()) tag_faces(builder, buildings, type);
      if (!trees.empty()) tag_faces(builder, trees, type);
      builder.finilize();
    }

    template<
    typename Builder,
    typename Urban_object>
    void add_constraints(
      Builder& builder,
      const std::vector<Urban_object>& objects,
      const Reconstruction_type type) const {

      switch (type) {
        case Reconstruction_type::LOD0: {
          for (const auto& object : objects)
            builder.add_constraints(object->edges0);
          return;
        }
        case Reconstruction_type::LOD1: {
          for (const auto& object : objects)
            builder.add_constraints(object->edges1);
          return;
        }
        case Reconstruction_type::LOD2: {
          for (const auto& object : objects)
            builder.add_constraints(object->edges2);
          return;
        }
        default: {
          return;
        }
      }
    }

    template<
    typename Builder,
    typename Urban_object>
    void tag_faces(
      Builder& builder,
      const std::vector<Urban_object>& objects,
      const Reconstruction_type type) const {

      switch (type) {
        case Reconstruction_type::LOD0: {
          for (const auto& object : objects)
            builder.tag_faces(*object, object->base0.triangulation.delaunay);
          return;
        }
        case Reconstruction_type::LOD1: {
          for (const auto& object : objects)
            builder.tag_faces(*object, object->base1.triangulation.delaunay);
          return;
        }
        case Reconstruction_type::LOD2: {
          for (const auto& object : objects)
            builder.tag_faces(*object, object->base2.triangulation.delaunay);
          return;
        }
        default: {
          return;
        }
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod12(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type type,
      const FT ground_precision) const {

      std::vector<Tree_ptr> trees;
      m_trees.get_trees(trees);
      std::vector<Building_ptr> buildings;
      m_buildings.get_buildings(buildings);

      // Create ground.
      Ground_base ground_base;
      m_ground.initialize(ground_base);
      auto neighbor_query_ptr = m_ground.make_neighbor_query();
      auto builder_ptr = m_ground.make_smooth_builder(
        ground_base, *neighbor_query_ptr, ground_precision);
      add_footprints(*builder_ptr, trees, buildings, type);

      // Output.
      Indexer indexer;
      std::size_t num_vertices = 0;

      if (!trees.empty()) output_urban_objects(
        ground_base.triangulation,
        trees, "trees", type,
        indexer, num_vertices, vertices, faces);
      if (!buildings.empty()) output_urban_objects(
        ground_base.triangulation,
        buildings, "buildings", type,
        indexer, num_vertices, vertices, faces);
      return ground_base.output_for_lod12(indexer, num_vertices, vertices, faces);
    }

    template<
    typename Urban_object,
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_urban_objects(
      const Triangulation& tri,
      const std::vector<Urban_object>& objects,
      const std::string name,
      const Reconstruction_type type,
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      if (m_data.verbose)
        std::cout << "* adding " << name << std::endl;

      switch (type) {
        case Reconstruction_type::LOD1: {
          for (const auto& object : objects)
            object->output_lod1(tri, indexer, num_vertices, vertices, faces);
          return;
        }
        case Reconstruction_type::LOD2: {
          for (const auto& object : objects)
            object->output_lod2(tri, indexer, num_vertices, vertices, faces);
          return;
        }
        default: {
          return;
        }
      }
    }

    template<
    typename Urban_object,
    typename OutputIterator>
    void output_objects_wire(
      const Triangulation& tri,
      const std::vector<Urban_object>& objects,
      const Reconstruction_type type,
      OutputIterator output) const {

      switch (type) {
        case Reconstruction_type::LOD1: {
          for (const auto& object : objects)
            object->output_lod1_wire(tri, output, true);
          return;
        }
        case Reconstruction_type::LOD2: {
          for (const auto& object : objects)
            object->output_lod2_wire(tri, output, true);
          return;
        }
        default: {
          return;
        }
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_RECONSTRUCTION_H
