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

    using Trees = internal::Trees<Data_structure>;
    using Buildings = internal::Buildings<Data_structure>;
    using Ground = internal::Ground<Data_structure>;
    
    using Ground_base = typename Ground::Ground_base;
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
    output(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      Reconstruction_type lod_type) const {

      if (empty())
        return boost::none;

      switch (lod_type) {
        case Reconstruction_type::LOD0: {
          return lod0(vertices, faces); }
        case Reconstruction_type::LOD1: {
          return lod1(vertices, faces); }
        case Reconstruction_type::LOD2: {
          return lod2(vertices, faces); }
        default: {
          return boost::none; }
      }
    }

    bool empty() const {
      return ( m_ground.empty() && m_trees.empty() && m_buildings.empty() );
    }

  private:
    const Data_structure& m_data;
    const Ground& m_ground;
    const Trees& m_trees;
    const Buildings& m_buildings;

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod0(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      // Create ground.
      Ground_base ground_base;
      m_ground.initialize(ground_base);
      auto builder_ptr = m_ground.make_planar_builder(ground_base);
      add_footprints(*builder_ptr, Reconstruction_type::LOD0);

      // Output.
      return ground_base.output(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod1(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      return lod12(vertices, faces, Reconstruction_type::LOD1);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod2(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      return lod12(vertices, faces, Reconstruction_type::LOD2);
    }

    template<typename Builder>
    void add_footprints(
      Builder& builder, 
      const Reconstruction_type type) const {

      const auto& trees = m_trees.trees();
      const auto& buildings = m_buildings.buildings();

      builder.initialize();
      if (trees) add_urban_objects(builder, *trees, type);
      if (buildings) add_urban_objects(builder, *buildings, type);
      if (trees) tag_faces(builder, *trees);
      if (buildings) tag_faces(builder, *buildings);
      builder.finilize();
    }

    template<
    typename Builder,
    typename Urban_object>
    void add_urban_objects(
      Builder& builder, 
      const std::vector<Urban_object>& objects,
      const Reconstruction_type type) const {
      for (const auto& object : objects)
        builder.add_urban_object(object, type);
    }

    template<
    typename Builder,
    typename Urban_object>
    void tag_faces(
      Builder& builder, 
      const std::vector<Urban_object>& objects) const {
      for (const auto& object : objects)
        builder.tag_faces(object);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lod12(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type type) const {

      // Create ground.
      Ground_base ground_base;
      m_ground.initialize(ground_base);
      auto neighbor_query_ptr = m_ground.make_neighbor_query();
      auto builder_ptr = m_ground.make_smooth_builder(
        ground_base, *neighbor_query_ptr);
      add_footprints(*builder_ptr, type);

      // Output.
      Indexer indexer;
      std::size_t num_vertices = 0;
  
      const auto& trees = m_trees.trees();
      if (trees) output_urban_objects(*trees, type,
        indexer, num_vertices, vertices, faces);
      
      const auto& buildings = m_buildings.buildings();
      if (buildings) output_urban_objects(*buildings, type,
        indexer, num_vertices, vertices, faces);

      return ground_base.output(vertices, faces);
    }

    template<
    typename Urban_object,
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_urban_objects(
      const std::vector<Urban_object>& objects,
      const Reconstruction_type type,
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      switch (type) {
        case Reconstruction_type::LOD1: {
          for (const auto& object : objects)
            object.output_lod1(indexer, num_vertices, vertices, faces);
          return;
        }
        case Reconstruction_type::LOD2: {
          for (const auto& object : objects)
            object.output_lod2(indexer, num_vertices, vertices, faces);
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
