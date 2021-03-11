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

#ifndef CGAL_LEVELS_OF_DETAIL_H
#define CGAL_LEVELS_OF_DETAIL_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <utility>

// Boost includes.
#include <boost/optional/optional.hpp>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/property_map.h>

// Internal components.
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/parameters.h>
#include <CGAL/Levels_of_detail/internal/Trees/Trees.h>
#include <CGAL/Levels_of_detail/internal/Ground/Ground.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Buildings.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/Reconstruction.h>

namespace CGAL {
namespace Levels_of_detail {

  /*!
    \ingroup PkgLevelsOfDetailRef

    \brief Given as input a point set or polygon soup, which represents urban
    environment, this class turns this input into a 3D model with Levels Of Detail (LODs).

    The class handles

    - items labeled as ground to turn them into either planar or smooth
    represenation of the ground;
    - items labeled as vegetation to turn them into trees represented either
    as discs or cylinders or trunks + crowns;
    - items labeled as buildings to turn them into buildings represented either
    as polygons or extruded polygons or walls + roofs.

    All above objects can be combined into LOD0, LOD1, or LOD2 complete 3D model
    of the input.

    \tparam GeomTraits
    must be a model of `Kernel`.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam ItemMap
    must be an `LvaluePropertyMap` whose key type is the value type of the input
    range and value type is `Kernel::Point_3` or `std::vector<Kernel::Point_3>`.
    The latter represents a convex polygon in 3D.

    \tparam SemanticMap
    must be an `LvaluePropertyMap` whose key type is the value type of the
    input range and value type is `CGAL::Levels_of_detail::Semantic_label`.

    \tparam VisibilityMap
    must be an `LvaluePropertyMap` whose key type is the value type of the
    input range and value type is `double`.
    %Default is `CGAL::Levels_of_detail::Visibility_from_semantic_map`.

    \tparam Verbose
    must be either `CGAL::Tag_true` or `CGAL::Tag_false`.
    %Default is `CGAL::Tag_false`.
  */
  template<
    typename GeomTraits,
    typename InputRange,
    typename ItemMap,
    typename SemanticMap,
    typename VisibilityMap = Visibility_from_semantic_map<SemanticMap>,
    typename Verbose = CGAL::Tag_false>
	class Levels_of_detail {

	public:

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Item_map = ItemMap;
    using Semantic_map = SemanticMap;
    using Visibility_map = VisibilityMap;
    /// \endcond

    /// \name Types
    /// @{

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// @}

    /// \cond SKIP_IN_MANUAL
    using Data_structure = internal::Data_structure<
    Traits,
    Input_range,
    Item_map,
    Semantic_map,
    Visibility_map>;

    using Parameters = internal::Parameters<FT>;
    using Ground = internal::Ground<Data_structure>;
    using Buildings = internal::Buildings<Data_structure>;
    using Trees = internal::Trees<Data_structure>;
    using LODs = internal::Reconstruction<Data_structure>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param input_range
      an instance of `InputRange` with 3D points or convex polygons

      \param item_map
      an instance of `ItemMap` that maps an item from `input_range`
      to `Kernel::Point_3` or `std::vector<Kernel::Point_3>`

      \param semantic_map
      an instance of `SemanticMap` that maps an item from `input_range`
      to `CGAL::Levels_of_detail::Semantic_label`

      \param parameters
      a struct with LOD parameters.

      \param visibility_map
      an instance of `VisibilityMap` that maps an item from `input_range`
      to a visibility value in the range [0,1]

      \param input_type
      any of `CGAL::Levels_of_detail::Input_type`

      \pre `input_range.size() > 0`
    */
    Levels_of_detail(
      const InputRange& input_range,
      const ItemMap item_map,
      const SemanticMap semantic_map,
      const Parameters& parameters,
      const VisibilityMap visibility_map = VisibilityMap(),
      const Input_type input_type = Input_type::POINT_SET) :
    m_input_range(input_range),
    m_item_map(item_map),
    m_semantic_map(semantic_map),
    m_parameters(parameters),
    m_visibility_map(visibility_map),
    m_data(
      m_input_range,
      m_item_map,
      m_semantic_map,
      m_parameters,
      m_visibility_map,
      Verbose::value ? true : false),
    m_ground(m_data),
    m_buildings(m_data),
    m_trees(m_data),
    m_lods(m_data, m_ground, m_trees, m_buildings) {
      CGAL_precondition(m_input_range.size() > 0);
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    ~Levels_of_detail() {
      if (Verbose::value)
        std::cout << std::endl;
    }

    Data_structure& data() {
      return m_data;
    }

    const Data_structure& data() const {
      return m_data;
    }
    /// \endcond

    /// \name Complete Generation
    /// @{

    /*!
      \brief builds all available objects.

      This method:

      - builds all trees;

      - builds all buildings.
    */
    void build() {
      m_trees.make_trees();
      m_buildings.make_buildings();
    }

    /*!
      \brief builds all trees.
    */
    void build_trees() {
      m_trees.make_trees();
    }

    /*!
      \brief builds all buildings.
    */
    void build_buildings() {
      m_buildings.make_buildings();
    }

    /// @}

    /// \name Trees
    /// @{

    /*!
      \brief initializes trees.

      This method:

      - finds all connected components represented as clusters of points
      labeled as `CGAL::Levels_of_detail::Semantic_label::VEGETATION`;

      - initializes all internal data structures.
    */
    void initialize_trees() {
      m_trees.initialize();
    }

    /*!
      \brief computes tree footprints.

      The trees, after this step, are planar discs.

      This method:

      - extracts points, which form potential trees;

      - estimates each tree center, radius, and height;

      - creates tree footprints.

      \warning `initialize_trees()` should be called before calling this method
    */
    void compute_tree_footprints() {
      m_trees.compute_footprints();
    }

    /*!
      \brief extrudes tree footprints.

      The trees, after this step, are cylinders with a planar top.

      \warning `compute_tree_footprints()` should be called before calling this method
    */
    void extrude_tree_footprints() {
      m_trees.extrude_footprints();
    }

    /*!
      \brief computes tree crowns.

      The trees, after this step, consist of trunks and crowns.

      This method:

      - creates crown icons;

      - fits these icons to trees.

      \warning `extrude_tree_footprints()` should be called before calling this method
    */
    void compute_tree_crowns() {
      m_trees.compute_crowns();
    }

    /// @}

    /// \name Buildings
    /// @{

    /*!
      \brief initializes buildings.

      This method:

      - finds all connected components represented as clusters of points
      labeled as `CGAL::Levels_of_detail::Semantic_label::BUILDING_INTERIOR` and
      `CGAL::Levels_of_detail::Semantic_label::BUILDING_BOUNDARY`;

      - initializes all internal data structures.
    */
    void initialize_buildings() {
      m_buildings.initialize();
    }

    /*!
      \brief detects building boundaries.

      This method:

      - computes the alpha shape of the points labeled as
        `CGAL::Levels_of_detail::Semantic_label::BUILDING_INTERIOR` (if any) and
        `CGAL::Levels_of_detail::Semantic_label::BUILDING_BOUNDARY` (if any)
        and extracts the boundary points of this alpha shape;

      - downsamples this union of points using a regular grid;

      - contracts all points towards the union's skeleton;

      - detects line segments using the region growing approach.

      \warning `initialize_buildings()` should be called before calling this method
    */
    void detect_building_boundaries() {
      m_buildings.detect_boundaries();
    }

    /*!
      \brief computes building footprints.

      The buildings, after this step, are planar polygons.

      This method:

      - creates the partitioning by extending initial building boundary
        segments until the defined number of intersections with other segments
        is reached;

      - applies the visibility computation that assignes to each polygon face
        of the partitioning a visibility value in the range [0,1], where 0
        means certainly outside and 1 means certainly inside;

      - corrects the visibility estimations by applying a graphcut;

      - tags subsets of all polygon faces with the visibility value >= 0.5
        that form separate buildings.

      \warning `detect_building_boundaries()` should be called before calling this method
    */
    void compute_building_footprints() {
      m_buildings.compute_footprints();
    }

    /*!
      \brief extrudes building footprints.

      The buildings, after this step, are extruded polygons with a planar top.

      \warning `compute_building_footprints()` should be called before calling this method
    */
    void extrude_building_footprints() {
      m_buildings.extrude_footprints();
    }

    /*!
      \brief detects building roofs.

      This method:

      - detects chunks of 3D points that form planes using the region growing
        approach on all points labeled as `CGAL::Levels_of_detail::Semantic_label::BUILDING_INTERIOR`;

      - filters out all chunks, which do not fit to such criteria as
        verticality, size, etc;

      - creates convex polygons, which approximate all left chunks.

      \warning `extrude_building_footprints()` should be called before calling this method
    */
    void detect_building_roofs() {
      m_buildings.detect_roofs();
    }

    /*!
      \brief computes building roofs.

      The buildings, after this step, consist of walls and roofs.

      This method:

      - creates the partitioning by extending initial building walls, roofs, and
        base represented as polygons until the defined number of
        intersections with other polygons is reached;

      - applies the visibility computation that assignes to each polyhedral
        facet of the partitioning a visibility value in the range [0,1], where 0
        means certainly outside and 1 means certainly inside;

      - corrects the visibility estimations by applying a graphcut;

      - extracts polygons, which represent building walls and roofs.

      \warning `detect_building_roofs()` should be called before calling this method
    */
    void compute_building_roofs() {
      m_buildings.compute_roofs();
    }

    /// @}

    /// \name Access Final Objects
    /// @{

    /*!
      \brief returns ground as triangle soup.

      This method computes either a planar or a smooth representation
      of the ground.

      For a planar ground:
      The plane is estimated through the principal component analysis on
      the points semantically labeled as `CGAL::Levels_of_detail::Semantic_label::GROUND`.

      For a smooth ground:
      The ground is represented as a Delaunay triangulation with
      associated ground heights, which is computed upon the points
      semantically labeled as `CGAL::Levels_of_detail::Semantic_label::GROUND`.

      \tparam VerticesOutputIterator
      must be an output iterator whose value type is `Kernel::Point_3`.

      \tparam FacesOutputIterator
      must be an output iterator whose value type is
      `std::pair<std::vector<std::size_t>, CGAL::Levels_of_detail::Urban_object_type::GROUND>`.

      \param vertices
      an output iterator with vertices

      \param faces
      an output iterator with faces

      \param lod_type
      either `Reconstruction_type::PLANAR_GROUND` or `Reconstruction_type::SMOOTH_GROUND`

      \param ground_precision
      max distance between input points and a reconstructed ground

      \pre  `lod_type == Reconstruction_type::PLANAR_GROUND` ||
            `lod_type == Reconstruction_type::SMOOTH_GROUND`

      \return a pair of the past-the-end iterators
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    ground(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type lod_type,
      const FT ground_precision) const {

      CGAL_precondition(
        lod_type == Reconstruction_type::PLANAR_GROUND ||
        lod_type == Reconstruction_type::SMOOTH_GROUND);

      return m_ground.output_ground(vertices, faces, lod_type, ground_precision);
    }

    /*!
      \brief returns trees as triangle soup.

      This method returns different levels of detail for trees.

      \warning the corresponding tree-related functions should be called.

      \tparam VerticesOutputIterator
      must be an output iterator whose value type is `Kernel::Point_3`.

      \tparam FacesOutputIterator
      must be an output iterator whose value type is
      `std::pair<std::vector<std::size_t>, CGAL::Levels_of_detail::Urban_object_type::TREE_TRUNK>` or
      `std::pair<std::vector<std::size_t>, CGAL::Levels_of_detail::Urban_object_type::TREE_CROWN>`

      \param vertices
      an output iterator with vertices

      \param faces
      an output iterator with faces

      \param lod_type
      either `Reconstruction_type::TREES0` or `Reconstruction_type::TREES1` or
      `Reconstruction_type::TREES2`

      \pre  `lod_type == Reconstruction_type::TREES0` ||
            `lod_type == Reconstruction_type::TREES1` ||
            `lod_type == Reconstruction_type::TREES2`

      \return a pair of the past-the-end iterators
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    trees(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type lod_type) const {

      CGAL_precondition(
        lod_type == Reconstruction_type::TREES0 ||
        lod_type == Reconstruction_type::TREES1 ||
        lod_type == Reconstruction_type::TREES2);

      return m_trees.output_trees(vertices, faces, lod_type);
    }

    /*!
      \brief returns buildings as triangle soup.

      This method returns different levels of detail for buildings.

      \warning the corresponding building-related functions should be called.

      \tparam VerticesOutputIterator
      must be an output iterator whose value type is `Kernel::Point_3`.

      \tparam FacesOutputIterator
      must be an output iterator whose value type is
      `std::pair<std::vector<std::size_t>, CGAL::Levels_of_detail::Urban_object_type::BUILDING_WALL>` or
      `std::pair<std::vector<std::size_t>, CGAL::Levels_of_detail::Urban_object_type::BUILDING_ROOF>`

      \param vertices
      an output iterator with vertices

      \param faces
      an output iterator with faces

      \param lod_type
      either `Reconstruction_type::BUILDINGS0` or `Reconstruction_type::BUILDINGS1` or
      `Reconstruction_type::BUILDINGS2`

      \pre  `lod_type == Reconstruction_type::BUILDINGS0` ||
            `lod_type == Reconstruction_type::BUILDINGS1` ||
            `lod_type == Reconstruction_type::BUILDINGS2`

      \return a pair of the past-the-end iterators
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    buildings(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type lod_type) const {

      CGAL_precondition(
        lod_type == Reconstruction_type::BUILDINGS0 ||
        lod_type == Reconstruction_type::BUILDINGS1 ||
        lod_type == Reconstruction_type::BUILDINGS2);

      return m_buildings.output_buildings(vertices, faces, lod_type);
    }

    /*!
      \brief returns full levels of detail as triangle soup.

      This method constructs LODs by merging all available urban objects with
      the corresponding ground.

      \tparam VerticesOutputIterator
      must be an output iterator whose value type is `Kernel::Point_3`.

      \tparam FacesOutputIterator
      must be an output iterator whose value type is
      `std::pair< std::vector<std::size_t>, CGAL::Levels_of_detail::Urban_object_type>`.

      \param vertices
      an output iterator with vertices

      \param faces
      an output iterator with faces

      \param lod_type
      either `Reconstruction_type::LOD0` or `Reconstruction_type::LOD1` or
      `Reconstruction_type::LOD2`

      \param ground_precision
      max distance between input points and a reconstructed ground

      \pre  `lod_type == Reconstruction_type::LOD0` ||
            `lod_type == Reconstruction_type::LOD1` ||
            `lod_type == Reconstruction_type::LOD2`

      \return a pair of the past-the-end iterators
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    lods(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type lod_type,
      const FT ground_precision) const {

      CGAL_precondition(
        lod_type == Reconstruction_type::LOD0 ||
        lod_type == Reconstruction_type::LOD1 ||
        lod_type == Reconstruction_type::LOD2);

      return m_lods.output_lod(vertices, faces, lod_type, ground_precision);
    }

    /// @}

    /// \name Access Intermediate Steps
    /// @{

    /*!
      \brief returns a point set.

      This method returns data related to the intermediate steps of the
      algorithm, which can be stored as a point set.

      \tparam OutputIterator
      must be an output iterator whose value type is
      `std::pair<Kernel::Point_3, std::size_t>` or
      `std::pair<Kernel::Point_3, CGAL::Levels_of_detail::Semantic_label>`
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> points(
      OutputIterator output,
      const Intermediate_step step) const {

      CGAL_precondition(
        step == Intermediate_step::INPUT_GROUND_POINTS ||
        step == Intermediate_step::INPUT_VEGETATION_POINTS ||
        step == Intermediate_step::TREE_CLUSTERS ||
        step == Intermediate_step::TREE_POINTS ||
        step == Intermediate_step::INPUT_BUILDING_BOUNDARY_POINTS ||
        step == Intermediate_step::INPUT_BUILDING_INTERIOR_POINTS ||
        step == Intermediate_step::BUILDING_CLUSTERS ||
        step == Intermediate_step::BUILDING_BOUNDARY_POINTS ||
        step == Intermediate_step::BUILDING_WALL_POINTS ||
        step == Intermediate_step::BUILDING_POINTS ||
        step == Intermediate_step::BUILDING_ROOF_POINTS);

      switch (step) {
        case Intermediate_step::INPUT_GROUND_POINTS: {
          return m_data.get_points(Semantic_label::GROUND, output);
        }
        case Intermediate_step::INPUT_VEGETATION_POINTS: {
          return m_data.get_points(Semantic_label::VEGETATION, output);
        }
        case Intermediate_step::TREE_CLUSTERS: {
          return m_trees.get_tree_clusters(output);
        }
        case Intermediate_step::TREE_POINTS: {
          return m_trees.get_tree_points(output);
        }
        case Intermediate_step::INPUT_BUILDING_BOUNDARY_POINTS: {
          return m_data.get_points(Semantic_label::BUILDING_BOUNDARY, output);
        }
        case Intermediate_step::INPUT_BUILDING_INTERIOR_POINTS: {
          return m_data.get_points(Semantic_label::BUILDING_INTERIOR, output);
        }
        case Intermediate_step::BUILDING_CLUSTERS: {
          return m_buildings.get_building_clusters(output);
        }
        case Intermediate_step::BUILDING_BOUNDARY_POINTS: {
          return m_buildings.get_boundary_points(output);
        }
        case Intermediate_step::BUILDING_WALL_POINTS: {
          return m_buildings.get_wall_points(output);
        }
        case Intermediate_step::BUILDING_POINTS: {
          return m_buildings.get_building_points(output);
        }
        case Intermediate_step::BUILDING_ROOF_POINTS: {
          return m_buildings.get_roof_points(output);
        }
        default: return boost::none;
      }
    }

    /*!
      \brief returns polylines.

      This method returns data related to the intermediate steps of the
      algorithm, which can be stored as a set of polylines.

      \tparam OutputIterator
      must be an output iterator whose value type is
      `Kernel::Segment_3` or
      `std::pair<Kernel::Segment_3, std::size_t>`
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> polylines(
      OutputIterator output,
      const Intermediate_step step) const {

      CGAL_precondition(
        step == Intermediate_step::TREE_BOUNDARIES ||
        step == Intermediate_step::BUILDING_APPROXIMATE_BOUNDARIES ||
        step == Intermediate_step::BUILDING_BOUNDARIES);

      switch (step) {
        case Intermediate_step::TREE_BOUNDARIES: {
          return m_trees.get_tree_boundaries(output);
        }
        case Intermediate_step::BUILDING_APPROXIMATE_BOUNDARIES: {
          return m_buildings.get_approximate_boundaries(output);
        }
        case Intermediate_step::BUILDING_BOUNDARIES: {
          return m_buildings.get_building_boundaries(output);
        }
        default: return boost::none;
      }
    }

    /*!
      \brief returns mesh.

      This method returns data related to the intermediate steps of the
      algorithm, which can be stored as a mesh.

      \tparam VerticesOutputIterator
      must be an output iterator whose value type is `Kernel::Point_3`.

      \tparam FacesOutputIterator
      must be an output iterator whose value type is
      `std::pair<std::vector<std::size_t>, std::size_t>` or
      `std::pair<std::vector<std::size_t>, CGAL::Levels_of_detail::Visibility_label>` or
      `std::pair<std::vector<std::size_t>, std::size_t>`
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    mesh(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Intermediate_step step) const {

      CGAL_precondition(
        step == Intermediate_step::TREE_FOOTPRINTS ||
        step == Intermediate_step::EXTRUDED_TREE_BOUNDARIES ||
        step == Intermediate_step::EXTRUDED_TREE_FOOTPRINTS ||
        step == Intermediate_step::TREE_TRUNKS ||
        step == Intermediate_step::TREE_CROWNS ||
        step == Intermediate_step::BUILDING_PARTITIONING_2 ||
        step == Intermediate_step::BUILDING_FOOTPRINTS ||
        step == Intermediate_step::EXTRUDED_BUILDING_BOUNDARIES ||
        step == Intermediate_step::EXTRUDED_BUILDING_FOOTPRINTS ||
        step == Intermediate_step::APPROXIMATE_BUILDING_BOUNDS ||
        step == Intermediate_step::BUILDING_PARTITIONING_3 ||
        step == Intermediate_step::BUILDING_WALLS ||
        step == Intermediate_step::BUILDING_ROOFS);

      switch (step) {
        case Intermediate_step::TREE_FOOTPRINTS: {
          return m_trees.get_tree_footprints(vertices, faces);
        }
        case Intermediate_step::EXTRUDED_TREE_BOUNDARIES: {
          return m_trees.get_extruded_tree_boundaries(vertices, faces);
        }
        case Intermediate_step::EXTRUDED_TREE_FOOTPRINTS: {
          return m_trees.get_extruded_tree_footprints(vertices, faces);
        }
        case Intermediate_step::TREE_TRUNKS: {
          return m_trees.get_tree_trunks(vertices, faces);
        }
        case Intermediate_step::TREE_CROWNS: {
          return m_trees.get_tree_crowns(vertices, faces);
        }
        case Intermediate_step::BUILDING_PARTITIONING_2: {
          return m_buildings.get_partitioning_2(vertices, faces);
        }
        case Intermediate_step::BUILDING_FOOTPRINTS: {
          return m_buildings.get_building_footprints(vertices, faces);
        }
        case Intermediate_step::EXTRUDED_BUILDING_BOUNDARIES: {
          return m_buildings.get_extruded_building_boundaries(vertices, faces);
        }
        case Intermediate_step::EXTRUDED_BUILDING_FOOTPRINTS: {
          return m_buildings.get_extruded_building_footprints(vertices, faces);
        }
        case Intermediate_step::APPROXIMATE_BUILDING_BOUNDS: {
          return m_buildings.get_building_approximate_bounds(vertices, faces);
        }
        case Intermediate_step::BUILDING_PARTITIONING_3: {
          return m_buildings.get_building_partitioning_3(vertices, faces);
        }
        case Intermediate_step::BUILDING_WALLS: {
          return m_buildings.get_building_walls(vertices, faces);
        }
        case Intermediate_step::BUILDING_ROOFS: {
          return m_buildings.get_building_roofs(vertices, faces);
        }
        default: return boost::none;
      }
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    template<typename OutputIterator>
    boost::optional<OutputIterator> wire(
      OutputIterator output,
      const Wire_type step,
      const FT ground_precision) const {

      CGAL_precondition(
        step == Wire_type::PLANAR_GROUND_WIRE ||
        step == Wire_type::SMOOTH_GROUND_WIRE ||
        step == Wire_type::TREES_WIRE0 ||
        step == Wire_type::TREES_WIRE1 ||
        step == Wire_type::TREES_WIRE2 ||
        step == Wire_type::BUILDINGS_WIRE0 ||
        step == Wire_type::BUILDINGS_WIRE1 ||
        step == Wire_type::BUILDINGS_WIRE2 ||
        step == Wire_type::LOD_WIRE0 ||
        step == Wire_type::LOD_WIRE1 ||
        step == Wire_type::LOD_WIRE2);

      switch (step) {
        case Wire_type::PLANAR_GROUND_WIRE: {
          return m_ground.output_wire0(output);
        }
        case Wire_type::SMOOTH_GROUND_WIRE: {
          return m_ground.output_wire12(output, ground_precision);
        }
        case Wire_type::TREES_WIRE0: {
          return m_trees.output_wire0(output);
        }
        case Wire_type::TREES_WIRE1: {
          return m_trees.output_wire1(output);
        }
        case Wire_type::TREES_WIRE2: {
          return m_trees.output_wire2(output);
        }
        case Wire_type::BUILDINGS_WIRE0: {
          return m_buildings.output_wire0(output);
        }
        case Wire_type::BUILDINGS_WIRE1: {
          return m_buildings.output_wire1(output);
        }
        case Wire_type::BUILDINGS_WIRE2: {
          return m_buildings.output_wire2(output);
        }
        case Wire_type::LOD_WIRE0: {
          return m_lods.output_wire0(output);
        }
        case Wire_type::LOD_WIRE1: {
          return m_lods.output_wire1(output, ground_precision);
        }
        case Wire_type::LOD_WIRE2: {
          return m_lods.output_wire2(output, ground_precision);
        }
        default: return boost::none;
      }
    }
    /// \endcond

  private:
    const Input_range& m_input_range;
    const Item_map m_item_map;
    const Semantic_map m_semantic_map;
    const Parameters& m_parameters;
    const Visibility_map m_visibility_map;

    Data_structure m_data;
    const Ground m_ground;
    Buildings m_buildings;
    Trees m_trees;
    const LODs m_lods;

  }; // Levels_of_detail

} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_H
