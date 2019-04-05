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
#include <CGAL/Levels_of_detail/internal/Trees/Trees.h>
#include <CGAL/Levels_of_detail/internal/Ground/Ground.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Buildings.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/Reconstruction.h>

namespace CGAL {
namespace Levels_of_detail {

  /*!
    \ingroup PkgLevelsOfDetailRef

    \brief Given a point cloud, reconstructs its model with Levels Of Detail (LOD).

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam PointMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `Kernel::Point_3`.

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
    typename PointMap,
    typename SemanticMap,
    typename VisibilityMap = Visibility_from_semantic_map<SemanticMap>,
    typename Verbose = CGAL::Tag_false>
	class Levels_of_detail {

	public:

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
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
    Point_map, 
    Semantic_map, 
    Visibility_map>;

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
      an instance of `InputRange` with 3D points

      \param point_map
      an instance of `PointMap` that maps an item from `input_range` 
      to `Kernel::Point_3`

      \param semantic_map
      an instance of `SemanticMap` that maps an item from `input_range` 
      to `CGAL::Levels_of_detail::Semantic_label`

      \param visibility_map
      an instance of `VisibilityMap` that maps an item from `input_range`
      to a value in the range [0,1]

      \pre `input_range.size() > 0`
    */
    Levels_of_detail(
      const InputRange& input_range,
      const PointMap point_map,
      const SemanticMap semantic_map,
      const VisibilityMap visibility_map = VisibilityMap()) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_semantic_map(semantic_map),
    m_visibility_map(visibility_map), 
    m_data(
      m_input_range, 
      m_point_map, 
      m_semantic_map, 
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
    /// \endcond

    /// \name Complete Generation
    /// @{

    /*!
      \brief builds all available urban objects.

      This method:

      - builds all trees;

      - builds all buildings.

      \param scale
      scale parameter

      \param noise_level
      noise level parameter
    */
    void build(
      const FT scale, 
      const FT noise_level) { 
      
      m_data.parameters.scale = scale;
      m_data.parameters.noise_level = noise_level;
      m_data.parameters.update_dependent();
      
      m_trees.make_trees();
      m_buildings.make_buildings();
    }

    /*!
      \brief builds all trees.

      \param scale
      scale parameter

      \param noise_level
      noise level parameter
    */
    void build_trees(
      const FT scale, 
      const FT noise_level) { 
      
      m_data.parameters.scale = scale;
      m_data.parameters.noise_level = noise_level;
      m_data.parameters.update_dependent();
      
      m_trees.make_trees();
    }

    /*!
      \brief builds all buildings.

      \param scale
      scale parameter

      \param noise_level
      noise level parameter
    */
    void build_buildings(
      const FT scale, 
      const FT noise_level) { 
      
      m_data.parameters.scale = scale;
      m_data.parameters.noise_level = noise_level;
      m_data.parameters.update_dependent();
      
      m_trees.make_buildings();
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

      \param scale
      scale parameter

      \param noise_level
      noise level parameter

      \param cluster_scale
      cluster scale parameter
    */
    void initialize_trees(
      const FT scale, 
      const FT noise_level,
      const FT cluster_scale) {
      
      m_data.parameters.scale = scale;
      m_data.parameters.noise_level = noise_level;
      m_data.parameters.update_dependent();
      m_data.parameters.trees.cluster_scale = cluster_scale;

      m_trees.initialize();
    }

    /*!
      \brief computes tree footprints.

      This method:

      - extracts points, which form potential trees;

      - estimates each tree center, radius, and height;

      - creates tree footprints.

      \warning `initialize_trees()` should be called before 
      calling this method

      \param grid_cell_width_2
      fixed width of a cell in a 2D regular grid

      \param min_height
      min height of a tree

      \param min_radius_2
      min radius of a circle that approximates the boundary of all tree points
      projected in 2D

      \param min_faces_per_footprint
      min number of faces in the tree footprint
    */
    void compute_tree_footprints(
      const FT grid_cell_width_2, 
      const FT min_height, 
      const FT min_radius_2,
      const std::size_t min_faces_per_footprint) {
      
      m_data.parameters.trees.grid_cell_width_2 = grid_cell_width_2;
      m_data.parameters.trees.min_height = min_height;
      m_data.parameters.trees.min_radius_2 = min_radius_2;
      m_data.parameters.trees.min_faces_per_footprint = min_faces_per_footprint;
      
      m_trees.compute_footprints();
    }

    /*!
      \brief extrudes tree footprints.

      The trees, after extrusion, are cylinder models with a planar top.

      \warning `compute_tree_footprints()` should be called before 
      calling this method

      \param extrusion_type
      any of `CGAL::Levels_of_detail::Extrusion_type`
    */
    void extrude_tree_footprints(
      const Extrusion_type extrusion_type) {
      
      m_data.parameters.trees.extrusion_type = extrusion_type;
      m_trees.extrude_footprints();
    }

    /*!
      \brief computes tree crowns.

      This method:

      - creates crown icons;

      - fits these icons to trees.

      \warning `extrude_tree_footprints()` should be called before 
      calling this method
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

      \param scale
      scale parameter

      \param noise_level
      noise level parameter

      \param cluster_scale
      cluster scale parameter
    */
    void initialize_buildings(
      const FT scale, 
      const FT noise_level,
      const FT cluster_scale) {
      
      m_data.parameters.scale = scale;
      m_data.parameters.noise_level = noise_level;
      m_data.parameters.update_dependent();
      m_data.parameters.buildings.cluster_scale = cluster_scale;
      
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

      \warning `initialize_buildings()` should be called before 
      calling this method

      \param alpha_shape_size_2
      alpha value from `CGAL::Alpha_shape_2`

      \param grid_cell_width_2
      fixed width of a cell in a 2D regular grid

      \param region_growing_scale_2
      region growing scale

      \param region_growing_noise_level_2
      region growing noise level

      \param region_growing_angle_2
      region growing angle

      \param region_growing_min_length_2
      min accepted length of each detected boundary segment
    */
    void detect_building_boundaries(
      const FT alpha_shape_size_2,
      const FT grid_cell_width_2,
      const FT region_growing_scale_2,
      const FT region_growing_noise_level_2,
      const FT region_growing_angle_2,
      const FT region_growing_min_length_2) {

      m_data.parameters.buildings.alpha_shape_size_2 = alpha_shape_size_2;
      m_data.parameters.buildings.grid_cell_width_2 = grid_cell_width_2;
      m_data.parameters.buildings.region_growing_scale_2 = region_growing_scale_2;
      m_data.parameters.buildings.region_growing_noise_level_2 = region_growing_noise_level_2;
      m_data.parameters.buildings.region_growing_angle_2 = region_growing_angle_2;
      m_data.parameters.buildings.region_growing_min_length_2 = region_growing_min_length_2;

      m_buildings.detect_boundaries();
    }

    /*!
      \brief computes building footprints.

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

      \warning `detect_building_boundaries()` should be called before 
      calling this method

      \param kinetic_min_face_width_2
      min width of each detected polygon face

      \param kinetic_max_intersections_2
      max number of intersections between propagating segments

      \param min_faces_per_footprint
      min number of faces in the building footprint

      \param visibility_scale_2
      visibility scale 2

      \param graphcut_beta_2
      a graph cut precision parameter in the range [0,1], where 0 means
      keep items and 1 means remove them
    */
    void compute_building_footprints(
      const FT kinetic_min_face_width_2,
      const std::size_t kinetic_max_intersections_2,
      const std::size_t min_faces_per_footprint,
      const FT visibility_scale_2,
      const FT graphcut_beta_2) {

      m_data.parameters.buildings.kinetic_min_face_width_2 = kinetic_min_face_width_2;
      m_data.parameters.buildings.kinetic_max_intersections_2 = kinetic_max_intersections_2;
      m_data.parameters.buildings.min_faces_per_footprint = min_faces_per_footprint;
      m_data.parameters.buildings.visibility_scale_2 = visibility_scale_2;
      m_data.parameters.buildings.graphcut_beta_2 = graphcut_beta_2;

      m_buildings.compute_footprints();
    }

    /*!
      \brief extrudes building footprints.
        
      The buildings, after extrusion, are box models with a planar top.

      \warning `compute_building_footprints()` should be called before 
      calling this method

      \param extrusion_type
      any of `CGAL::Levels_of_detail::Extrusion_type`
    */
    void extrude_building_footprints(
      const Extrusion_type extrusion_type) {
      
      m_data.parameters.buildings.extrusion_type = extrusion_type;
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

      \warning `extrude_building_footprints()` should be called before 
      calling this method

      \param region_growing_scale_3
      region growing scale

      \param region_growing_noise_level_3
      region growing noise level

      \param region_growing_angle_3
      region growing angle

      \param region_growing_min_area_3
      min accepted area of each detected roof polygon 

      \param region_growing_distance_to_line_3
      distance to line
    */
    void detect_building_roofs(
      const FT region_growing_scale_3,
      const FT region_growing_noise_level_3,
      const FT region_growing_angle_3,
      const FT region_growing_min_area_3,
      const FT region_growing_distance_to_line_3) {

      m_data.parameters.buildings.region_growing_scale_3 = region_growing_scale_3;
      m_data.parameters.buildings.region_growing_noise_level_3 = region_growing_noise_level_3;
      m_data.parameters.buildings.region_growing_angle_3 = region_growing_angle_3;
      m_data.parameters.buildings.region_growing_min_area_3 = region_growing_min_area_3;
      m_data.parameters.buildings.region_growing_distance_to_line_3 = region_growing_distance_to_line_3;

      m_buildings.detect_roofs();
    }

    /*!
      \brief computes building roofs.

      This method:

      - creates the partitioning by extending initial building walls, roofs, and 
        base represented as polygons until the defined number of 
        intersections with other polygons is reached;

      - applies the visibility computation that assignes to each polyhedral 
        facet of the partitioning a visibility value in the range [0,1], where 0 
        means certainly outside and 1 means certainly inside;

      - corrects the visibility estimations by applying a graphcut;

      - extracts polygons, which represent building roofs.

      \warning `detect_building_roofs()` should be called before 
      calling this method

      \param kinetic_max_intersections_3
      max number of intersections between propagating polygons

      \param visibility_scale_3
      visibility scale 3

      \param graphcut_beta_3
      a graph cut precision parameter in the range [0,1], where 0 means
      keep items and 1 means remove them
    */
    void compute_building_roofs(
      const std::size_t kinetic_max_intersections_3,
      const FT visibility_scale_3,
      const FT graphcut_beta_3) {

      m_data.parameters.buildings.kinetic_max_intersections_3 = kinetic_max_intersections_3;
      m_data.parameters.buildings.visibility_scale_3 = visibility_scale_3;
      m_data.parameters.buildings.graphcut_beta_3 = graphcut_beta_3;

      m_buildings.compute_roofs();
    }

    /// @}

    /// \name Access
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

      \param ground_type
      either `Reconstruction_type::PLANAR_GROUND` or `Reconstruction_type::SMOOTH_GROUND`

      \param ground_precision
      max distance between input points and a reconstructed ground

      \pre `ground_type == Reconstruction_type::PLANAR_GROUND ||
            ground_type == Reconstruction_type::SMOOTH_GROUND`

      \return a pair of the past-the-end iterators
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    ground(
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const Reconstruction_type ground_type,
      const FT ground_precision) const {
      
      CGAL_precondition(
        ground_type == Reconstruction_type::PLANAR_GROUND ||
        ground_type == Reconstruction_type::SMOOTH_GROUND);
      
      return m_ground.output_ground(vertices, faces, ground_type, ground_precision);
    }

    /*!
      \brief returns trees as triangle soup.

      This method returns different levels of detail for trees.

      \warning the corresponding tree-related functions above 
      should be called.

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `Kernel::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair<std::vector<std::size_t>, CGAL::Levels_of_detail::Urban_object_type::TREE>`.

      \param vertices
      an output iterator with vertices

      \param faces
      an output iterator with faces

      \param lod_type
      either `Reconstruction_type::TREES0` or `Reconstruction_type::TREES1` or 
      `Reconstruction_type::TREES2`

      \pre `lod_type == Reconstruction_type::TREES0 ||
            lod_type == Reconstruction_type::TREES1 ||
            lod_type == Reconstruction_type::TREES2`

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

      \warning the corresponding building-related functions above 
      should be called.

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `Kernel::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair<std::vector<std::size_t>, CGAL::Levels_of_detail::Urban_object_type::BUILDING>`.

      \param vertices
      an output iterator with vertices

      \param faces
      an output iterator with faces

      \param lod_type
      either `Reconstruction_type::BUILDINGS0` or `Reconstruction_type::BUILDINGS1` or 
      `Reconstruction_type::BUILDINGS2`

      \pre `lod_type == Reconstruction_type::BUILDINGS0 ||
            lod_type == Reconstruction_type::BUILDINGS1 ||
            lod_type == Reconstruction_type::BUILDINGS2`

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

      \pre `lod_type == Reconstruction_type::LOD0 ||
            lod_type == Reconstruction_type::LOD1 ||
            lod_type == Reconstruction_type::LOD2`

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

    /// \name Intermediate Steps
    /// @{

    /*!
      \brief returns a point set.
        
      This method returns data related to the intermediate steps of the 
      algorithm, which can be stored as a point set.

      \tparam OutputIterator 
      must be an output iterator whose value type is 
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> points(
      OutputIterator output,
      const Intermediate_step step) const {

      CGAL_precondition(
        step == Intermediate_step::TREE_CLUSTERS ||
        step == Intermediate_step::TREE_POINTS ||
        step == Intermediate_step::BUILDING_CLUSTERS ||
        step == Intermediate_step::BUILDING_BOUNDARY_POINTS ||
        step == Intermediate_step::BUILDING_WALL_POINTS ||
        step == Intermediate_step::BUILDING_POINTS ||
        step == Intermediate_step::BUILDING_ROOF_POINTS);

      switch (step) {
        case Intermediate_step::TREE_CLUSTERS: {
          return m_trees.get_tree_clusters(output);
        }
        case Intermediate_step::TREE_POINTS: {
          return m_trees.get_tree_points(output);
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

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Semantic_map m_semantic_map;
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
