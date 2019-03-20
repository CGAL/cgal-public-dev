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

#ifndef CGAL_LEVELS_OF_DETAIL_H
#define CGAL_LEVELS_OF_DETAIL_H

#include <CGAL/license/Levels_of_detail.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/property_map.h>

// Internal components.
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/Ground/Ground.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Buildings.h>
#include <CGAL/Levels_of_detail/internal/Trees/Trees.h>

#include <CGAL/Levels_of_detail/internal/Reconstruction/LOD0.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/LOD1.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/LOD2.h>

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
    range and value type is `GeomTraits::Point_3`.

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
    using Point_3 = typename Traits::Point_3;

    using Data_structure = internal::Data_structure<
    Traits, 
    Input_range, 
    Point_map, 
    Semantic_map, 
    Visibility_map>;

    using Ground = internal::Ground<Data_structure>;
    using Buildings = internal::Buildings<Data_structure>;
    using Trees = internal::Trees<Data_structure>;

    using LOD0 = internal::LOD0;
    using LOD1 = internal::LOD1;
    using LOD2 = internal::LOD2;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal structures.
      
      \param input_range
      an instance of `InputRange` with 3D points.

      \param point_map
      an instance of `PointMap` that maps an item from `input_range` 
      to `GeomTraits::Point_3`.

      \param semantic_map
      an instance of `SemanticMap` that maps an item from `input_range` 
      to `CGAL::Levels_of_detail::Semantic_label`.

      \param visibility_map
      an instance of `VisibilityMap` that maps an item from `input_range`
      to a value in the range [0,1].

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
    m_trees(m_data) { 

      CGAL_precondition(m_input_range.size() > 0);

      if (Verbose::value)
        std::cout << "Initializing LOD with:" << std::endl
          << "* " << m_data.ground_points().size() 
          << " ground point(s)" << std::endl
          << "* " << m_data.building_boundary_points().size() 
          << " building boundary point(s)" << std::endl
          << "* " << m_data.building_interior_points().size() 
          << " building interior point(s)" << std::endl
          << "* " << m_data.vegetation_points().size() 
          << " vegetation point(s)" << std::endl;
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
      \brief reconstructs a model of type `CGAL::Levels_of_detail::Reconstruction_type`.

      \param scale
      scale parameter

      \param noise_level
      noise level parameter

      \param ground_precision
      mas distance between input points and reconstructed ground

      \param reconstruction_type
      any of `CGAL::Levels_of_detail::Reconstruction_type`
    */
    void build(
      const FT scale, 
      const FT noise_level, 
      const FT ground_precision,
      const Reconstruction_type reconstruction_type) { 
      
      m_data.parameters.scale = scale;
      m_data.parameters.noise_level = noise_level;
      m_data.parameters.ground_precision = ground_precision;
      m_data.parameters.reconstruction_type = reconstruction_type;
    }

    /// @}

    /// \name Step by Step Generation
    /// @{

    /*!
      \brief computes a planar representation of the ground.

      The plane is estimated through principal component analysis on
      the points semantically labeled as `CGAL::Levels_of_detail::Semantic_label::GROUND`.
    */
    void compute_planar_ground() {
      m_ground.make_planar();
    }

    /*!
      \brief computes a smooth representation of the ground.
      
      The ground is represented as Delaunay triangulation with
      associated ground heights, which is computed upon the points
      semantically labeled as `CGAL::Levels_of_detail::Semantic_label::GROUND`.

      \param ground_precision
      max distance between input points and reconstructed ground
    */
    void compute_smooth_ground(const FT ground_precision) {
      
      m_data.parameters.ground_precision = ground_precision;
      m_ground.make_smooth(m_data.parameters.ground_precision);
    }

    /*!
      \brief detects building boundaries projected on the ground plane.

      This method:

      - computes the alpha shape of the points labeled as
        `CGAL::Levels_of_detail::Semantic_label::BUILDING_INTERIOR` (if any) and
        `CGAL::Levels_of_detail::Semantic_label::BUILDING_BOUNDARY` (if any) 
        and extracts the boundary points of this alpha shape;

      - downsamples this union of points using a regular grid;

      - detects line segments using the region growing approach.

      \warning `compute_planar_ground()` 
      should be called before calling this method

      \param alpha_shape_size
      alpha value from `CGAL::Alpha_shape_2`

      \param grid_cell_width
      fixed width of a regular grid cell

      \param region_growing_scale
      max distance from a query point such that all points within this distance are 
      assumed to be related to this query point

      \param region_growing_noise_level
      max distance from the boundary segment such that all points within this
      distance are assumed to be related to this segment

      \param region_growing_angle
      max angle in degrees between the point normal and the boundary segment normal

      \param region_growing_min_length
      min accepted length of each detected boundary segment
    */
    void detect_building_boundaries(
      const FT alpha_shape_size,
      const FT grid_cell_width,
      const FT region_growing_scale,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_min_length) {

        /*
        m_buildings.detect_boundaries(
          alpha_shape_size,
          grid_cell_width,
          region_growing_scale,
          region_growing_noise_level,
          region_growing_angle,
          region_growing_min_length); */
    }

    /*!
      \brief computes building footprints projected on the ground plane.

      This method:

      - creates the partitioning by extending initial building boundary
        segments until the defined number of intersections with other segments
        is reached;

      - applies the visibility computation that assignes to each polygon face 
        of the partitioning a visibility value in the range [0,1], where 0 
        means certainly outside and 1 means certainly inside;

      - tags subsets of all polygon faces with the visibility value >= 0.5
        that form separate buildings.

      \warning `detect_building_boundaries()` 
      should be called before calling this method

      \param kinetic_face_min_width
      min width of each detected polygon face

      \param kinetic_max_intersections
      max number of intersections between propagating segments

      \param min_faces_per_building
      min number of faces in the building footprint
    */
    void compute_building_footprints(
      const FT kinetic_face_min_width,
      const std::size_t kinetic_max_intersections,
      const std::size_t min_faces_per_building) {

      /*
      m_buildings.compute_footprints(
        kinetic_face_min_width,
        kinetic_max_intersections,
        min_faces_per_building); */
    }

    /*!
      \brief computes tree footprints projected on the ground plane.

      This method:

      - clusters all points labeled as 
      `CGAL::Levels_of_detail::Semantic_label::VEGETATION` that form potential trees;

      - estimates tree center point, radius, and height;

      - creates tree footprints.

      \warning `compute_planar_ground()` 
      should be called before calling this method

      \param grid_cell_width
      fixed width of a regular grid cell

      \param min_height
      min height of a tree

      \param min_radius
      min radius of a circle centered at the center of mass of a tree projected 
      on the ground plane

      \param min_faces_per_tree
      min number of faces in the tree footprint
    */
    void compute_tree_footprints(
      const FT grid_cell_width, 
      const FT min_height, 
      const FT min_radius,
      const std::size_t min_faces_per_tree) {
      
      /*
      m_trees.compute_footprints(
        grid_cell_width,
        min_height,
        min_radius, 
        min_faces_per_tree); */
    }

    /*!
      \brief extrudes the footprints to generate 3D buildings.
        
      The buildings are box models with a planar roof.

      \warning `compute_building_footprints()` 
      should be called before calling this method

      \param extrusion_type
      any of `CGAL::Levels_of_detail::Extrusion_type`
    */
    void extrude_building_footprints(const Extrusion_type extrusion_type) {
      // m_buildings.extrude_footprints(extrusion_type);
    }

    /*!
      \brief extrudes the footprints to generate 3D trees.
        
      The trees are cylinder models with a planar top.

      \warning `compute_tree_footprints()` 
      should be called before calling this method

      \param extrusion_type
      any of `CGAL::Levels_of_detail::Extrusion_type`
    */
    void extrude_tree_footprints(const Extrusion_type extrusion_type) {
      // m_trees.extrude_footprints(extrusion_type);
    }

    /*!
      \brief detects building roofs.

      This method:

      - detects chunks of 3D points that form planes using the region growing 
        approach on all points labeled as 
        `CGAL::Levels_of_detail::Semantic_label::BUILDING_INTERIOR`;

      - filters out all chunks, which do not fit to such criteria as 
        verticality, size, etc;

      - creates convex polygons, which approximate all left chunks.

      \warning `compute_building_footprints()` 
      should be called before calling this method

      \param region_growing_scale
      max distance from a query point such that all points within this distance are 
      assumed to be related to this query point

      \param region_growing_noise_level
      max distance from the roof such that all points within this
      distance are assumed to be related to this roof

      \param region_growing_angle
      max angle in degrees between the point normal and the roof normal

      \param region_growing_min_area
      min accepted area of each detected roof

      \param min_roof_size
      min size of each roof
    */
    void detect_building_roofs(
      const FT region_growing_scale,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_min_area,
      const FT min_roof_size) {

      /*
      m_buildings.detect_roofs(
        region_growing_scale,
        region_growing_noise_level,
        region_growing_angle,
        region_growing_min_area,
        min_roof_size); */
    }

    /*!
      \brief computes building roofs.

      This method:

      - creates the partitioning by extending initial building walls, roofs, and 
        ground represented as polygons until the defined number of 
        intersections with other polygons is reached;

      - applies the visibility computation that assignes to each polyhedral 
        facet of the partitioning a visibility value in the range [0,1], where 0 
        means certainly outside and 1 means certainly inside;

      - corrects the visibility estimations by applying a 3D graphcut;

      - extracts 3D polygons that represent exact roofs for each building.

      \warning `detect_building_roofs()` 
      should be called before calling this method

      \param kinetic_max_intersections
      max number of intersections between propagating polygons

      \param graph_cut_beta_3
      a graph cut precision parameter in the range [0,1], where 0 means
      keep all items and 1 means remove all of them
    */
    void compute_building_roofs(
      const std::size_t kinetic_max_intersections,
      const FT graph_cut_beta_3) {

      /*
      m_buildings.compute_roofs(
        kinetic_max_intersections,
        graph_cut_beta_3); */
    }

    /*!
      \brief computes tree models.

      This method:

      - creates tree icons;

      - fit these icons to all detected trees.

      \warning `compute_tree_footprints()` 
      should be called before calling this method

      \param precision
      max distance between points of the tree and its reconstructed model
    */
    void compute_tree_models(const FT precision) {
      // m_trees.compute_models(precision);
    }

    /// @}

    /// \name Output
    /// @{

    /*!
      \brief returns an estimated planar ground.

      \warning `compute_planar_ground()` 
      should be called before calling this method

      \tparam OutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.
    */
    template<typename OutputIterator>
    void output_ground_as_polygon(OutputIterator output) const {
      // m_ground.return_as_polygon(output);
    }

    /*!
      \brief returns an estimated smooth ground as triangle soup.

      \warning `compute_smooth_ground()` 
      should be called before calling this method

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is `cpp11::array<std::size_t, 3>`.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_ground_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      // m_ground.return_as_triangle_soup(output_vertices, output_faces);
    }

    /*!
      \brief returns points used for detecting building boundaries.

      All points are 3D points located on the estimated ground plane.

      \warning `detect_building_boundaries()` 
      should be called before calling this method
        
      \tparam OutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.
    */
    template<typename OutputIterator>
    void output_points_along_building_boundary(OutputIterator output) const {
      // m_buildings.return_boundary_points(output);
    }

    /*!
      \brief returns points along detected building walls.

      All points are 3D points located on the estimated ground plane.

      Detecting building boundaries creates a segmentation of the
      points: each point is associated to an index identifying a
      detected boundary segment or in other words a building wall.
      This index matches the order of segments given by 
      `output_building_boundaries_as_polylines()`. Points not associated to 
      any segment are given the index `std::size_t(-1)`.

      \warning `detect_building_boundaries()` 
      should be called before calling this method
        
      \tparam OutputIterator 
      must be an output iterator whose value type is 
      `std::pair<GeomTraits::Point_3, std::size_t>`.
    */
    template<typename OutputIterator>
    void output_points_along_building_walls(OutputIterator output) const {
      // m_buildings.return_wall_points(output);
    }

    /*!
      \brief returns polylines, which approximate building walls, or exact
      building boundary edges when available.

      All polylines are 3D segments located on the estimated ground plane.

      \warning `detect_building_boundaries()` 
      for approximate boundaries and `compute_building_footprints()` for exact 
      boundaries should be called before calling this method
        
      \tparam OutputIterator 
      must be an output iterator whose value type is `GeomTraits::Segment_3`.
    */
    template<typename OutputIterator>
    void output_building_boundaries_as_polylines(OutputIterator output) const {
      
      /*
      if (m_buildings.has_exact_boundaries())
        m_buildings.return_exact_boundary_edges(output);
      else
        m_buildings.return_approximate_boundary_edges(output); */
    }

    /*!
      \brief returns the partitionning based on boundary edges of all buildings.

      Each output face of the partitioning is a polygon.
        
      All vertices are 3D points located on the estimated ground plane.

      \warning `compute_building_footprints()` 
      should be called before calling this method

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair<std::vector<std::size_t>, CGAL::Levels_of_detail::Visibility_label>`,
      where the first item in the pair holds indices of the face vertices 
      and the second item is the visibility label.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_partitioning_as_polygon_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {

      // m_buildings.return_partitioning_2(output_vertices, output_faces);
    }

    /*!
      \brief returns footprints of all buildings as a triangle soup.
        
      Each triangle is associated to the index of the corresponding building.

      All vertices are 3D points located on the estimated ground
      plane or on the plane through the corresponding building height 
      if `extrude = true`.

      \warning `compute_building_footprints()` 
      should be called before calling this method and `extrude_building_footprints()` 
      if `extrude = true`

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair<cpp11::array<std::size_t, 3>, std::size_t>`,where the first item 
      in the pair holds indices of the face vertices and the second item is the 
      building index. All buildings are sorted by the index.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_footprints_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      const bool extruded = false) const {

      // m_buildings.return_footprints(output_vertices, output_faces, extruded);
    }

    /*!
      \brief returns clustered vegetation points used for detecting trees.

      \warning `compute_tree_footprints()` 
      should be called before calling this method
        
      \tparam OutputIterator 
      must be an output iterator whose value type is 
      `std::pair<GeomTraits::Point_3, std::size_t>`, where the first item in the 
      pair is the point and the second is the index of the corresponding cluster.
    */
    template<typename OutputIterator>
    void output_clustered_vegetation_points(OutputIterator output) const {
      // m_trees.return_clusters(output);
    }

    /*!
      \brief returns footprints of all trees as a triangle soup.
        
      Each triangle is associated to the index of the corresponding tree.

      All vertices are 3D points located on the estimated ground
      plane or on the plane through the corresponding tree height 
      if `extrude = true`.

      \warning `compute_tree_footprints()` 
      should be called before calling this method and `extrude_tree_footprints()` 
      if `extrude = true`

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair<cpp11::array<std::size_t, 3>, std::size_t>`, where the first 
      item in the pair holds indices of the face vertices and the second item 
      is the tree index. All trees are sorted by the index.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_tree_footprints_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      const bool extruded = false) const {
      // m_trees.return_footprints(output_vertices, output_faces, extruded);
    }

    /*!
      \brief returns polylines, which represent tree boundary edges.

      All polylines are 3D segments located on the estimated ground plane.

      \warning `compute_tree_footprints()` 
      should be called before calling this method
        
      \tparam OutputIterator 
      must be an output iterator whose value type is `GeomTraits::Segment_3`.
    */
    template<typename OutputIterator>
    void output_tree_boundaries_as_polylines(OutputIterator output) const {
      // m_trees.return_boundary_edges(output);
    }

    /*!
      \brief returns points along detected building roofs.

      Detecting roofs for each building creates a segmentation of the
      points: each point is associated to an index identifying a
      detected building roof.

      \warning `detect_building_roofs()` 
      should be called before calling this method.
        
      \tparam OutputIterator 
      must be an output iterator whose value type is 
      `std::tuple<GeomTraits::Point_3, std::size_t, std::size_t>`, where
      the first item in the tuple is the point, the second is the building index,
      and the third is the roof index.
    */
    template<typename OutputIterator>
    void output_points_along_building_roofs(OutputIterator output) const {
      // m_buildings.return_roof_points(output);
    }

    /*!
      \brief returns either approximate or exact building roofs.

      \warning `detect_building_roofs()` 
      for approximate roofs and `compute_building_roofs()` for exact roofs 
      should be called before calling this method

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::tuple<std::vector<std::size_t>, std::size_t, std::size_t>`, where 
      the first item in the tuple holds indices of the face vertices, 
      the second item is the building index, and the third item is 
      the roof index.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_roofs_as_polygon_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      
      /*
      if (m_buildings.has_exact_roofs())
        m_buildings.return_exact_roofs(output_vertices, output_faces);
      else
        m_buildings.return_approximate_roofs(output_vertices, output_faces); */
    }

    /*!
      \brief returns input to the partitioning algorithm.

      \warning `compute_building_roofs()` 
      should be called before calling this method

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair<std::vector<std::size_t>, std::size_t>`, where the first item 
      in the pair holds indices of the face vertices and the second item 
      is the building index.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_partitioning_in_3_as_polygon_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const { 
      // m_buildings.return_partitioning_input_3(output_vertices, output_faces);
    }

    /*!
      \brief returns output of the partitioning algorithm.

      \warning `compute_building_roofs()` 
      should be called before calling this method

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair<std::vector<std::size_t>, std::size_t>`, where the first item 
      in the pair holds indices of the face vertices and the second item 
      is the building index.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_partitioning_out_3_as_polygon_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      // m_buildings.return_partitioning_output_3(output_vertices, output_faces, false);
    }

    /*!
      \brief returns all building models as a polygon soup.

      \warning `compute_building_roofs()` 
      should be called before calling this method.

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair<std::vector<std::size_t>, std::size_t>`, where the first item 
      in the pair holds indices of the face vertices and the second item 
      is the building model index.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_models_as_polygon_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      // m_buildings.return_partitioning_output_3(output_vertices, output_faces, true);
    }

    /*!
      \brief returns all tree models as a triangle soup.
        
      Each triangle is associated to the index of the corresponding tree model.

      \warning `compute_tree_models()` 
      should be called before calling this method

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair<cpp11::array<std::size_t, 3>, std::size_t>`, where the first item 
      in the pair holds indices of the face vertices and the second item 
      is the tree model index.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_tree_models_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      // m_trees.return_models(output_vertices, output_faces);
    }

    /*!
      \brief returns LOD.

      \warning `build()` 
      should be called before calling this method

      \tparam VerticesOutputIterator
      must be an output iterator whose value type is `GeomTraits::Point_3`. 

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair< cpp11::array<std::size_t, 3>, Urban_object_type>`, where the 
      first item in the pair holds indices of the face vertices and second
      item is the type of the urban object this face belongs to.

      \pre 
      `lod_type == Reconstruction_type::LOD0 ||
       lod_type == Reconstruction_type::LOD1 ||
       lod_type == Reconstruction_type::LOD2`
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_LOD_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      const Reconstruction_type lod_type,
      const FT ground_precision = FT(0)) {
      
      /*
      CGAL_precondition(
        lod_type == Reconstruction_type::LOD0 ||
        lod_type == Reconstruction_type::LOD1 ||
        lod_type == Reconstruction_type::LOD2);

      switch (lod_type) {

        case Reconstruction_type::LOD0 : {
          
          CGAL_precondition(ground_precision == FT(0));
          LOD0 lod0(m_data_structure);

          lod0.reconstruct();
          lod0.output_as_triangle_soup(output_vertices, output_faces);

          return;
        }

        case Reconstruction_type::LOD1 : {

          CGAL_precondition(ground_precision != FT(0));
          LOD1 lod1(m_data_structure, ground_precision);

          lod1.reconstruct();
          lod1.output_as_triangle_soup(output_vertices, output_faces);

          return;
        }

        case Reconstruction_type::LOD2 : {

          CGAL_precondition(ground_precision != FT(0));
          LOD2 lod2(m_data_structure, ground_precision);

          lod2.reconstruct();
          lod2.output_as_triangle_soup(output_vertices, output_faces);

          return;
        }

        default: return;
      } */
    }

    /// @}

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Semantic_map m_semantic_map;
    const Visibility_map m_visibility_map;

    Data_structure m_data;
    Ground m_ground;
    Buildings m_buildings;
    Trees m_trees;

  }; // Levels_of_detail

} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_H
