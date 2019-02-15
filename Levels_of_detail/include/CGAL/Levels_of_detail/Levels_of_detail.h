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

// LOD includes.
#include <CGAL/Levels_of_detail/enumerations.h>
#include <CGAL/Levels_of_detail/property_maps.h>
#include <CGAL/Levels_of_detail/internal/utilities.h>
#include <CGAL/Levels_of_detail/internal/Data_structure.h>

// Internal components.
#include <CGAL/Levels_of_detail/internal/Ground.h>
#include <CGAL/Levels_of_detail/internal/Buildings.h>
#include <CGAL/Levels_of_detail/internal/Vegetation.h>

namespace CGAL {
namespace Levels_of_detail {

  /*!
    \ingroup PkgLevelsOfDetailRef

    \brief The Levels Of Detail algorithm, constructs Levels Of Detail (LOD) from an input point cloud.

    \tparam GeomTraits A model of \cgal `Kernel`.

    \tparam InputRange A range with points. 
    A model of `ConstRange`. The value type of its iterator is the key type of `PointMap`.

    \tparam PointMap Returns a point from `InputRange`. 
    A model of `ReadablePropertyMap` whose key type is the value type of the iterator of `InputRange` 
    and value type is `CGAL::Point_3`.

    \tparam SemanticMap Maps a point from `InputRange` to a semantic class from `SemanticLabel`. 
    A model of `ReadablePropertyMap` whose key type is the value type of the iterator of `InputRange` 
    and value type is `Semantic_label`.

    \tparam VisibilityMap Maps a point from `InputRange` to a visibility value in the range [0,1].
    A model of `ReadablePropertyMap` whose key type is the value type of the iterator of `InputRange` 
    and value type is `GeomTraits::FT`.

    \tparam Verbose Use if you want to print extra information about execution of the algorithm.
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

    /// \name Types
    /// @{
      
    using Traits = GeomTraits;
    ///< A traits class with geometric constructors and predicates.

    using Input_range = InputRange;
    ///< A point range in 3D.

    using Point_map = PointMap;
    ///< A map that returns a point from `Input_range`.

    using Semantic_map = SemanticMap;
    ///< A map that returns a semantic class from `Semantic_label` for each point in `Input_range`.

    using Visibility_map = VisibilityMap;
    ///< A map that returns a visibility value [0,1] for each point from `Input_range`.
      
    /// @}

    /// \cond SKIP_IN_MANUAL

    using FT = typename Traits::FT;

    using Data_structure = internal::Data_structure<
    Traits, 
    Input_range, 
    Point_map, 
    Semantic_map, 
    Visibility_map>;

    using Ground = internal::Ground<Data_structure>;
    using Buildings = internal::Buildings<Data_structure>;
    using Vegetation = internal::Vegetation<Data_structure>;

    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief Initializes data structures for computing Levels Of Detail, 
      given an input range with 3D points, a point, semantic, and visibility map.
    */
    Levels_of_detail(
      const Input_range& input_range,
      Point_map point_map,
      Semantic_map semantic_map,
      Visibility_map visibility_map = VisibilityMap()) : 
    m_data_structure(
      input_range, 
      point_map, 
      semantic_map, 
      visibility_map,
      Verbose::value ? true : false),
    m_ground(m_data_structure),
    m_buildings(m_data_structure),
    m_vegetation(m_data_structure) { 

      if (Verbose::value)
        std::cout << "Initializing LOD with:" << std::endl
          << "* " << m_data_structure.ground_points().size() 
          << " ground point(s)" << std::endl
          << "* " << m_data_structure.building_boundary_points().size() 
          << " building boundary point(s)" << std::endl
          << "* " << m_data_structure.building_interior_points().size() 
          << " building interior point(s)" << std::endl
          << "* " << m_data_structure.vegetation_points().size() 
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

    void build(
      const FT scale, 
      const FT noise_level, 
      const FT ground_precision,
      const Reconstruction_type reconstruction_type) { 

    }

    /// @}

    /// \name Step by Step Generation
    /// @{

    /*!
      \brief Computes a planar representation of the ground.

      The plane is estimated through principal component analysis on
      the points semantically labeled as `Semantic_label::GROUND`.
    */
    void compute_planar_ground() {
      m_ground.make_planar();
    }

    /*!
      \brief Computes a smooth representation of the ground.
      
      The ground is represented as Delaunay triangulation with
      associated ground heights, which is computed upon the points
      semantically labeled as `Semantic_label::GROUND`.
    */
    void compute_smooth_ground() {
      m_ground.make_smooth();
    }

    /*!
      \brief Detects building boundaries projected on the ground plane.

      This method:

      - computes the alpha shape of the points labeled as
        `BUILDING_INTERIOR` (if any) and `BUILDING_BOUNDARY` (if any) 
        and extracts the boundary points of this alpha shape;

      - downsamples this union of points using a regular grid;

      - detects line segments using the region growing approach;

      \warning `compute_planar_ground()` should be called before
      calling this method.
    */
    void detect_building_boundaries(
      const FT alpha_shape_size,
      const FT grid_cell_width,
      const FT region_growing_search_size,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_min_length) {
        
        m_buildings.detect_boundaries(
          alpha_shape_size,
          grid_cell_width,
          region_growing_search_size,
          region_growing_noise_level,
          region_growing_angle,
          region_growing_min_length);
    }

    /*!
      \brief Computes building footprints projected on the ground plane.

      This method:

      - creates the partitioning by extending initial building boundary
        segments until the defined number of intersections with other segments
        is reached;

      - applies the visibility computation that assignes to each polygon face 
        of the partitioning a visibility value in the range [0, 1], where 0 
        means certainly outside and 1 means certainly inside.

      - tags subsets of all polygon faces with the visibility value >= 0.5
        that form separate buildings.

      \warning `detect_building_boundaries()` should be called
      before calling this method.
    */
    void compute_building_footprints(
      const FT kinetic_min_face_width,
      const std::size_t kinetic_max_intersections,
      const std::size_t min_faces_per_building) {

        m_buildings.compute_footprints(
          kinetic_min_face_width,
          kinetic_max_intersections,
          min_faces_per_building);
    }

    /*!
      \brief Computes tree footprints projected on the ground plane.

      This method:

      - clusters vegetation points that form potential trees;

      - estimates tree center point, radius, and height;

      - creates tree footprints.

      \warning `compute_planar_ground()` should be called 
      before calling this method.
    */
    void compute_tree_footprints(
      const FT grid_cell_width, 
      const FT min_height, 
      const FT min_radius,
      const std::size_t min_faces_per_tree) {
      
      m_vegetation.compute_tree_footprints(
        grid_cell_width,
        min_height,
        min_radius, 
        min_faces_per_tree);
    }

    /*!
      \brief Extrudes the footprints to generate 3D buildings.
        
      The buildings are shoebox models with a planar roof.

      \warning `compute_building_footprints()` should be called before
      calling this method.
    */
    void extrude_building_footprints(const Extrusion_type extrusion_type) {
      m_buildings.extrude_footprints(extrusion_type);
    }

    /*!
      \brief Extrudes the footprints to generate 3D trees.
        
      The trees are shoebox models with a planar top.

      \warning `compute_tree_footprints()` should be called before
      calling this method.
    */
    void extrude_tree_footprints(const Extrusion_type extrusion_type) {
      m_vegetation.extrude_tree_footprints(extrusion_type);
    }

    /*!
      \brief Detects building roofs.

      This method:

      - detects planar point subsets, using the region growing approach, on 
        all points labeled as `BUILDING_INTERIOR`;

      - filters out all subsets, which do not fit to such criteria as 
        verticality, size, etc;

      - creates convex polygons, which approximate all the rest point subsets.

      \warning `compute_building_footprints()` should be called 
      before calling this method.
    */
    void detect_building_roofs(
      const FT region_growing_search_size,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_min_area,
      const FT min_size) {
      
      m_buildings.detect_roofs(
        region_growing_search_size,
        region_growing_noise_level,
        region_growing_angle,
        region_growing_min_area,
        min_size);
    }

    /// @}

    /// \name Output
    /// @{

    /*!
      \brief Returns an estimated planar ground.

      \warning `compute_planar_ground()` should be called
      before calling this method.

      \tparam OutputIterator is a model of `OutputIterator`
      that holds `CGAL::Point_3` objects.

      \param output an iterator with polygon vertices given as 3D points.
    */
    template<typename OutputIterator>
    void output_ground_as_polygon(OutputIterator output) const {
      m_ground.return_as_polygon(output);
    }

    /*!
      \brief Returns an estimated smooth ground as triangle soup.

      \warning `compute_smooth_ground()` should be called
      before calling this method.

      \tparam VerticesOutputIterator is a model of `OutputIterator`
      that holds `Point_3` objects.

      \tparam FacesOutputIterator is a model of `OutputIterator`
      that holds `cpp11::array<std::size_t, 3>` objects.

      \param output_vertices an iterator with all vertices of the triangle soup.

      \param output_faces an iterator with all faces of the triangle soup
      given as arrays of indices in `output_vertices`.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_ground_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      
      m_ground.return_as_triangle_soup(output_vertices, output_faces);
    }

    /*!
      \brief Returns points used for detecting building boundaries.

      All points are 3D points located on the estimated ground
      plane (see `ground_plane()`).

      \warning `detect_building_boundaries()` should be called
      before calling this method.
        
      \tparam OutputIterator is a model of `OutputIterator`
      that holds `CGAL::Point_3` objects.

      \param output an iterator with 3D points.
    */
    template<typename OutputIterator>
    void output_points_along_building_boundary(OutputIterator output) const {
      m_buildings.return_boundary_points(output);
    }

    /*!
      \brief Returns points along detected building walls.

      All points are 3D points located on the estimated ground
      plane (see `ground_plane()`).

      Detecting building boundaries creates a segmentation of the
      points: each point is associated to an index identifying a
      detected boundary segment or in other words a building wall.
      This index matches the order of segments given by 
      `output_building_boundaries_as_polylines()`. Points not associated to 
      any segment are given the index `-1`.

      \warning `detect_building_boundaries()` should be called
      before calling this method.
        
      \tparam OutputIterator is a model of `OutputIterator`
      that holds `std::pair<Point_3, long>` objects.

      \param output an iterator with points and assigned to them ids of
      the detected building walls.
    */
    template<typename OutputIterator>
    void output_points_along_building_walls(OutputIterator output) const {
      m_buildings.return_wall_points(output);
    }

    /*!
      \brief Returns polylines, which approximate building walls, or exact
      building boundary edges when available.

      All polylines are 3D segments located on the estimated ground
      plane (see `ground_plane()`).

      \warning `detect_building_boundaries()` for approximate boundaries and
      `compute_building_footprints()` for exact boundaries should be called
      before calling this method.
        
      \tparam OutputIterator is a model of `OutputIterator`
      that holds `CGAL::Segment_3` objects.

      \param output an iterator with 3D segments.
    */
    template<typename OutputIterator>
    void output_building_boundaries_as_polylines(OutputIterator output) const {
      
      if (m_buildings.has_exact_boundaries())
        m_buildings.return_exact_boundary_edges(output);
      else
        m_buildings.return_approximate_boundary_edges(output);
    }

    /*!
      \brief Returns the partitionning based on boundary edges of all buildings.

      Each output face of the partitioning is a polygon.
        
      All vertices are 3D points located on the estimated ground
      plane (see `ground_plane()`).

      \warning `compute_building_footprints()` should be called before
      calling this method.

      \tparam VerticesOutputIterator is a model of `OutputIterator`
      that holds `Point_3` objects.

      \tparam FacesOutputIterator is a model of `OutputIterator`
      that holds an `std::pair<std::vector<std::size_t>, Levels_of_detail::Visibility_label>` 
      objects, where the first item in the pair holds indices of the face vertices 
      and the second item is the visibility label.

      \param output_vertices an iterator with all vertices of the polygon soup.
      
      \param output_faces an iterator with all faces of the polygon soup
      given as vectors of indices in `output_vertices` and the corresponding 
      visibility labels.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_partitioning_as_polygon_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {

      m_buildings.return_partitioning(output_vertices, output_faces);
    }

    /*!
      \brief Returns footprints of all buildings as a triangle soup.
        
      Each triangle is associated to the index of the corresponding
      building.

      All vertices are 3D points located on the estimated ground
      plane (see `ground_plane()`) or on the plane through the corresponding 
      building height if `extrude = true`.

      \warning `compute_building_footprints()` should be called before
      calling this method and `extrude_building_footprints()` if `extrude = true`.

      \tparam VerticesOutputIterator is a model of `OutputIterator`
      that holds `Point_3` objects.

      \tparam FacesOutputIterator is a model of `OutputIterator`
      that holds `std::pair<cpp11::array<std::size_t, 3>, std::size_t>` objects,
      where the first item in the pair holds indices of the face vertices and
      the second item is the building index. All buildings are sorted by the index.

      \param output_vertices an iterator with all vertices of the triangle soup.

      \param output_faces an iterator with all faces of the triangle soup
      given as arrays of indices in `output_vertices` and the corresponding
      building indices.

      \param extruded should be false, which is default, if no extrusion was 
      made prior to calling this function, otherwise can be true to output each 
      footprint brought to the height of the corresponding building.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_footprints_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      const bool extruded = false) const {

      m_buildings.return_footprints(output_vertices, output_faces, extruded);
    }

    /*!
      \brief Returns clustered vegetation points used for detecting trees.

      All points are 3D points located on the estimated ground
      plane (see `ground_plane()`).

      \warning `compute_tree_footprints()` should be called
      before calling this method.
        
      \tparam OutputIterator is a model of `OutputIterator`
      that holds `std::pair<CGAL::Point_3, std::size_t>` objects.

      \param output an iterator with 3D points and the corresponding
      cluster indices.
    */
    template<typename OutputIterator>
    void output_clustered_vegetation_points(OutputIterator output) const {
      m_vegetation.return_clustered_points(output);
    }

    /*!
      \brief Returns footprints of all trees as a triangle soup.
        
      Each triangle is associated to the index of the corresponding
      tree.

      All vertices are 3D points located on the estimated ground
      plane (see `ground_plane()`) or on the plane through the corresponding 
      tree height if `extrude = true`.

      \warning `compute_tree_footprints()` should be called before
      calling this method and `extrude_tree_footprints()` if `extrude = true`.

      \tparam VerticesOutputIterator is a model of `OutputIterator`
      that holds `Point_3` objects.

      \tparam FacesOutputIterator is a model of `OutputIterator`
      that holds `std::pair<cpp11::array<std::size_t, 3>, std::size_t>` objects,
      where the first item in the pair holds indices of the face vertices and
      the second item is the tree index. All trees are sorted by the index.

      \param output_vertices an iterator with all vertices of the triangle soup.

      \param output_faces an iterator with all faces of the triangle soup
      given as arrays of indices in `output_vertices` and the corresponding
      tree indices.

      \param extruded should be false, which is default, if no extrusion was 
      made prior to calling this function, otherwise can be true to output each 
      footprint brought to the height of the corresponding tree.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_tree_footprints_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      const bool extruded = false) const {

      m_vegetation.return_tree_footprints(
        output_vertices, output_faces, extruded);
    }

    /*!
      \brief Returns polylines, which represent tree boundary edges.

      All polylines are 3D segments located on the estimated ground
      plane (see `ground_plane()`).

      \warning `compute_tree_footprints()` should be called
      before calling this method.
        
      \tparam OutputIterator is a model of `OutputIterator`
      that holds `CGAL::Segment_3` objects.

      \param output an iterator with 3D segments.
    */
    template<typename OutputIterator>
    void output_tree_boundaries_as_polylines(OutputIterator output) const {
      m_vegetation.return_tree_boundary_edges(output);
    }

    /*!
      \brief Returns points along detected building roofs.

      Detecting roofs for each building creates a segmentation of the
      points: each point is associated to an index identifying a
      detected building roof.

      \warning `detect_building_roofs()` should be called
      before calling this method.
        
      \tparam OutputIterator is a model of `OutputIterator`
      that holds `std::tuple<Point_3, long, long>` objects.

      \param output an iterator with points, assigned to them ids of
      the corresponding buildings, and detected roofs.
    */
    template<typename OutputIterator>
    void output_points_along_building_roofs(OutputIterator output) const {
      m_buildings.return_roof_points(output);
    }

    /*!
      \brief Returns either approximate or exact building roofs.

      \warning `detect_building_roofs()` for approximate roofs and
      `compute_building_roofs()` for exact roofs should be called
      before calling this method.

      \tparam VerticesOutputIterator is a model of `OutputIterator`
      that holds `Point_3` objects.

      \tparam FacesOutputIterator is a model of `OutputIterator`
      that holds `std::tuple<std::vector<std::size_t>, std::size_t, std::size_t>` 
      objects, where the first item in the tuple holds indices of the face 
      vertices, the second item is the building index, and the third item is 
      the roof index.

      \param output_vertices an iterator with all vertices of the polygon soup.

      \param output_faces an iterator with all faces of the polygon soup
      given as vectors of indices in `output_vertices`, the corresponding
      building indices, and roof indices.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_roofs_as_polygon_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      
      if (m_buildings.has_exact_roofs())
        m_buildings.return_exact_roofs(output_vertices, output_faces);
      else
        m_buildings.return_approximate_roofs(output_vertices, output_faces);
    }

    /// @}

  private:
    Data_structure m_data_structure;
    Ground m_ground;
    Buildings m_buildings;
    Vegetation m_vegetation;

  }; // end of class

} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_H
