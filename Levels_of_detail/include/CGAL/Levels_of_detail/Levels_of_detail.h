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
      const FT region_growing_minimum_length) {
        
        m_buildings.detect_building_boundaries(
          alpha_shape_size,
          grid_cell_width,
          region_growing_search_size,
          region_growing_noise_level,
          region_growing_angle,
          region_growing_minimum_length);
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

      \param output iterator with polygon vertices given as 3D points.
    */
    template<typename OutputIterator>
    void output_ground_as_polygon(OutputIterator output) const {
      m_ground.return_as_polygon(output);
    }

    /*!
      \brief Returns points used for detecting building boundaries.

      All points are 3D points located on the estimated ground
      plane (see `ground_plane()`).

      \warning `detect_building_boundaries()` should be called
      before calling this method.
        
      \tparam OutputIterator is a model of `OutputIterator`
      that holds `CGAL::Point_3` objects.

      \param output iterator with 3D points.
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
      that holds `std::pair<Point_3, int>` objects.

      \param output iterator with points and assigned to them ids of
      the detected building walls.
    */
    template<typename OutputIterator>
    void output_points_along_building_walls(OutputIterator output) const {
      m_buildings.return_wall_points(output);
    }

    /*!
      \brief Returns polylines that approximate building walls.

      All polylines are 3D segments located on the estimated ground
      plane (see `ground_plane()`).

      \warning `detect_building_boundaries()` should be called
      before calling this method.
        
      \tparam OutputIterator is a model of `OutputIterator`
      that holds `CGAL::Segment_3` objects.

      \param output iterator with 3D segments.
    */
    template<typename OutputIterator>
    void output_building_boundaries_as_polylines(OutputIterator output) const {
      m_buildings.return_boundary_edges(output);
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
