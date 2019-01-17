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
#include <CGAL/Levels_of_detail/internal/Utilities.h>
#include <CGAL/Levels_of_detail/internal/Data_structure.h>

// LOD components.
#include <CGAL/Levels_of_detail/internal/Ground.h>

namespace CGAL {

	namespace Levels_of_detail {

    /*!
      \ingroup PkgLevelsOfDetailRef

      \brief The Levels Of Detail algorithm, constructs levels of detail (LOD) from an input point cloud.

      \tparam GeometricTraits A model of \cgal `Kernel`.

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
      and value type is `GeometricTraits::FT`.

      \tparam Verbose Use if you want to print extra information about execution of the algorithm.
    */
    template<typename GeometricTraits,
             typename InputRange,
             typename PointMap,
             typename SemanticMap,
             typename VisibilityMap = Visibility_from_semantic_map<SemanticMap>,
             typename Verbose = CGAL::Tag_false>
		class Levels_of_detail {

		public:

      /// \name Types
      /// @{
      
      using Traits = GeometricTraits;
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

      /// \endcond

      /// \name Constructor
      /// @{

      /*!
        \brief Initializes data structures for computing levels of detail, 
        given an input range with 3D points, a point, semantic, and visibility map.
      */
      Levels_of_detail(
        const Input_range &input_range,
        Point_map point_map,
        Semantic_map semantic_map,
        Visibility_map visibility_map = VisibilityMap()) : 
      m_data_structure(
        input_range, 
        point_map, 
        semantic_map, 
        visibility_map,
        Verbose::value ? true : false),
      m_ground(m_data_structure) {

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
          `BUILDING_INTERIOR` and extracts the boundary points of this
          alpha shape;

        - uses the union of these boundary points with the points labeled as
          `BUILDING_BOUNDARY` (if any);

        - downsamples this union of points using a regular grid;

        - detects line segments using the region growing approach;

        \warning `compute_planar_ground()` should be called before
        calling this method.
      */
      void detect_building_boundaries(
        FT alpha_shape_size,
        FT grid_cell_width,
        FT region_growing_scale,
        FT region_growing_noise_level,
        FT region_growing_normal_threshold,
        FT region_growing_minimum_length) {
        
        

      }

      /// @}

      /// \name Output
      /// @{

      /*!
        \brief Returns an estimated planar ground.

        \tparam VerticesOutputIterator An output iterator.
        
        \return a planar polygon.
      */
      template<typename VerticesOutputIterator>
      void return_ground_as_polygon(VerticesOutputIterator vertices) const {
        m_ground.return_as_polygon(vertices);
      }

      /// @}

    private:
      Data_structure m_data_structure;
      Ground m_ground;

    }; // end of class

	} // Levels_of_detail

} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_H
