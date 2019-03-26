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
      \brief reconstructs a model of type `CGAL::Levels_of_detail::Reconstruction_type`.

      \param scale
      scale parameter

      \param noise_level
      noise level parameter

      \param ground_precision
      max distance between input points and a reconstructed ground

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
      m_data.parameters.ground.precision = ground_precision;
      m_data.parameters.reconstruction_type = static_cast<std::size_t>(
        reconstruction_type);
    }

    /// @}

    /// \name Ground
    /// @{

    /*!
      \brief computes a planar representation of the ground.

      The plane is estimated through the principal component analysis on
      the points semantically labeled as `CGAL::Levels_of_detail::Semantic_label::GROUND`.
    */
    void compute_planar_ground() {
      m_ground.make_planar();
    }

    /*!
      \brief computes a smooth representation of the ground.
      
      The ground is represented as a Delaunay triangulation with
      associated ground heights, which is computed upon the points
      semantically labeled as `CGAL::Levels_of_detail::Semantic_label::GROUND`.

      \param ground_precision
      max distance between input points and a reconstructed ground
    */
    void compute_smooth_ground(const FT ground_precision) {
      m_data.parameters.ground.precision = ground_precision;
      m_ground.make_smooth();
    }

    /*!
      \brief computes tree footprints.

      This method:

      - extracts points, which form potential trees, from all points labeled as
      `CGAL::Levels_of_detail::Semantic_label::VEGETATION`;

      - estimates tree center point, radius, and height;

      - creates tree footprints.

      \param cluster_scale
      cluster scale parameter

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
      const FT cluster_scale,
      const FT grid_cell_width_2, 
      const FT min_height, 
      const FT min_radius_2,
      const std::size_t min_faces_per_footprint) {
      m_data.parameters.trees.cluster_scale = cluster_scale;
      m_data.parameters.trees.grid_cell_width_2 = grid_cell_width_2;
      m_data.parameters.trees.min_height = min_height;
      m_data.parameters.trees.min_radius_2 = min_radius_2;
      m_data.parameters.trees.min_faces_per_footprint = min_faces_per_footprint;
      m_trees.compute_footprints();
    }

    /// @}

    /// \name Output
    /// @{

    /*!
      \brief returns ground as triangle soup.

      \warning `compute_planar_ground()` or `compute_smooth_ground()` should be 
      called before calling this method

      \tparam VerticesOutputIterator 
      must be an output iterator whose value type is `Kernel::Point_3`.

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair<std::vector<std::size_t>, CGAL::Levels_of_detail::Urban_object_type::GROUND>`.

      \pre `ground_type == Reconstruction_type::PLANAR_GROUND ||
            ground_type == Reconstruction_type::SMOOTH_GROUND`

      \return a pair of the past-the-end iterators
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_ground_as_triangle_soup(
      const Reconstruction_type ground_type,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {

      CGAL_precondition(
        ground_type == Reconstruction_type::PLANAR_GROUND ||
        ground_type == Reconstruction_type::SMOOTH_GROUND);
      return m_ground.output(vertices, faces, ground_type);
    }

    /*!
      \brief returns LOD as triangle soup.

      \warning `build()` 
      should be called before calling this method

      \tparam VerticesOutputIterator
      must be an output iterator whose value type is `Kernel::Point_3`. 

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
      `std::pair< std::vector<std::size_t>, CGAL::Levels_of_detail::Urban_object_type>`.

      \pre `lod_type == Reconstruction_type::LOD0 ||
            lod_type == Reconstruction_type::LOD1 ||
            lod_type == Reconstruction_type::LOD2`

      \return a pair of the past-the-end iterators
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> >
    output_LOD_as_triangle_soup(
      const Reconstruction_type lod_type,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      
      CGAL_precondition(
        lod_type == Reconstruction_type::LOD0 ||
        lod_type == Reconstruction_type::LOD1 ||
        lod_type == Reconstruction_type::LOD2);
      return m_lods.output(vertices, faces, lod_type);
    }

    /*!
      \brief returns a point set.
        
      Returns data related to the intermediate steps of the algorithm, which
      can be stored as a point set.

      \tparam OutputIterator 
      must be an output iterator whose value type is 
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> output_points(
      const Intermediate_step step,
      OutputIterator output) const {
      
      switch (step) {
        case Intermediate_step::TREE_CLUSTERS: {
          return m_trees.get_tree_clusters(output);
        }
        case Intermediate_step::TREE_POINTS: {
          return m_trees.get_tree_points(output);
        }
        default: return boost::none;
      }
    }

    /*!
      \brief returns polylines.
        
      Returns data related to the intermediate steps of the algorithm, which
      can be stored as a set of polylines.

      \tparam OutputIterator 
      must be an output iterator whose value type is 
    */
    template<typename OutputIterator>
    boost::optional<OutputIterator> output_polylines(
      const Intermediate_step step,
      OutputIterator output) const {
      
      switch (step) {
        case Intermediate_step::TREE_BOUNDARIES: {
          return m_trees.get_tree_boundaries(output);
        }
        default: return boost::none;
      }
    }

    /*!
      \brief returns mesh.
        
      Returns data related to the intermediate steps of the algorithm, which
      can be stored as a mesh.

      \tparam VerticesOutputIterator
      must be an output iterator whose value type is `Kernel::Point_3`. 

      \tparam FacesOutputIterator 
      must be an output iterator whose value type is 
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> > 
    output_mesh(
      const Intermediate_step step,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces) const {
      
      switch (step) {
        case Intermediate_step::TREE_FOOTPRINTS: {
          return m_trees.get_tree_footprints(vertices, faces);
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
    Ground m_ground;
    Buildings m_buildings;
    Trees m_trees;
    LODs m_lods;

  }; // Levels_of_detail

} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_H
