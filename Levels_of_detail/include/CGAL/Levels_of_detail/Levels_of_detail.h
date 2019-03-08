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

#include <CGAL/Levels_of_detail/internal/Reconstruction/LOD0.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/LOD1.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/LOD2.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_3.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_3_fuzzy_sphere_connectivity.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_3_empty_conditions.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Roof_cleaner.h>

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
    using Point_3 = typename Traits::Point_3;

    using Data_structure = internal::Data_structure<
    Traits, 
    Input_range, 
    Point_map, 
    Semantic_map, 
    Visibility_map>;

    using Ground = internal::Ground<Data_structure>;
    using Buildings = internal::Buildings<Data_structure>;
    using Vegetation = internal::Vegetation<Data_structure>;

    using LOD0 = internal::LOD0<Data_structure>;
    using LOD1 = internal::LOD1<Data_structure>;
    using LOD2 = internal::LOD2<Data_structure>;

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
    m_semantic_map(semantic_map),
    m_visibility_map(visibility_map), 
    m_data_structure(
      input_range, 
      point_map, 
      semantic_map, 
      visibility_map,
      Verbose::value ? true : false),
    m_ground(m_data_structure),
    m_buildings(m_data_structure),
    m_vegetation(m_data_structure),
    m_is_components_based(false),
    m_save(false) { 

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

    template<typename Saver>
    void initialize(
      const FT scale, 
      Saver& saver, 
      const std::string path,
      const std::string gi,
      const std::string bi,
      const std::string ii,
      const std::string vi) {

      m_components.clear();
      m_structures.clear();
      m_com_buildings.clear();

      using Label_map = typename Input_range:: template Property_map<int>;
      const auto building_points = m_data_structure.building_interior_points();
      
      Label_map b_labels = 
      m_data_structure.input_range. template property_map<int>("label").first;
      
      bool success;
      Label_map p_labels;

      Input_range point_range;
      boost::tie(p_labels, success) = point_range. template add_property_map<int>("label", -1);
      
      for (const auto& item : building_points) {
        auto it = point_range.insert(get(m_data_structure.point_map, item));
        p_labels[*it] = b_labels[item];
      }

      using Points_connectivity_3 = 
      internal::Points_3_fuzzy_sphere_connectivity<Traits, Input_range, Point_map>;
      using Normals_estimator_3 = 
      internal::Estimate_normals_3<Traits, Input_range, Point_map, Points_connectivity_3>;
      using Points_conditions_3 = 
      internal::Points_3_empty_conditions<Traits, Input_range, Point_map>;
      using Points_region_growing_3 = 
      internal::Region_growing<Points_connectivity_3, Points_conditions_3>;
      using Roof_cleaner = 
      internal::Roof_cleaner<Traits, Input_range, Point_map>;

      Points_connectivity_3 connectivity(
        point_range,
        point_range.point_map(), 
        scale);

      Normals_estimator_3 estimator(
        point_range, 
        point_range.point_map(),
        connectivity);

      Points_conditions_3 conditions(
        point_range,
        point_range.point_map());

      std::vector<std::size_t> indices(point_range.size());
      for (std::size_t i = 0; i < point_range.size(); ++i)
        indices[i] = i;
      
      Points_region_growing_3 region_growing(
        indices,
        connectivity,
        conditions);

      std::vector< std::vector<std::size_t> > regions;
      region_growing.detect(regions);
        
      const Roof_cleaner cleaner(
        point_range, 
        point_range.point_map(),
        estimator.normals(),
        scale);

      cleaner.clean(regions);
      for (const auto& region : regions) {

        Label_map c_labels;
        Input_range component;
        boost::tie(c_labels, success) = component. template add_property_map<int>("label", -1);

        for (const std::size_t idx : region) {
          auto it = component.insert(
            get(point_range.point_map(), *(point_range.begin() + idx)));
          c_labels[*it] = p_labels[*(point_range.begin() + idx)];
        }

        FT minx = internal::max_value<FT>(); FT miny = internal::max_value<FT>();
        FT maxx = -internal::max_value<FT>(); FT maxy = -internal::max_value<FT>();

        for (const std::size_t idx : region) {
          const Point_3& p = 
          get(point_range.point_map(), *(point_range.begin() + idx));
          
          minx = CGAL::min(minx, p.x()); miny = CGAL::min(miny, p.y());
          maxx = CGAL::max(maxx, p.x()); maxy = CGAL::max(maxy, p.y());
        }

        for (auto pit = m_data_structure.input_range.begin(); pit != m_data_structure.input_range.end(); ++pit) {
          const Point_3& p = get(m_data_structure.point_map, *pit);

          if (p.x() >= minx && p.x() <= maxx && 
              p.y() >= miny && p.y() <= maxy &&  
              (b_labels[*pit] == std::stoi(gi) || 
               b_labels[*pit] == std::stoi(vi))) {
            
            auto it = component.insert(p);
            c_labels[*it] = b_labels[*pit];
          }
        }
        m_components.push_back(component);
      }

      std::cout << std::endl << m_components.size() << " components found" << std::endl;
      for (const auto& component : m_components) {

        Label_map c_labels = 
        component. template property_map<int>("label").first;
        Semantic_map semantic_map(c_labels, gi, bi, ii, vi);
        Visibility_map visibility_map(semantic_map);

        Data_structure structure(
          component, 
          component.point_map(), 
          semantic_map, 
          visibility_map,
          false);
        m_structures.push_back(structure);

        if (m_save)
          saver.export_point_set(m_components[0], path);
      }

      for (auto& structure : m_structures) {
        Buildings buildings(structure);
        m_com_buildings.push_back(buildings);
      }

      CGAL_assertion(m_structures.size() == m_components.size());
      CGAL_assertion(m_com_buildings.size() == m_structures.size());

      m_is_components_based = true;
    }

    template<typename OutputIterator>
    void output_components(OutputIterator output) const {

      long i = 0;
      for (const auto& component : m_components) {  
        for (const auto& item : component)
          *(output++) = std::make_pair(
            get(component.point_map(), item), i);
        ++i;
      }
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

      if (m_is_components_based) {
        for (auto& structure : m_structures)
          structure.planar_ground = m_data_structure.planar_ground;
      }
    }

    /*!
      \brief Computes a smooth representation of the ground.
      
      The ground is represented as Delaunay triangulation with
      associated ground heights, which is computed upon the points
      semantically labeled as `Semantic_label::GROUND`.
    */
    void compute_smooth_ground(const FT ground_precision) {
      m_ground.make_smooth(ground_precision);

      if (m_is_components_based) {
        for (auto& structure : m_structures)
          structure.smooth_ground = m_data_structure.smooth_ground;
      }
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
        
        if (m_is_components_based) {
          for (std::size_t i = 0; i < m_com_buildings.size(); ++i) {
            
            m_com_buildings[i].detect_boundaries(
              alpha_shape_size,
              grid_cell_width,
              region_growing_search_size,
              region_growing_noise_level,
              region_growing_angle,
              region_growing_min_length);
          }
          
          if (m_save) {
            m_data_structure.building_boundary_points_2 = m_structures[0].building_boundary_points_2;
            m_data_structure.building_boundary_indices_2 = m_structures[0].building_boundary_indices_2;
            m_data_structure.buildings = m_structures[0].buildings;
            m_data_structure.building_polygon_faces_2 = m_structures[0].building_polygon_faces_2;
            m_data_structure.building_clusters = m_structures[0].building_clusters;
          }
          return;
        }

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
        means certainly outside and 1 means certainly inside;

      - tags subsets of all polygon faces with the visibility value >= 0.5
        that form separate buildings.

      \warning `detect_building_boundaries()` should be called
      before calling this method.
    */
    void compute_building_footprints(
      const FT kinetic_min_face_width,
      const std::size_t kinetic_max_intersections,
      const std::size_t min_faces_per_building) {

      if (m_is_components_based) {
        for (std::size_t i = 0; i < m_com_buildings.size(); ++i) {
          m_com_buildings[i].compute_footprints(
            kinetic_min_face_width,
            kinetic_max_intersections,
            min_faces_per_building);
        }
        
        if (m_save) {
          m_data_structure.building_boundary_points_2 = m_structures[0].building_boundary_points_2;
          m_data_structure.planar_ground.plane = m_structures[0].planar_ground.plane;
          m_data_structure.building_boundary_indices_2 = m_structures[0].building_boundary_indices_2;
          m_data_structure.buildings = m_structures[0].buildings;
          m_data_structure.building_polygon_faces_2 = m_structures[0].building_polygon_faces_2;
          m_data_structure.building_clusters = m_structures[0].building_clusters;
        }
        return;
      }

      m_buildings.compute_footprints(
        kinetic_min_face_width,
        kinetic_max_intersections,
        min_faces_per_building);
    }

    void finilize_lod0() {

      if (!m_is_components_based) 
        return;

      if (m_save) {
        m_data_structure.building_boundary_points_2 = m_structures[0].building_boundary_points_2;
        m_data_structure.planar_ground.plane = m_structures[0].planar_ground.plane;
        m_data_structure.building_boundary_indices_2 = m_structures[0].building_boundary_indices_2;
        m_data_structure.buildings = m_structures[0].buildings;
        m_data_structure.building_polygon_faces_2 = m_structures[0].building_polygon_faces_2;
        m_data_structure.building_clusters = m_structures[0].building_clusters;

        return;
      }

      m_data_structure.buildings.clear();
      for (const auto& structure : m_structures) {
        for (const auto& building : structure.buildings)
          m_data_structure.buildings.push_back(building);
      }
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
      
      if (m_is_components_based) {
        for (std::size_t i = 0; i < m_com_buildings.size(); ++i)
          m_com_buildings[i].extrude_footprints(extrusion_type);

        if (m_save) {
          m_data_structure.building_boundary_points_2 = m_structures[0].building_boundary_points_2;
          m_data_structure.planar_ground.plane = m_structures[0].planar_ground.plane;
          m_data_structure.building_boundary_indices_2 = m_structures[0].building_boundary_indices_2;
          m_data_structure.buildings = m_structures[0].buildings;
          m_data_structure.building_polygon_faces_2 = m_structures[0].building_polygon_faces_2;
          m_data_structure.building_clusters = m_structures[0].building_clusters;
        }
        return;
      }
      
      m_buildings.extrude_footprints(extrusion_type);
    }

    void finilize_lod1() {

      if (!m_is_components_based) 
        return;

      if (m_save) {
        m_data_structure.building_boundary_points_2 = m_structures[0].building_boundary_points_2;
        m_data_structure.planar_ground.plane = m_structures[0].planar_ground.plane;
        m_data_structure.building_boundary_indices_2 = m_structures[0].building_boundary_indices_2;
        m_data_structure.buildings = m_structures[0].buildings;
        m_data_structure.building_polygon_faces_2 = m_structures[0].building_polygon_faces_2;
        m_data_structure.building_clusters = m_structures[0].building_clusters;

        return;
      }

      m_data_structure.buildings.clear();
      for (const auto& structure : m_structures) {
        for (const auto& building : structure.buildings)
          m_data_structure.buildings.push_back(building);
      }
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
      
      if (m_is_components_based) {
        for (std::size_t i = 0; i < m_com_buildings.size(); ++i) {
          m_com_buildings[i].detect_roofs(
            region_growing_search_size,
            region_growing_noise_level,
            region_growing_angle,
            region_growing_min_area,
            min_size);
        }

        if (m_save) {
          m_data_structure.building_boundary_points_2 = m_structures[0].building_boundary_points_2;
          m_data_structure.planar_ground.plane = m_structures[0].planar_ground.plane;
          m_data_structure.building_boundary_indices_2 = m_structures[0].building_boundary_indices_2;
          m_data_structure.buildings = m_structures[0].buildings;
          m_data_structure.building_polygon_faces_2 = m_structures[0].building_polygon_faces_2;
          m_data_structure.building_clusters = m_structures[0].building_clusters;
        }
        return;
      }

      m_buildings.detect_roofs(
        region_growing_search_size,
        region_growing_noise_level,
        region_growing_angle,
        region_growing_min_area,
        min_size);
    }

    /*!
      \brief Computes building roofs.

      This method:

      - creates the partitioning by extending initial building walls, roofs, and 
        ground represented as polygons until the defined number of 
        intersections with other polygons is reached;

      - applies the visibility computation that assignes to each polyhedral 
        facet of the partitioning a visibility value in the range [0, 1], where 0 
        means certainly outside and 1 means certainly inside;

      - corrects the visibility estimations by applying a 3D graphcut;

      - extracts 3D polygons that represent exact roofs for each building.

      \warning `detect_building_roofs()` should be called 
      before calling this method.
    */
    void compute_building_roofs(
      const std::size_t kinetic_max_intersections,
      const FT graph_cut_beta_3) {
      
      if (m_is_components_based) {
        for (std::size_t i = 0; i < m_com_buildings.size(); ++i) {
          m_com_buildings[i].compute_roofs(
            kinetic_max_intersections,
            graph_cut_beta_3);
        }

        if (m_save) {
          m_data_structure.building_boundary_points_2 = m_structures[0].building_boundary_points_2;
          m_data_structure.planar_ground.plane = m_structures[0].planar_ground.plane;
          m_data_structure.building_boundary_indices_2 = m_structures[0].building_boundary_indices_2;
          m_data_structure.buildings = m_structures[0].buildings;
          m_data_structure.building_polygon_faces_2 = m_structures[0].building_polygon_faces_2;
          m_data_structure.building_clusters = m_structures[0].building_clusters;
        }
        return;
      }

      m_buildings.compute_roofs(
        kinetic_max_intersections,
        graph_cut_beta_3);
    }

    void finilize_lod2() {

      if (!m_is_components_based) 
        return;

      if (m_save) {
        m_data_structure.building_boundary_points_2 = m_structures[0].building_boundary_points_2;
        m_data_structure.planar_ground.plane = m_structures[0].planar_ground.plane;
        m_data_structure.building_boundary_indices_2 = m_structures[0].building_boundary_indices_2;
        m_data_structure.buildings = m_structures[0].buildings;
        m_data_structure.building_polygon_faces_2 = m_structures[0].building_polygon_faces_2;
        m_data_structure.building_clusters = m_structures[0].building_clusters;

        return;
      }

      m_data_structure.buildings.clear();
      for (const auto& structure : m_structures) {
        for (const auto& building : structure.buildings)
          m_data_structure.buildings.push_back(building);
      }
    }

    /*!
      \brief Creates 3D tree icons.

      This method:

      - creates tree icons.

      \warning `compute_tree_footprints()` should be called 
      before calling this method.
    */
    void fit_tree_models(const FT precision) {
      m_vegetation.fit_tree_models(precision);
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

    /*!
      \brief Returns input to the partitioning algorithm.

      \warning `compute_building_roofs()` should be called
      before calling this method.

      \tparam VerticesOutputIterator is a model of `OutputIterator`
      that holds `Point_3` objects.

      \tparam FacesOutputIterator is a model of `OutputIterator`
      that holds `std::pair<std::vector<std::size_t>, std::size_t>` 
      objects, where the first item in the pair holds indices of the face 
      vertices and the second item is the building index.

      \param output_vertices an iterator with all vertices of the polygon soup.

      \param output_faces an iterator with all faces of the polygon soup
      given as vectors of indices in `output_vertices` and the corresponding
      building indices.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_partitioning_in_3_as_polygon_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      
      m_buildings.return_partitioning_input_3(output_vertices, output_faces);
    }

    /*!
      \brief Returns output of the partitioning algorithm.

      \warning `compute_building_roofs()` should be called
      before calling this method.

      \tparam VerticesOutputIterator is a model of `OutputIterator`
      that holds `Point_3` objects.

      \tparam FacesOutputIterator is a model of `OutputIterator`
      that holds `std::pair<std::vector<std::size_t>, std::size_t>` 
      objects, where the first item in the pair holds indices of the face 
      vertices and the second item is the building index.

      \param output_vertices an iterator with all vertices of the polygon soup.

      \param output_faces an iterator with all faces of the polygon soup
      given as vectors of indices in `output_vertices` and the corresponding
      building indices.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_partitioning_out_3_as_polygon_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      
      m_buildings.return_partitioning_output_3(
        output_vertices, output_faces, false);
    }

    /*!
      \brief Returns polygons, which bound buildings.

      \warning `compute_building_roofs()` should be called
      before calling this method.

      \tparam VerticesOutputIterator is a model of `OutputIterator`
      that holds `Point_3` objects.

      \tparam FacesOutputIterator is a model of `OutputIterator`
      that holds `std::pair<std::vector<std::size_t>, std::size_t>` 
      objects, where the first item in the pair holds indices of the face 
      vertices and the second item is the building index.

      \param output_vertices an iterator with all vertices of the polygon soup.

      \param output_faces an iterator with all faces of the polygon soup
      given as vectors of indices in `output_vertices` and the corresponding
      building indices.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_building_bounds_as_polygon_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      
      m_buildings.return_partitioning_output_3(
        output_vertices, output_faces, true);
    }

    /*!
      \brief Returns all trees as a triangle soup.
        
      Each triangle is associated to the index of the corresponding
      tree.

      \warning `fit_tree_models()` should be called before
      calling this method.

      \tparam VerticesOutputIterator is a model of `OutputIterator`
      that holds `Point_3` objects.

      \tparam FacesOutputIterator is a model of `OutputIterator`
      that holds `std::pair<cpp11::array<std::size_t, 3>, std::size_t>` objects,
      where the first item in the pair holds indices of the face vertices and
      the second item is the tree index.

      \param output_vertices an iterator with all vertices of the triangle soup.

      \param output_faces an iterator with all faces of the triangle soup
      given as arrays of indices in `output_vertices` and the corresponding
      tree indices.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_trees_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {
      
      m_vegetation.return_trees(output_vertices, output_faces);
    }

    /*!
      \brief Returns LOD.

      \warning `build()` should be called before calling this method.

      \tparam VerticesOutputIterator is a model of `OutputIterator`
      that holds `Point_3` objects.

      \tparam FacesOutputIterator is a model of `OutputIterator`
      that holds `std::pair< cpp11::array<std::size_t, 3>, Urban_object_type>` objects,
      where the first item in the pair holds indices of the face vertices and second
      item is the type of the urban object this face belongs to.

      \param output_vertices an iterator with all vertices of the triangle soup.

      \param output_faces an iterator with all faces of the triangle soup
      given as arrays of indices in `output_vertices` and the corresponding urban 
      object types.

      \param lod_type type of the LOD output. Can be LOD0, LOD1, or LOD2.
    */
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_LOD_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      const Reconstruction_type lod_type,
      const FT ground_precision = FT(0)) {
      
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
      }
    }

    /// @}

  private:
    Semantic_map m_semantic_map;
    Visibility_map m_visibility_map;
    Data_structure m_data_structure;
    Ground m_ground;
    Buildings m_buildings;
    Vegetation m_vegetation;
    
    std::vector<Data_structure> m_structures;
    std::vector<Input_range> m_components;

    bool m_is_components_based;
    std::vector<Buildings> m_com_buildings;
    bool m_save;

  }; // end of class

} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_H
