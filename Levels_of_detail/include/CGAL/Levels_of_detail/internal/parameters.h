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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_PARAMETERS_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_PARAMETERS_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <string>

namespace CGAL {
namespace Levels_of_detail {  
namespace internal {

  template<typename FT>
  struct Building_parameters {

    // Constructor.
    Building_parameters(const FT scale_, const FT noise_level_) :
    // Clustering.
    cluster_scale(scale_),
    // Filtering.
    alpha_shape_size_2(scale_ / FT(2)),
    grid_cell_width_2(scale_ / FT(4)),
    // Region growing 2.
    region_growing_scale_2(FT(12)),
    region_growing_noise_level_2(noise_level_),
    region_growing_angle_2(FT(25)),
    region_growing_min_length_2(scale_),
    // Kinetic partitioning 2.
    kinetic_min_face_width_2(scale_ / FT(2)),
    kinetic_max_intersections_2(2),
    // Tagging buildings.
    min_faces_per_footprint(2),
    // Graph cut 2.
    graph_cut_beta_2(FT(1) / FT(10)),
    // Region growing 3.
    region_growing_scale_3(region_growing_scale_2),
    region_growing_noise_level_3(region_growing_noise_level_2),
    region_growing_angle_3(region_growing_angle_2),
    region_growing_min_area_3(scale_),
    // Cleaning roofs.
    min_roof_scale(scale_ / FT(2)),
    // Kinetic partitioning 3.
    kinetic_max_intersections_3(kinetic_max_intersections_2),
    // Graph cut 3.
    graph_cut_beta_3(graph_cut_beta_2)
    { }


    // STEP1: Detecting building boundaries.

    // Clustering.
    FT cluster_scale; // meters

    // Filtering.
    FT alpha_shape_size_2; // meters
    FT grid_cell_width_2; // meters

    // Region growing 2.
    FT region_growing_scale_2; // meters / number of points
    FT region_growing_noise_level_2; // meters
    FT region_growing_angle_2; // degrees
    FT region_growing_min_length_2; // meters


    // STEP2: Computing building footprints.

    // Kinetic partitioning 2.
    FT kinetic_min_face_width_2; // meters
    std::size_t kinetic_max_intersections_2; // number

    // Tagging buildings.
    std::size_t min_faces_per_footprint; // number

    // Graph cut 2.
    FT graph_cut_beta_2; // floating in [0, 1]


    // STEP3: Detecting building roofs.

    // Region growing 3.
    FT region_growing_scale_3; // meters / number of points
    FT region_growing_noise_level_3; // meters
    FT region_growing_angle_3; // degrees
    FT region_growing_min_area_3; // meters

    // Cleaning roofs.
    FT min_roof_scale; // meters


    // STEP4: Computing building roofs.

    // Kinetic partitioning 3.
    std::size_t kinetic_max_intersections_3; // number

    // Graph cut 3.
    FT graph_cut_beta_3; // floating in [0, 1]


    // Update all parameters, which depend on scale and noise_level.
    void update_dependent(const FT scale_, const FT noise_level_) {
      cluster_scale = scale_;
      alpha_shape_size_2 = scale_ / FT(2);
      grid_cell_width_2 = scale_ / FT(4);
      region_growing_noise_level_2 = noise_level_;
      region_growing_min_length_2 = scale_;
      kinetic_min_face_width_2 = scale_ / FT(2);
      region_growing_noise_level_3 = noise_level_;
      region_growing_min_area_3 = scale_;
      min_roof_scale = scale_ / FT(2);
    }
  };

  template<typename FT>
  struct Tree_parameters {

    // Constructor.
    Tree_parameters(const FT scale_, const FT noise_level_) :
    // Clustering.
    cluster_scale(scale_),
    grid_cell_width_2(scale_),
    min_height(noise_level_ * FT(3) / FT(2)),
    // Estimation.
    min_radius_2(noise_level_),
    // Creation.
    min_faces_per_footprint(12),
    // Fitting tree models.
    precision(scale_)
    { }
    

    // STEP1: Computing tree footprints.

    // Clustering.
    FT cluster_scale; // meters
    FT grid_cell_width_2; // meters
    FT min_height; // meters

    // Estimation.
    FT min_radius_2; // meters

    // Creation.
    std::size_t min_faces_per_footprint; // number


    // STEP2: Computing trees.

    // Fitting tree models.
    FT precision; // meters


    // Update all parameters, which depend on scale and noise_level.
    void update_dependent(const FT scale_, const FT noise_level_) {
      cluster_scale = scale_;
      grid_cell_width_2 = scale_;
      min_height = noise_level_ * FT(3) / FT(2);
      min_radius_2 = noise_level_;
      precision = scale_;
    }
  };

  template<typename FT>
  struct Ground_parameters {

    // Constructor.
    Ground_parameters(const FT scale_) :
    // Estimating smooth ground.
    precision(scale_)
    { }

    // Estimating smooth ground.
    FT precision; // meters

    // Update all parameters, which depend on scale and noise_level.
    void update_dependent(const FT scale_) {
      precision = scale_;
    }
  };

  template<typename FT>
  struct Parameters {

    // Path to the input data file.
    std::string data;

    // Label indices defined in the ply header: 
    // ground (gi), 
    // building boundary (bi), 
    // building interior (ii), 
    // vegetation (vi).
    std::string gi, bi, ii, vi;

    // Main parameters.
    FT scale; // meters
    FT noise_level; // meters

    std::size_t min_cluster_size; // fixed number

    // Object parameters.
    Building_parameters<FT> buildings;
    Tree_parameters<FT> trees;
    Ground_parameters<FT> ground;

    // Extrusion type.
    std::size_t extrusion_type; // see enum.h

    // Reconstruction type.
    std::size_t reconstruction_type; // see enum.h

    // Constructor.
    Parameters() : 
    data(""),
    gi("0"), bi("1"), ii("2"), vi("3"),
    scale(FT(4)),
    noise_level(FT(2)),
    min_cluster_size(10),
    buildings(scale, noise_level),
    trees(scale, noise_level),
    ground(scale),
    extrusion_type(2), // max
    reconstruction_type(10) // LOD2
    { }

    // Update all parameters, which depend on scale and noise_level.
    void update_dependent() {
      buildings.update_dependent(scale, noise_level);
      trees.update_dependent(scale, noise_level);
      ground.update_dependent(scale);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_PARAMETERS_H
