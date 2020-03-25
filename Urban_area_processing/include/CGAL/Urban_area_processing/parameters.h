// Copyright (c) 2020 SARL GeometryFactory (France).
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

#ifndef CGAL_URBAN_AREA_PROCESSING_PARAMETERS_H
#define CGAL_URBAN_AREA_PROCESSING_PARAMETERS_H

// #include <CGAL/license/Urban_area_processing.h>

#include <CGAL/Urban_area_processing/enum.h>

namespace CGAL {
namespace Urban_area_processing {

  /*!
    \ingroup PkgUrbanAreaProcessingRefGround

    \brief a struct with common parameters, which are used by the classes for 
    reconstructing ground.

    \tparam FT 
    must be a model of `FieldNumberType`.
  */
  template<typename FT>
  struct Ground_parameters {

    /// \name Initialization
    /// @{

    /*!
      \brief sets all parameters from the input `user_scale`.

      \param user_scale
      a user-defined scale

      \pre `user_scale > 0`
    */
    Ground_parameters(
      const FT user_scale) :
    precision(user_scale),
    make_convex_hull(false)
    { }

    /// @}

    /// \name Parameters
    /// @{

    /// Precision of the ground reconstruction in meters.
    FT precision;

    /// If true a convex hull instead of a bouding box is returned. %Default is false.
    bool make_convex_hull;

    /// @}
  };

  /*!
    \ingroup PkgUrbanAreaProcessingRefTrees

    \brief a struct with common parameters, which are used by the classes for 
    reconstructing trees.

    \tparam FT 
    must be a model of `FieldNumberType`.
  */
  template<typename FT>
  struct Tree_parameters {

    /// \name Initialization
    /// @{

    /*!
      \brief sets all parameters from the input `user_scale` and `user_noise`.

      \param user_scale
      a user-defined scale

      \param user_noise
      a user-defined noise level

      \pre `user_scale > 0`
      \pre `user_noise > 0`
    */
    Tree_parameters(
      const FT user_scale, const FT user_noise) :
    scale(user_scale),
    noise(user_noise)
    { }

    /// @}

    /// \name Parameters
    /// @{

    /// Scale in meters.
    FT scale;

    /// Noise level in meters.
    FT noise;

    /// Extrusion type. %Default is maximum height.
    Extrusion_type extrusion_type;

    /// @}
  };

  /*!
    \ingroup PkgUrbanAreaProcessingRefBuildings

    \brief a struct with common parameters, which are used by the classes for 
    reconstructing buildings.

    \tparam FT 
    must be a model of `FieldNumberType`.
  */
  template<typename FT>
  struct Building_parameters {

    /// \name Initialization
    /// @{

    /*!
      \brief sets all parameters from the input `user_scale` and `user_noise`.

      \param user_scale
      a user-defined scale

      \param user_noise
      a user-defined noise level

      \pre `user_scale > 0`
      \pre `user_noise > 0`
    */
    Building_parameters(
      const FT user_scale, 
      const FT user_noise) :
    scale(user_scale),
    noise(user_noise)
    { }

    /// @}

    /// \name Parameters
    /// @{

    /// Scale in meters.
    FT scale;

    /// Noise level in meters.
    FT noise;

    /// Extrusion type. %Default is maximum height.
    Extrusion_type extrusion_type;

    /// @}
  };

  /// \cond SKIP_IN_MANUAL

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
    FT noise; // meters

    // Object parameters.
    Building_parameters<FT> buildings;
    Tree_parameters<FT> trees;
    Ground_parameters<FT> ground;

    // Constructor.
    Parameters() : 
    data(""),
    gi("0"), bi("1"), ii("2"), vi("3"),
    scale(FT(4)),
    noise(FT(2)),
    buildings(scale, noise),
    trees(scale, noise),
    ground(scale)
    { }
  };

  /// \endcond

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_PARAMETERS_H
