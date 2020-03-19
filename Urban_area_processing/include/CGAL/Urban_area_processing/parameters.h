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

namespace CGAL {
namespace Urban_area_processing {

  template<typename FT>
  struct Building_parameters {

    // Constructor.
    Building_parameters(
      const FT scale_, const FT noise_) :
    scale(scale_),
    noise(noise_)
    { }

    FT scale; // meters
    FT noise; // meters
  };

  template<typename FT>
  struct Tree_parameters {

    // Constructor.
    Tree_parameters(
      const FT scale_, const FT noise_) :
    scale(scale_),
    noise(noise_)
    { }

    FT scale; // meters
    FT noise; // meters
  };

  template<typename FT>
  struct Ground_parameters {

    // Constructor.
    Ground_parameters(const FT scale_) :
    precision(scale_)
    { }

    FT precision; // meters
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

} // Urban_area_processing
} // CGAL

#endif // CGAL_URBAN_AREA_PROCESSING_PARAMETERS_H
