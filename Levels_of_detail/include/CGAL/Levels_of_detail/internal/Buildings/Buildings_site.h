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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDINGS_SITE_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDINGS_SITE_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <memory>
#include <vector>
#include <utility>

// Boost includes.
#include <boost/optional/optional.hpp>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_builder.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Buildings_site {

  public:
    using Data_structure = DataStructure;

    using Traits = typename Data_structure::Traits;
    using Points = std::vector<std::size_t>;
    
    using Building = internal::Building<Traits>;
    using Building_ptr = std::shared_ptr<Building>;
    using Building_builder = internal::Building_builder<Traits>;

    Buildings_site(
      const Data_structure& data,
      const Points& points,
      const std::size_t site_index) : 
    m_data(data),
    m_points(points),
    m_site_index(site_index) { 
      CGAL_precondition(m_points.size() > 0);
    }

  private:
    const Data_structure& m_data;
    const Points& m_points;
    const std::size_t m_site_index;

    std::vector<Building> m_buildings;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDINGS_SITE_H
