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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_BUILDER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_BUILDER_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <cmath>
#include <vector>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/number_utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Building_builder {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Triangle_3 = typename Traits::Triangle_3;

    void add_lod0() const {

      create_bottom_z();
      create_edges0();
      create_base0();
    }

    void add_lod1() const {

      create_top_z();
      create_edges1();
      create_base1();
      create_walls1();
      create_roofs1();
    }

    void add_lod2() const {
    
      create_edges2();
      create_base2();
      create_walls2();
      create_roofs2();
    }

  private:
    void create_bottom_z() const {
    
    }

    void create_top_z() const {

    }

    void create_top_z_avg() const {

    }

    void create_top_z_max() const {

    }

    void create_edges0() const {

    }

    void create_edges1() const {
      
    }

    void create_edges2() const {
      
    }

    void create_base0() const {
      
    }

    void create_base1() const {
      
    }

    void create_base2() const {
      
    }

    void create_walls1() const {

    }

    void create_walls2() const {

    }

    void create_roofs1() const {
    
    }

    void create_roofs2() const {

    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_BUILDER_H
