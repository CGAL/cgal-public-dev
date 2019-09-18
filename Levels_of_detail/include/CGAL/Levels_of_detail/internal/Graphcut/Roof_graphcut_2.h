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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_GRAPHCUT_2_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_GRAPHCUT_2_H

// CGAL includes.
#include <CGAL/assertions.h>

#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Roof_graphcut_2 {

  public:
    using Traits = GeomTraits;
    
		using Partition_2 = internal::Partition_2<Traits>;

    using FT = typename Traits::FT;

    using Face = typename Partition_2::Face;
    using Edge = typename Partition_2::Edge;

		using Size_pair = std::pair<std::size_t, std::size_t>;
		using Alpha_expansion = CGAL::internal::Alpha_expansion_graph_cut_boost;

    Roof_graphcut_2(const FT graphcut_beta) : 
    m_beta(graphcut_beta)
    { }

    void apply(Partition_2& partition) const {
      
      if (partition.empty()) return;
      auto& pfaces = partition.faces;
      auto& pedges = partition.edges;

    }

  private:
    const FT m_beta;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_GRAPHCUT_2_H
