// Copyright (c) 2008  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
// Author(s)     : Ophir Setter       <ophirset@post.tau.ac.il>

#ifndef CGAL_ENVElOPE_VORONOI_2_H
#define CGAL_ENVElOPE_VORONOI_2_H

#include <CGAL/Envelope_voronoi_2/Voronoi_diagram_2.h>
#include <CGAL/envelope_3.h>

namespace CGAL {

template <typename InputIterator, typename GeomTraits_, typename TopTraits_>
  void voronoi_2(InputIterator begin, InputIterator end,
                 Envelope_voronoi_2::Voronoi_diagram_2<GeomTraits_,
                 TopTraits_>& min_diagram)
{
  using VD = Envelope_voronoi_2::Voronoi_diagram_2<GeomTraits_, TopTraits_>;
  using Base = typename VD::Base;
  lower_envelope_3 (begin, end, *static_cast<Base*>(&min_diagram));
}

template <typename InputIterator, typename GeomTraits_, typename TopTraits_>
  void farthest_voronoi_2(InputIterator begin, InputIterator end,
                          Envelope_voronoi_2::Voronoi_diagram_2<GeomTraits_,
                          TopTraits_>& max_diagram)
{
  using VD = Envelope_voronoi_2::Voronoi_diagram_2<GeomTraits_, TopTraits_>;
  using Base = typename VD::Base;
  upper_envelope_3 (begin, end, *static_cast<Base*>(&max_diagram));
}

} //namespace CGAL

#endif
