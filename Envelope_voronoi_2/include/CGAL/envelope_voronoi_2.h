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
//
// Author(s)     : Ophir Setter       <ophirset@post.tau.ac.il>

#ifndef CGAL_ENVElOPE_VORONOI_2_H
#define CGAL_ENVElOPE_VORONOI_2_H

#include <CGAL/Envelope_voronoi_2/Voronoi_diagram_2.h>
#include <CGAL/envelope_3.h>


namespace CGAL {

template <class InputIterator, class GeomTraits_, class TopTraits_>
  void voronoi_2 (InputIterator begin, InputIterator end,
                  Envelope_voronoi_2::Voronoi_diagram_2<GeomTraits_, 
                  TopTraits_>& min_diagram)
{
  typedef Envelope_voronoi_2::Voronoi_diagram_2<GeomTraits_, TopTraits_>
    VD;
  typedef typename VD::Base                                   Base;
  lower_envelope_3 (begin, end, *static_cast<Base*>(&min_diagram));
}

template <class InputIterator, class GeomTraits_, class TopTraits_>
  void farthest_voronoi_2 (InputIterator begin, InputIterator end,
                           Envelope_voronoi_2::Voronoi_diagram_2<GeomTraits_, 
                           TopTraits_>& max_diagram)
{
  typedef Envelope_voronoi_2::Voronoi_diagram_2<GeomTraits_, TopTraits_>
    VD;
  typedef typename VD::Base                                   Base;
  upper_envelope_3 (begin, end, *static_cast<Base*>(&max_diagram));
}


} //namespace CGAL

#endif
