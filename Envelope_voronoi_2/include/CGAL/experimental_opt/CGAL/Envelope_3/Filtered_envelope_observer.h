// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     :  Ophir Setter           <ophirset@post.tau.ac.il>

/*! \file Filtered_envelope_observer.h
*/

#ifndef CGAL_FILTERED_ENVELOPE_OBSERVER_H
#define CGAL_FILTERED_ENVELOPE_OBSERVER_H

#include <CGAL/Arr_observer.h>

#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>


namespace CGAL {

template <class Arrangement_2_>
class Filtered_envelope_observer : Arr_observer<Arrangement_2_>
{
 public:
  typedef Arrangement_2_                               Arrangement_2;
  typedef Arr_observer<Arrangement_2>                  Base;
  typedef Filtered_envelope_observer<Arrangement_2>    Self;

  typedef typename Base::Vertex_handle                 Vertex_handle;
  typedef typename Base::Halfedge_handle               Halfedge_handle;
  typedef typename Base::X_monotone_curve_2            X_monotone_curve_2;

  typedef typename Arrangement_2::Halfedge_around_vertex_circulator
    Halfedge_around_vertex_circulator;
  typedef typename Arrangement_2::Vertex               Vertex;
  typedef typename Arrangement_2::Halfedge             Halfedge;
  typedef typename Vertex::Surfaces_container          Surfaces_container;
  typedef typename Surfaces_container::value_type      Surface;

  Filtered_envelope_observer(Arrangement_2 &arr, Surface s1, Surface s2) 
    : Arr_observer<Arrangement_2_>(arr), _s1(s1), _s2(s2)
  {}
  
  void print_vertex(Vertex_handle v)
  {
    typedef typename Vertex::Surfaces_container       Surfaces_container;

    if (v->is_at_infinity())
    {
      std::cerr << "surfaces of fictitious vertex are:" << std::endl;
    }
    else
    {
      std::cerr << "surfaces of " << v->point() << " are:" << std::endl;
    }
    typename Surfaces_container::iterator it;
    for (it = v->surfaces().begin(); it != v->surfaces().end(); ++it)
    {
      std::cerr << *it << " ";
    }
  }
  
  virtual void after_create_vertex (Vertex_handle v)
  {
    v->surfaces().insert(_s1);
    v->surfaces().insert(_s2);
  }

  // The function adds surfaces to the given face. First it added the
  // current surfaces (probably the edges of e) and also the surfaces
  // of the neighboring edges of e - as the neighboring edges of e are
  // bisectors of of either _s1 with another surface or _s2 with another
  // surface.
  void add_surfaces_to_vertex(Vertex_handle v, Halfedge_handle e)
  {
    CGAL_envelope_voronoi_assertion(v == e->target());

    v->surfaces().insert(_s1);
    v->surfaces().insert(_s2);

    // add also the surfaces of the neighboring edges of e
    const Surfaces_container &sur1 = e->next()->surfaces();
    v->surfaces().insert(sur1.begin(), sur1.end());

    const Surfaces_container &sur2 = e->twin()->prev()->surfaces();
    v->surfaces().insert(sur2.begin(), sur2.end()); 
  }

  virtual void after_create_edge (Halfedge_handle e)
  {
/*     std::cerr << "before creating edge surfaces are: " << std::endl; */
/*     print_vertex(e->source()); */
/*     std::cerr << std::endl; */
/*     print_vertex(e->target()); */
/*     std::cerr << std::endl; */

    e->surfaces().insert(_s1);
    e->surfaces().insert(_s2);
    e->twin()->surfaces().insert(_s1);
    e->twin()->surfaces().insert(_s2);

    // add surfaces to the vetrices of the edge
    add_surfaces_to_vertex(e->target(), e);
    add_surfaces_to_vertex(e->source(), e->twin());

/*     std::cerr << "after creating edge surfaces are: " << std::endl; */
/*     print_vertex(e->source()); */
/*     std::cerr << std::endl; */
/*     print_vertex(e->target()); */
/*     std::cerr << std::endl; */
  }

  virtual void after_split_edge (Halfedge_handle e1, Halfedge_handle e2)
  {
    // The code assumes that e1 is the original halfedge that is split and
    // e2 is the new halfedge that was created from the split.

    const Surfaces_container &sur = e1->surfaces();
    e2->surfaces().insert(sur.begin(), sur.end());
    e2->twin()->surfaces().insert(sur.begin(), sur.end());

    // todo: check if this is really needed
    Vertex_handle v = e1->target();
    v->surfaces().insert(sur.begin(), sur.end());
  }

 protected:
  Surface _s1;
  Surface _s2;
};

} //namespace CGAL

#endif
