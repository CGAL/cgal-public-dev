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
// $URL:
// $Id:
//
// Author(s)     : Ophir Setter <ophir.setter@cs.tau.ac.il

#ifndef CGAL_MWA_OVERLAY_TRAITS_2_H
#define CGAL_MWA_OVERLAY_TRAITS_2_H

/*! \file
 * Definition of overlay-traits for computing minimum width annulus using
 * Vornoi diagrmas.
 */

namespace CGAL {

#define MWA_OVERLAY_TRAITS_FUNC(func_name, type1, type2, type3)         \
  virtual void func_name (type1 v1, type2 v2, type3 v) const            \
  {                                                                     \
    v->first.add_data(v1->begin_data(), v1->end_data());                \
    v->second.add_data(v2->begin_data(), v2->end_data());               \
  }

/*! \class Min_width_annulus_overlay_2
 * Class used when computing the minimum width annulus. For each new created
 * feature, we save the info from both arrangements, as well as from which
 * arrangement it came from.
 */
template <typename ArrangementA, class ArrangementB, class ArrangementR>
class MWA_overlay_traits_2 {
public:
  using Vertex_handle_A = typename ArrangementA::Vertex_const_handle;
  using Halfedge_handle_A = typename ArrangementA::Halfedge_const_handle;
  using Face_handle_A = typename ArrangementA::Face_const_handle;

  using Vertex_handle_B = typename ArrangementB::Vertex_const_handle;
  using Halfedge_handle_B = typename ArrangementB::Halfedge_const_handle;
  using Face_handle_B = typename ArrangementB::Face_const_handle;

  using Vertex_handle_R = typename ArrangementR::Vertex_handle;
  using Halfedge_handle_R = typename ArrangementR::Halfedge_handle;
  using Face_handle_R = typename ArrangementR::Face_handle;

  /*! Destructor. */
  virtual ~MWA_overlay_traits_2() {}

  MWA_OVERLAY_TRAITS_FUNC(create_vertex, Vertex_handle_A, Vertex_handle_B, \
                          Vertex_handle_R);
  MWA_OVERLAY_TRAITS_FUNC(create_vertex, Vertex_handle_A, Halfedge_handle_B, \
                          Vertex_handle_R);
  MWA_OVERLAY_TRAITS_FUNC(create_vertex, Halfedge_handle_A, Vertex_handle_B, \
                          Vertex_handle_R);
  MWA_OVERLAY_TRAITS_FUNC(create_vertex, Vertex_handle_A, Face_handle_B, \
                          Vertex_handle_R);
  MWA_OVERLAY_TRAITS_FUNC(create_vertex, Face_handle_A, Vertex_handle_B, \
                          Vertex_handle_R);
  MWA_OVERLAY_TRAITS_FUNC(create_vertex, Halfedge_handle_A, Halfedge_handle_B, \
                          Vertex_handle_R);

  MWA_OVERLAY_TRAITS_FUNC(create_edge, Halfedge_handle_A, Halfedge_handle_B, \
                          Halfedge_handle_R);
  MWA_OVERLAY_TRAITS_FUNC(create_edge, Halfedge_handle_A, Face_handle_B, \
                          Halfedge_handle_R);
  MWA_OVERLAY_TRAITS_FUNC(create_edge, Face_handle_A, Halfedge_handle_B, \
                          Halfedge_handle_R);

  MWA_OVERLAY_TRAITS_FUNC(create_face, Face_handle_A, Face_handle_B,    \
                          Face_handle_R);

};

} //namespace CGAL

#endif // CGAL_MWA_OVERLAY_TRAITS_2_H
