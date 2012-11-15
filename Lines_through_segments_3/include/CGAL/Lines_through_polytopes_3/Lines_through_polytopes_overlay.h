// Copyright (c) 2010  Tel-Aviv University (Israel).
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
// Author(s)     : Asaf Porat          <asafpor1@post.tau.ac.il>

#ifndef LINE_THROUGH_PLOYTOPES_OVERLAY_H
#define LINE_THROUGH_PLOYTOPES_OVERLAY_H

#include <CGAL/basic.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

namespace CGAL {

template <typename Arrangement_2>
class Overlay_traits_update_face_count : 
    public CGAL::Arr_default_overlay_traits<Arrangement_2>
{
  typedef typename Arrangement_2::Face_const_handle     Face_handle_A;
  typedef typename Arrangement_2::Face_const_handle     Face_handle_B;
  typedef typename Arrangement_2::Face_handle           Face_handle_R;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_handle_A;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_handle_B;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle_R;


public:
  /* constructs the face f induced by the an overlap between the faces
     f1 and f2. */
  virtual void create_face (Face_handle_A f1,
                            Face_handle_B f2,
                            Face_handle_R f) const
  {
    if (f1->segs_begin() != f1->segs_end())
    {
      f->add_segment(*f1->segs_begin());
    }
    if (f2->segs_begin() != f2->segs_end())
    {
      f->add_segment(*f2->segs_begin());
    }
           
//     if (!f1->is_unbounded() && !f2->is_unbounded())
//     {
// #if LTS_POLY_OVERLAY_DEBUG
//       typename Arrangement_2::Ccb_halfedge_const_circulator chc = f1->outer_ccb();
//       typename Arrangement_2::Ccb_halfedge_const_circulator chc1 = chc;
//       do
//       {

//         std::cout << "F1 edge = " << chc1->curve() << std::endl;
//         chc1++;
//       }
//       while (chc1 != chc);
              
//       chc = f2->outer_ccb();
//       chc1 = chc;
//       do
//       {
//         std::cout << "F2 edge = " << chc1->curve() << std::endl;
//         chc1++;
//       }
//       while (chc1 != chc);
// #endif              
//       /* The face represents a polytope. */
//       f->set_plane_face_count(1);
//     }
//     else
//     {
//       f->set_plane_face_count(0);
//     }
  }

  /*!
   * Create an edge e that matches the overlap between e1 and e2.
   */
  virtual void create_edge (Halfedge_handle_A e1,
                            Halfedge_handle_B e2,
                            Halfedge_handle_R e) const
  {
    e->add_segment(*e1->segs_begin());
    e->add_segment(*e2->segs_begin());
    e->twin()->add_segment(*e1->segs_begin());
    e->twin()->add_segment(*e2->segs_begin());
  }
            
  /*!
   * Create an edge e that matches the edge e1, contained in the face f2.
   */
  virtual void create_edge (Halfedge_handle_A e1,
                            Face_handle_B f2,
                            Halfedge_handle_R e) const
  {
     e->add_segment(*e1->segs_begin());
     e->twin()->add_segment(*e1->segs_begin());
  }
            
  /*!
   * Create an edge e that matches the edge e2, contained in the face f1.
   */
  virtual void create_edge (Face_handle_A f1,
                            Halfedge_handle_B e2,
                            Halfedge_handle_R e) const
  {
    e->add_segment(*e2->segs_begin());
    e->twin()->add_segment(*e2->segs_begin());
  }
};

template <typename Arrangement_2,typename Observer>
class Lines_through_polytopes_overlay_add_bounded_faces
{
public:
  Lines_through_polytopes_overlay_add_bounded_faces()
  {
            
  }
         
  ~Lines_through_polytopes_overlay_add_bounded_faces()
  {
  }
  
  template <typename Traits_2_adapt, typename Polyhedron_3>
  void operator()(Traits_2_adapt& traits_2_adapt,
                  const Arrangement_2& arr_1,
                  const Arrangement_2& arr_2,
                  Arrangement_2& res_arr,
                  const Polyhedron_3& p)
  {
    Arrangement_2 overlay_arr;
    const Overlay_traits_update_face_count<Arrangement_2> OT;
    CGAL::overlay(arr_1, arr_2, overlay_arr, OT);
    std::list<typename Arrangement_2::Curve_2> arcs;
            
            
    typename Arrangement_2::Edge_iterator   eit;
    for (eit = overlay_arr.edges_begin(); 
         eit != overlay_arr.edges_end();
         ++eit)
    {
      if (eit->face()->num_of_overlap_plane_faces() == 0 ||
          eit->twin()->face()->num_of_overlap_plane_faces() == 0)
      {
        typename Arrangement_2::Curve_2 arc;
        /* Convert x monotone curve to regular curve. */
#if LTS_POLY_OVERLAY_DEBUG
        std::cout << "eit->curve() = " << eit->curve() << std::endl;
#endif
             
        traits_2_adapt.create_curve_on_plane_arr(arc,eit->curve());

        arcs.push_back(arc);
      }
    }

    insert(res_arr,arcs.begin(),arcs.end());
  }
      
};

template <typename Arrangement_2,typename Observer>
class Lines_through_polytopes_overlay_add_edges_in_box
{
public:
  Lines_through_polytopes_overlay_add_edges_in_box()
  {
  }
         
  ~Lines_through_polytopes_overlay_add_edges_in_box()
  {
  }
             
  template <typename Traits_2_adapt, typename Polyhedron_3>
  void operator()(Traits_2_adapt& traits_2_adapt,
                  const Arrangement_2& arr_1,
                  const Arrangement_2& arr_box,
                  Arrangement_2& res_arr,
                  const Polyhedron_3& p)
  {
    Arrangement_2 overlay_arr;
    const Overlay_traits_update_face_count<Arrangement_2> OT;
    CGAL::overlay(arr_1, arr_box, overlay_arr, OT);
    std::list<typename Arrangement_2::Curve_2> arcs;
            
    typename Arrangement_2::Edge_iterator   eit;
    for (eit = overlay_arr.edges_begin(); 
         eit != overlay_arr.edges_end();
         ++eit)
    {
      if (!eit->face()->is_unbounded() && !eit->twin()->face()->is_unbounded() &&
          /* The obj of the edges in the box is set to NULL. */
          *(eit->segs_begin()) == &p)
      {
        typename Arrangement_2::Curve_2 arc;
        /* Convert x monotone curve to regular curve. */
#if LTS_POLY_OVERLAY_DEBUG
        std::cout << "eit->curve() = " << eit->curve() << std::endl;
#endif
                  
        traits_2_adapt.create_curve_on_plane_arr(arc,eit->curve());
                  
        arcs.push_back(arc);
      }
    }
    insert(res_arr,arcs.begin(),arcs.end());
  }
};

} //namespace CGAL

#endif
