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

#ifndef LINE_THROUGH_SEGMENTS_ARR_OBSERVER_H
#define LINE_THROUGH_SEGMENTS_ARR_OBSERVER_H

#include <CGAL/basic.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Lines_through_segments_general_functions.h>

/*************************************************************
 * An arrangement observer.
 * Used to receive notifications of creation of new edges vertexes and faces.
 * For each new edge associate it with the line that created it.
 * For each new split face set the number of planes that created it.
 *
 **************************************************************/

namespace CGAL {

template <typename Arrangement_2>
class Lines_through_segments_arr_observer 
  : public CGAL::Arr_observer<Arrangement_2>
{
  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle     Halfedge_handle;
  typedef typename Arrangement_2::Face_handle         Face_handle;
  typedef typename Arrangement_2::X_monotone_curve_2  X_monotone_curve_2;
  typedef typename Arrangement_2::Point_2             Point_2;
  typedef typename Arrangement_2::Geometry_traits_2::Curve_data Ext_obj;
   
private:
  const X_monotone_curve_2* create_new_edge_curve; //save the last inserted curve at before_create_edge.
  Ext_obj m_last_inserted_segment;
  std::list<Ext_obj> m_last_splitted_edge_segment_list;
  bool m_is_last_plane;
      
public:
  Lines_through_segments_arr_observer():
    CGAL::Arr_observer<Arrangement_2> ()
  {
  }
      
  Lines_through_segments_arr_observer (Arrangement_2& arr) :
    CGAL::Arr_observer<Arrangement_2> (arr)
  {
    CGAL_precondition (arr.is_empty());
    m_is_last_plane = false;
  }

  void set_is_plane(bool _is_plane)
  {
    m_is_last_plane = _is_plane;
  }

  void set_last_inserted_segment(Ext_obj s)
  {
    m_last_inserted_segment = s;
  }
            
  /*
   * issued just before a new edge that corresponds to the x-monotone curve c
   * and connects the vertices v1 and v2 is created.
   */
  virtual void before_create_edge(const X_monotone_curve_2& c, 
                                  Vertex_handle v1, 
                                  Vertex_handle v2)
  {
     create_new_edge_curve = &c;
#if OBSERVER_PRINTS
    std::cout << change_color(CGAL_BLINK,"before_create_edge") << std::endl;
    std::cout << "c = " << c << std::endl;
    typename Arrangement_2::Geometry_traits_2::Data_container::const_iterator       dit;

     for (dit = c.data().begin(); dit != c.data().end();
          ++dit)
     {
        std::cout << "c.data() = " << **dit << std::endl;
     }
#endif
  }
   
  /* issued just before an edge e is modified to be associated 
     with the x-monotone curve c.
  */
  virtual void before_modify_edge(Halfedge_handle e, 
                                  const X_monotone_curve_2& c)
  {
     typename Arrangement_2::Geometry_traits_2::Data_container::const_iterator       dit;
     for (dit = c.data().begin(); dit != c.data().end();
          ++dit)
     {
        e->add_segment(*dit);
        e->twin()->add_segment(*dit);
     }

#if OBSERVER_PRINTS
    std::cout << change_color(CGAL_BLINK,"before_modify_edge") << std::endl;
    std::cout << "e = (" << e->curve().source() << "," << e->curve().target()
              << ")" << std::endl;
#endif
  }
   
  /*
   * issued immediately after an existing edge e has been modified.
   */
  virtual void after_modify_edge(Halfedge_handle e)
  {
#if OBSERVER_PRINTS
    std::cout << change_color(CGAL_BLINK,"after_modify_edge") << std::endl;
    std::cout << "e = (" << e->curve().source() << "," << e->curve().target()
              << ")" << std::endl;
#endif
  }
     
  /* 
   *  Issued immediately after a new edge e has been created. 
   */

  virtual void after_create_edge(Halfedge_handle e)
  {
#if OBSERVER_PRINTS
    std::cout << change_color(CGAL_BLINK,"after_create_edge") << std::endl;
    std::cout << "curve = (" << e->curve() << ")" << std::endl;
#endif
    typename Arrangement_2::Geometry_traits_2::Data_container::const_iterator       dit;
    for (dit = create_new_edge_curve->data().begin(); dit != create_new_edge_curve->data().end();
         ++dit)
    {
       e->add_segment(*dit);
       e->twin()->add_segment(*dit);
    }
    create_new_edge_curve = NULL;
  }

  /* issued just before an edge e is split into two edges that should be
   * associated with the x-monotone curves c1 and c2.
   * The vertex v corresponds to the split point, and will be used to separate
   * the two resulting edges. 
   */
  virtual void before_split_edge(Halfedge_handle e,
                                 Vertex_handle v,
                                 const X_monotone_curve_2& c1,
                                 const X_monotone_curve_2& c2)
  {
#if OBSERVER_PRINTS
    std::cout << change_color(CGAL_BLINK,"before split edge") << std::endl;
    std::cout << c1 << std::endl;
    typename Arrangement_2::Geometry_traits_2::Data_container::const_iterator       dit;

     for (dit = c1.data().begin(); dit != c1.data().end();
          ++dit)
     {
        std::cout << "c1.data() = " << **dit << std::endl;
     }

//    std::cout << *(c1.data())  << std::endl;
    std::cout << c2  << std::endl;
    dit;

     for (dit = c2.data().begin(); dit != c2.data().end();
          ++dit)
     {
        std::cout << "c2.data() = " << *dit << std::endl;
     }
//  std::cout << *(c2.data())  << std::endl;
#endif
    m_last_splitted_edge_segment_list = e->get_segments_list();
  }
      
  /* Issued immediately after an existing edge has been split into the two
   * given edges e1 and e2.
   */
  virtual void after_split_edge(Halfedge_handle e1, Halfedge_handle e2)
  {
#if OBSERVER_PRINTS
    std::cout << change_color(CGAL_BLINK,"after_split_edge") << std::endl;
#endif

    e1->set_segments_list(m_last_splitted_edge_segment_list);
    e1->twin()->set_segments_list(m_last_splitted_edge_segment_list);

    e2->set_segments_list(m_last_splitted_edge_segment_list);
    e2->twin()->set_segments_list(m_last_splitted_edge_segment_list);
  }
      
  /* Issued immediately after the existing face f1 has been split,
     such that a portion of it now forms a new face f2.
     The flag is_hole designates whether f2 forms a hole inside f1. */

  virtual void after_split_face(Face_handle f1, Face_handle f2, bool is_hole)
  {
#if OBSERVER_PRINTS
    std::cout << change_color(CGAL_CYAN,"after_split_face") << std::endl;
#endif
    if (m_is_last_plane)
    {
#if OBSERVER_PRINTS
       std::cout << "f2->add_segment " << *m_last_inserted_segment << std::endl;
#endif
       f2->add_segment(m_last_inserted_segment);
    }
  }
};

} //namespace CGAL

#endif //LINE_THROUGH_SEGMENTS_ARR_OBSERVER_H
