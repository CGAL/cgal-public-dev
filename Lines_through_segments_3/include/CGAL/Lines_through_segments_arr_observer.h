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

template <typename Ext_obj, typename Arrangement_2>
class Lines_through_segments_arr_observer : 
    public CGAL::Arr_observer<Arrangement_2>
{
  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle     Halfedge_handle;
  typedef typename Arrangement_2::Face_handle         Face_handle;
  typedef typename Arrangement_2::X_monotone_curve_2  X_monotone_curve_2;
  typedef typename Arrangement_2::Point_2             Point_2;
   
private:
   const X_monotone_curve_2* create_new_edge_curve; //save the last inserted curve at before_create_edge.
   const Ext_obj* m_last_inserted_segment;
  std::list<const Ext_obj*> m_last_splitted_edge_segment_list;
  bool m_is_last_plane;
      
public:
  ~Lines_through_segments_arr_observer()
  {
  }

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

  void set_last_inserted_segment(const Ext_obj* s)
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


  /*!
   * Notification before the creation of a new vertex.
   * \param p The point to be associated with the vertex.
   *          This point cannot lies on the surface boundaries.
   */
  virtual void before_create_vertex (const Point_2&  p)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification after the creation of a new vertex.
   * \param v A handle to the created vertex.
   */
  virtual void after_create_vertex (Vertex_handle v)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification before the creation of a new boundary vertex.
   * \param cv The curve incident to the surface boundary.
   * \param ind The relevant curve-end.
   * \param ps_x The boundary condition of the vertex in x.
   * \param ps_y The boundary condition of the vertex in y.
   */
  virtual void before_create_boundary_vertex (const X_monotone_curve_2& /*cv*/,
                                              Arr_curve_end /* ind */,
                                              Arr_parameter_space /* ps_x */,
                                              Arr_parameter_space /* ps_y */)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification after the creation of a new vertex at infinity.
   * \param v A handle to the created vertex.
   */
  virtual void after_create_boundary_vertex (Vertex_handle v)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}



  /*!
   * Notification before the modification of an existing vertex.
   * \param v A handle to the vertex to be updated.
   * \param p The point to be associated with the vertex.
   */
  virtual void before_modify_vertex (Vertex_handle v,
                                     const Point_2& p)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification after a vertex was modified.
   * \param v A handle to the updated vertex.
   */
  virtual void after_modify_vertex (Vertex_handle v)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}



  /*!
   * Notification before the splitting of a fictitious edge into two.
   * \param e A handle to one of the existing halfedges.
   * \param v A vertex representing the unbounded split point.
   */
  virtual void before_split_fictitious_edge (Halfedge_handle e,
                                             Vertex_handle v)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification after a fictitious edge was split.
   * \param e1 A handle to one of the twin halfedges forming the first edge.
   * \param e2 A handle to one of the twin halfedges forming the second edge.
   */
  virtual void after_split_fictitious_edge (Halfedge_handle /* e1 */,
                                            Halfedge_handle /* e2 */)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification before the splitting of a face into two.
   * \param f A handle to the existing face.
   * \param e The new edge whose insertion causes the face to split.
   */
  virtual void before_split_face (Face_handle /* f */,
                                  Halfedge_handle e)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification before the splitting of an outer CCB into two.
  //  * \param f A handle to the face that owns the outer CCB.
  //  * \param h A circulator representing the component boundary.
  //  * \param e The new edge whose removal causes the outer CCB to split.
  //  */
  // virtual void before_split_outer_ccb (Face_handle /* f */,
  //                                      Ccb_halfedge_circulator /* h */,
  //                                      Halfedge_handle e)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification after an outer CCB was split.
  //  * \param f A handle to the face that owns the outer CCBs.
  //  * \param h1 A circulator representing the boundary of the first component.
  //  * \param h2 A circulator representing the boundary of the second component.
  //  */
  // virtual void after_split_outer_ccb (Face_handle /* f */,
  //                                     Ccb_halfedge_circulator /* h1 */,
  //                                     Ccb_halfedge_circulator /* h2 */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification before the splitting of an inner CCB into two.
  //  * \param f A handle to the face containing the inner CCB.
  //  * \param h A circulator representing the component boundary.
  //  * \param e The new edge whose removal causes the inner CCB to split.
  //  */
  // virtual void before_split_inner_ccb (Face_handle /* f */,
  //                                      Ccb_halfedge_circulator /* h */,
  //                                      Halfedge_handle e)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification after an inner CCB was split.
  //  * \param f A handle to the face containing the inner CCBs.
  //  * \param h1 A circulator representing the boundary of the first component.
  //  * \param h2 A circulator representing the boundary of the second component.
  //  */
  // virtual void after_split_inner_ccb (Face_handle /* f */,
  //                                     Ccb_halfedge_circulator /* h1 */,
  //                                     Ccb_halfedge_circulator /* h2 */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification before the creation of a new outer CCB of a face.
   * \param f A handle to the face that owns the outer CCB.
   * \param e A halfedge along the new outer CCB.
   */
  virtual void before_add_outer_ccb (Face_handle /* f */,
                                     Halfedge_handle e)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification after an outer CCB was added to a face.
  //  * \param h A circulator representing the boundary of the new outer CCB.
  //  */
  // virtual void after_add_outer_ccb (Ccb_halfedge_circulator /* h */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification before the creation of a new inner CCB inside a face.
   * \param f A handle to the face containing the inner CCB.
   * \param e The new halfedge that forms the new inner CCB.
   */
  virtual void before_add_inner_ccb (Face_handle /* f */,
                                     Halfedge_handle e)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification after an inner CCB was created inside a face.
  //  * \param h A circulator representing the boundary of the new inner CCB.
  //  */
  // virtual void after_add_inner_ccb (Ccb_halfedge_circulator /* h */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification before the creation of a new isolated vertex inside a face.
   * \param f A handle to the face containing the isolated vertex.
   * \param v The isolated vertex.
   */
  virtual void before_add_isolated_vertex (Face_handle /* f */,
                                           Vertex_handle v)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification after an isolated vertex was created inside a face.
   * \param v The isolated vertex.
   */
  virtual void after_add_isolated_vertex (Vertex_handle v)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification before the merging of two edges.
   * \param e1 A handle to one of the halfedges forming the first edge.
   * \param e2 A handle to one of the halfedges forming the second edge.
   * \param c The x-monotone curve to be associated with the merged edge.
   */
  virtual void before_merge_edge (Halfedge_handle /* e1 */,
                                  Halfedge_handle /* e2 */,
                                  const X_monotone_curve_2& /* c */)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification after an edge was merged.
   * \param e A handle to one of the twin halfedges forming the merged edge.
   */
  virtual void after_merge_edge (Halfedge_handle e)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification before the merging of two fictitious edges.
   * \param e1 A handle to one of the halfedges forming the first edge.
   * \param e2 A handle to one of the halfedges forming the second edge.
   */
  virtual void before_merge_fictitious_edge (Halfedge_handle /* e1 */,
                                             Halfedge_handle /* e2 */)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification after a fictitious edge was merged.
   * \param e A handle to one of the twin halfedges forming the merged edge.
   */
  virtual void after_merge_fictitious_edge (Halfedge_handle e)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification before the merging of two faces.
   * \param f1 A handle to the first face.
   * \param f2 A handle to the second face.
   * \param e The edge whose removal causes the faces to merge.
   */
  virtual void before_merge_face (Face_handle /* f1 */,
                                  Face_handle /* f2 */,
                                  Halfedge_handle e)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification after a face was merged.
   * \param f A handle to the merged face.
   */
  virtual void after_merge_face (Face_handle /* f */)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification before the merging of two outer CCBs.
  //  * \param f A handle to the face that owns the outer CCBs.
  //  * \param h1 A circulator representing the boundary of the first component.
  //  * \param h2 A circulator representing the boundary of the second component.
  //  * \param e The edge whose insertion or removal causes the CCBs to merge.
  //  */
  // virtual void before_merge_outer_ccb (Face_handle /* f */,
  //                                      Ccb_halfedge_circulator /* h1 */,
  //                                      Ccb_halfedge_circulator /* h2 */,
  //                                      Halfedge_handle e)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification after an outer CCB was merged.
  //  * \param f A handle to the face that owns the outer CCBs.
  //  * \param h A circulator representing the boundary of the merged component.
  //  */
  // virtual void after_merge_outer_ccb (Face_handle /* f */,
  //                                     Ccb_halfedge_circulator /* h */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification before the merging of two inner CCBs (holes).
  //  * \param f A handle to the face that contains the inner CCBs.
  //  * \param h1 A circulator representing the boundary of the first component.
  //  * \param h2 A circulator representing the boundary of the second component.
  //  * \param e The edge whose insertion causes the inner CCBs to merge.
  //  */
  // virtual void before_merge_inner_ccb (Face_handle /* f */,
  //                                      Ccb_halfedge_circulator /* h1 */,
  //                                      Ccb_halfedge_circulator /* h2 */,
  //                                      Halfedge_handle e)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification after an inner CCB was merged.
  //  * \param f A handle to the face that contains the inner CCBs.
  //  * \param h A circulator representing the boundary of the merged component.
  //  */
  // virtual void after_merge_inner_ccb (Face_handle /* f */,
  //                                     Ccb_halfedge_circulator /* h */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification before an outer CCB is moved from one face to another.
  //  * \param from_f A handle to the face that currently owns the outer CCB.
  //  * \param to_f A handle to the face that should own the outer CCB.
  //  * \param h A circulator representing the boundary of the component.
  //  */
  // virtual void before_move_outer_ccb (Face_handle /* from_f */,
  //                                     Face_handle /* to_f */,
  //                                     Ccb_halfedge_circulator /* h */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification after an outer CCB is moved from one face to another.
  //  * \param h A circulator representing the boundary of the component.
  //  */
  // virtual void after_move_outer_ccb (Ccb_halfedge_circulator /* h */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}


  // /*!
  //  * Notification before an inner CCB is moved from one face to another.
  //  * \param from_f A handle to the face currently containing the inner CCB.
  //  * \param to_f A handle to the face that should contain the inner CCB.
  //  * \param h A circulator representing the boundary of the component.
  //  */
  // virtual void before_move_inner_ccb (Face_handle /* from_f */,
  //                                     Face_handle /* to_f */,
  //                                     Ccb_halfedge_circulator /* h */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification after an inner CCB is moved from one face to another.
  //  * \param h A circulator representing the boundary of the component.
  //  */
  // virtual void after_move_inner_ccb (Ccb_halfedge_circulator /* h */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification before an isolated vertex is moved from one face to another.
   * \param from_f A handle to the face currently containing the vertex.
   * \param to_f A handle to the face that should contain the vertex.
   * \param v The isolated vertex.
   */
  virtual void before_move_isolated_vertex (Face_handle /* from_f */,
                                            Face_handle /* to_f */,
                                            Vertex_handle v)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification after an isolated vertex is moved from one face to another.
   * \param v The isolated vertex.
   */
  virtual void after_move_isolated_vertex (Vertex_handle v)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notificaion before the removal of a vertex.
   * \param v A handle to the vertex to be deleted.
   */
  virtual void before_remove_vertex (Vertex_handle v)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notificaion after the removal of a vertex.
   */
  virtual void after_remove_vertex ()
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notification before the removal of an edge.
   * \param e A handle to one of the twin halfedges to be deleted.
   */
  virtual void before_remove_edge (Halfedge_handle e)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notificaion after the removal of an edge.
   */
  virtual void after_remove_edge ()
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification before the removal of an outer CCB.
  //  * \param f The face that owns the outer CCB.
  //  * \param h A circulator representing the boundary of the component.
  //  */
  // virtual void before_remove_outer_ccb (Face_handle /* f */,
  //                                       Ccb_halfedge_circulator /* h */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notificaion after the removal of an outer CCB.
   * \param f The face that used to own the outer CCB.
   */
  virtual void after_remove_outer_ccb (Face_handle /* f */)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  // /*!
  //  * Notification before the removal of an inner CCB.
  //  * \param f The face containing the inner CCB.
  //  * \param h A circulator representing the boundary of the component.
  //  */
  // virtual void before_remove_inner_ccb (Face_handle /* f */,
  //                                       Ccb_halfedge_circulator /* h */)
  // {/*std::cout << "line = " << __LINE__ << std::endl;*/}

  /*!
   * Notificaion after the removal of an inner CCB.
   * \param f The face that used to contain the inner CCB.
   */
  virtual void after_remove_inner_ccb (Face_handle /* f */)
  {/*std::cout << "line = " << __LINE__ << std::endl;*/}


};

} //namespace CGAL

#endif //LINE_THROUGH_SEGMENTS_ARR_OBSERVER_H
