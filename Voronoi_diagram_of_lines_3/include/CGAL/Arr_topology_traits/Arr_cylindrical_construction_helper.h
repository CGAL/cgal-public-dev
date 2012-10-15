// Copyright (c) 2006-2007  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_on_surface_2/include/CGAL/Arr_topology_traits/Arr_cylindrical_construction_helper.h $
// $Id: Arr_cylindrical_construction_helper.h 41118 2007-12-07 14:25:32Z efif $
// 
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_CYLINDRICAL_CONSTRUCTION_HELPER_H
#define CGAL_ARR_CYLINDRICAL_CONSTRUCTION_HELPER_H

/*! \file
 * Definition of the Arr_cylindrical_construction_helper class-template.
 */

#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

/*! \class Arr_cylindrical_construction_helper
 * A helper class for the construction sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <class Traits_, class Arrangement_, class Event_, class Subcurve_> 
class Arr_cylindrical_construction_helper
{
public:
  typedef Traits_                                         Traits_2;
  typedef Arrangement_                                    Arrangement_2;
  typedef Event_                                          Event;
  typedef Subcurve_                                       Subcurve;

  typedef typename Traits_2::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits_2::Point_2                      Point_2;

  typedef Sweep_line_empty_visitor<Traits_2, Subcurve, Event>
                                                          Base_visitor;

  typedef typename Arrangement_2::Vertex_handle           Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement_2::Face_handle             Face_handle;
  
  typedef typename Subcurve::Halfedge_indices_list        Indices_list;
  typedef Unique_hash_map<Halfedge_handle, Indices_list>  Halfedge_indices_map;

protected:
  typedef typename Arrangement_2::Topology_traits         Topology_traits;

  typedef typename Topology_traits::Vertex                DVertex;
  typedef typename Topology_traits::Halfedge              DHalfedge;

  // Data members:

  //! The topology-traits class
  Topology_traits * m_top_traits;

  //! An arrangement accessor
  Arr_accessor<Arrangement_2> m_arr_access;

  //! A pointer to a map of halfedges to indices lists
  // (stored in the visitor class)
  Halfedge_indices_map * m_he_ind_map_p;

  //! The face with no outer ccb 
  Face_handle m_top_face;
  
  // The current left fictitious halfedge (on x = -oo).
  Halfedge_handle   m_l_vh;
  // The current left fictitious halfedge (on x = +oo).
  Halfedge_handle   m_r_vh;
  
  // The previous event at x = -oo, 
  // needed to comfortable split of fictitious edges at the left 
  Event*            m_prev_minus_inf_x_event; 
  

public:
  /*! Constructor. */
  Arr_cylindrical_construction_helper(Arrangement_2 * arr) :
    m_top_traits(arr->topology_traits()),
    m_arr_access(*arr),
    m_he_ind_map_p(NULL),
    m_prev_minus_inf_x_event (NULL)
  {}

  /*! Destructor. */
  virtual ~Arr_cylindrical_construction_helper() {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  virtual void before_sweep(){ 
    CGAL_assertion(m_top_traits->cylindrical_face() != NULL);
    CGAL_assertion(m_top_traits->l_fictitious_face() != NULL);
    CGAL_assertion(m_top_traits->r_fictitious_face() != NULL);
    CGAL_assertion(m_top_traits->l_fictitious_vertex() != NULL);
    CGAL_assertion(m_top_traits->r_fictitious_vertex() != NULL);

    // Get the face with no outer ccb 
    m_top_face = Face_handle(m_top_traits->cylindrical_face());
    // init the current left halfedge 
    m_l_vh =  Halfedge_handle(m_top_traits->l_fictitious_vertex()->halfedge());
    CGAL_assertion(m_l_vh->direction() != m_l_vh->twin()->direction());
    if(m_l_vh->direction() != ARR_RIGHT_TO_LEFT)
      m_l_vh = m_l_vh->twin();
    CGAL_assertion(
        m_l_vh->twin()->face() == Face_handle(m_top_traits->l_fictitious_face()));

    // init the current right halfedge 
    m_r_vh =  Halfedge_handle(m_top_traits->r_fictitious_vertex()->halfedge());
    CGAL_assertion(m_r_vh->direction() != m_r_vh->twin()->direction());
    if(m_r_vh->direction() != ARR_LEFT_TO_RIGHT)
      m_r_vh = m_r_vh->twin();
    CGAL_assertion(
        m_r_vh->twin()->face() == Face_handle(m_top_traits->r_fictitious_face()));
   



    
    // there is no prev event yet, since we are before sweep .-)
    m_prev_minus_inf_x_event = NULL;
    
    return; 
  }

  /*! A notification invoked before the sweep-line starts handling the given
   * event.
   * 
   * AFAIKS: 
   * - The function has to do nothing for an event that is not on the boundary. 
   * - For an event at ARR_LEFT_BOUNDARY/ARR_RIGHT_BOUNDARY the function must 
   *   create a new vertex and split the corresponding fictitious edge 
   * - For an event at ARR_TOP_BOUNDARY/ARR_BOTTOM_BOUNDARY the function checks 
   *   whether the vertex already exists, if not it creates a new one and 
   *   puts it into the event. Note, there are no edges that we have to split 
   *   since we do not handle curves on the identification line.  
   */
  virtual void before_handle_event(Event * event){ 
    std::cout << " before_handle_event begin " << std::endl;

    std::cout << " printing event inforamtion: " << std::endl;
    
    std::cout << " event.is_closed : " << event->is_closed()<< std::endl;
    std::cout << " event.is_on_bou : " << event->is_on_boundary()<< std::endl;
    std::cout << " event.ps_x      : " << event->parameter_space_in_x()<< std::endl;
    std::cout << " event.ps_y      : " << event->parameter_space_in_y()<< std::endl;
    std::cout << " event.n_l_curves: " << event->number_of_left_curves()<< std::endl;
    std::cout << " event.n_r_curves: " << event->number_of_right_curves()<< std::endl;
    
    std::cout << " event.curve     : " << event->curve()<< std::endl;
    if (event->is_closed())
      std::cout << " event.point     : " << (event->point()) << std::endl;
    else
      std::cout << " event.point     : " << " NONE " << std::endl;
    
    if(!event->is_on_boundary()) return;
 
    Arr_parameter_space ps_x = event->parameter_space_in_x();
    Arr_parameter_space ps_y = event->parameter_space_in_y();
    Arr_curve_end ind = 
      (event->number_of_right_curves()==1) ? ARR_MIN_END : ARR_MAX_END;
    const X_monotone_curve_2&  xc = event->curve();

    CGAL_assertion(ps_x == ARR_INTERIOR || ps_y == ARR_INTERIOR);   
    CGAL_assertion( 
        event->number_of_left_curves()  +
        event->number_of_right_curves() == 1);
    
    if(ps_y == ARR_INTERIOR){
      CGAL_assertion(ps_x == 
          (ind == ARR_MIN_END)?ARR_LEFT_BOUNDARY:ARR_RIGHT_BOUNDARY);
      Vertex_handle v_at_inf =
        m_arr_access.create_boundary_vertex (xc, ind, ps_x, ps_y, false);
      
      if(ps_x == ARR_LEFT_BOUNDARY){
        m_arr_access.split_fictitious_edge(m_l_vh, v_at_inf);
        event->set_halfedge_handle(m_l_vh);
      
        // Update the incident halfedge of the previous vertex at x = -oo
        // (m_l_vh used to be incident to it, but now we have split it).
        if (m_prev_minus_inf_x_event != NULL)
          m_prev_minus_inf_x_event->set_halfedge_handle(m_l_vh->next());
        m_prev_minus_inf_x_event = event;
        return;
      }else{
        // The event lies on the right fictitious halfedge.
        m_arr_access.split_fictitious_edge(m_r_vh, v_at_inf);
        event->set_halfedge_handle(m_r_vh);
        m_r_vh = m_r_vh->next();
        return;
      } 
    }

    assert(event->parameter_space_in_y() != ARR_TOP_BOUNDARY); // TODO 
    assert(event->parameter_space_in_y() != ARR_BOTTOM_BOUNDARY); // TODO
    std::cout << " before_handle_event end " << std::endl;
    std::cout << " =============================== " << std::endl;
    return; 
  }

  /*! A notification invoked when a new subcurve is created. */
  virtual void add_subcurve(Halfedge_handle he, Subcurve * sc) { return; }

  /*! Collect a subcurve index that does not see any status-line from below.
   */
  void add_subcurve_in_top_face(unsigned int index){ return; }

  /*! A notification invoked before the given event it deallocated. */
  void before_deallocate_event(Event * event) { return; }
  //@} 
  
  /*! Set the map that maps each halfedge to the list of subcurve indices
   * that "see" the halfedge from below.
   */
  void set_halfedge_indices_map(Halfedge_indices_map & table)
  {
    m_he_ind_map_p = &table;
    return;
  }

  /*! Determine if we should swap the order of predecessor halfedges when
   * calling insert_at_vertices_ex() .
   */
  bool swap_predecessors (Event * event) const
  {
    // no idea what is going on here (Michael)
    return (event->parameter_space_in_x() == ARR_INTERIOR &&
            event->parameter_space_in_y() == ARR_TOP_BOUNDARY);
  }

  /*! Get the current top face. */
  Face_handle top_face() const {
    CGAL_assertion(m_top_face->number_of_outer_ccbs()==0);
    return m_top_face; 
  }
};

} //namespace CGAL

#endif
