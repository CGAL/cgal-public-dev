// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//

#ifndef CGAL_ARR_SPHERICAL_VERT_DECOMP_HELPER_H
#define CGAL_ARR_SPHERICAL_VERT_DECOMP_HELPER_H

/*! \file
 * Definition of the Arr_spherical_vert_decomp_helper class-template.
 */

namespace CGAL {

#include <CGAL/Sweep_line_empty_visitor.h>

/*! \class Arr_spherical_vert_decomp_helper
 * A helper class for the vertical decomposition sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <class Traits_, class Arrangement_>
class Arr_spherical_vert_decomp_helper
{
public:

  typedef Traits_                                      Traits_2;
  typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
  typedef Arrangement_                                 Arrangement_2;

  typedef typename Arrangement_2::Face_const_handle    Face_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle  Vertex_const_handle;

  typedef Sweep_line_empty_visitor<Traits_2>           Base_visitor;
  typedef typename Base_visitor::Event                 Event;
  typedef typename Base_visitor::Subcurve              Subcurve;

protected:

  typedef typename Arrangement_2::Topology_traits      Topology_traits;

  const Topology_traits*  m_top_traits;        // The topology traits.
  Vertex_const_handle     m_north_pole;        // The north pole.
  bool                    m_valid_north_pole;  // Is this a valid vertex.
  Face_const_handle       m_above_event_face;  // Current face above event.
  Vertex_const_handle     m_south_pole;        // The south pole.
  bool                    m_valid_south_pole;  // Is this a valid vertex.
  Face_const_handle       m_below_event_face;  // Current face below event.

public:

  /*! Constructor.
   * \param arr The arrangement.
   */
  Arr_spherical_vert_decomp_helper(const Arrangement_2* arr) :
    m_top_traits(arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /*! A notification issued before the sweep process starts. */
  void before_sweep();

  /*!
   * A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  void after_handle_event(Event* event);
  //@}

  /*! Get the current top object. */
  CGAL::Object top_object() const
  {
    return (m_valid_north_pole) ?
      CGAL::make_object(m_north_pole) : CGAL::make_object(m_above_event_face);
  }

  /*! Get the current bottom object. */
  CGAL::Object bottom_object() const
  {
    return (m_valid_south_pole) ?
      CGAL::make_object(m_south_pole) : CGAL::make_object(m_below_event_face);
  }
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <class Tr, class Arr> 
void Arr_spherical_vert_decomp_helper<Tr, Arr>::before_sweep()
{
  // Get the north pole and the face that intially contains it.
  m_valid_north_pole = (m_top_traits->north_pole() != NULL);
  if (m_valid_north_pole)
    m_north_pole = Vertex_const_handle(m_top_traits->north_pole());

  // initialize face above with top_face; it is updated at "after_handle_event"
  // TODO EBEF use "top_face()"
  m_above_event_face = Face_const_handle(m_top_traits->spherical_face());

  // Get the south pole and the face that intially contains it.
  m_valid_south_pole = (m_top_traits->south_pole() != NULL);
  if (m_valid_south_pole)
    m_south_pole = Vertex_const_handle (m_top_traits->south_pole());

  // initialize face below with bottom_face; it is updated at "after_handle_event"
  // TODO EBEF use "bottom_face()"
  m_below_event_face = Face_const_handle(m_top_traits->min_face());
}

//-----------------------------------------------------------------------------
// A notification invoked after the sweep-line finishes handling the given
// event.
///
template <class Tr, class Arr>
void
Arr_spherical_vert_decomp_helper<Tr, Arr>::after_handle_event(Event* event)
{
  // Ignore events that are not incident to the poles.
  if (event->parameter_space_in_y() == ARR_INTERIOR)
    return;

  // The is exactly one curve incident to an event with boundary conditions.
  // Obtain this curve and check whether it already exists in the arrangement.
  CGAL_assertion(((event->number_of_left_curves() == 0) &&
                  (event->number_of_right_curves() == 1)) ||
                 ((event->number_of_left_curves() == 1) &&
                  (event->number_of_right_curves() == 0)));

  const Arr_curve_end   ind =
    ((event->number_of_left_curves() == 0) &&
     (event->number_of_right_curves() == 1)) ? ARR_MIN_END : ARR_MAX_END;
  const X_monotone_curve_2& xc = (ind == ARR_MIN_END) ?
    (*(event->right_curves_begin()))->last_curve() :
    (*(event->left_curves_begin()))->last_curve();

  // update faces below and above event:

  if (event->parameter_space_in_y() == ARR_TOP_BOUNDARY)
  {
    // The event is incident to the north pole: update the face above the event
      m_above_event_face = (ind == ARR_MIN_END) ?
        xc.halfedge_handle()->twin()->face() :
        xc.halfedge_handle()->face();
  }
  else if (event->parameter_space_in_y() == ARR_BOTTOM_BOUNDARY)
  {
    // The event is incident to the south pole: update the face below the event
      m_below_event_face = (ind == ARR_MIN_END) ?
        xc.halfedge_handle()->face() :
        xc.halfedge_handle()->twin()->face();
  }
}

} //namespace CGAL

#endif
