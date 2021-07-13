// Copyright (c) 2007,2009,2010,2011,2013,2014 Max-Planck-Institute Saarbruecken (Germany), Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//                 Eric Berberich <eric.berberich@cgal.org>

#ifndef CGAL_ARR_TOROIDAL_BATCHED_PL_HELPER_H
#define CGAL_ARR_TOROIDAL_BATCHED_PL_HELPER_H

/*! \file
 * Definition of the Arr_toroidal_batched_pl_helper class-template.
 */

namespace CGAL {

/*! \class Arr_toroidal_batched_pl_helper
 * A helper class for the batched point-location sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class.
 */
template <class Traits_, class Arrangement_, typename Event_,
          typename Subcurve_>
class Arr_toroidal_batched_pl_helper
{
public:

  typedef Traits_                                      Traits_2;
  typedef Arrangement_                                 Arrangement_2;
  typedef Event_                                       Event;
  typedef Subcurve_                                    Subcurve;

  typedef typename Arrangement_2::Face_const_handle    Face_const_handle;

protected:

  typedef typename Arrangement_2::Topology_traits      Topology_traits;

  // Data members:
  //! The topology-traits class.
  const Topology_traits * m_top_traits;

  //! The unbounded arrangement face.
  Face_const_handle m_top_face;

public:
  /*! Constructor.
   * \param arr The arrangement.
   */
  Arr_toroidal_batched_pl_helper(const Arrangement_2 *arr) :
    m_top_traits(arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /*! A notification issued before the sweep process starts. */
  void before_sweep()
  {
    m_top_face = Face_const_handle(m_top_traits->reference_face());
  }

  /*! A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  void after_handle_event(Event * ) { return; }
  //@}

  /*! Get the current top face. */
  Face_const_handle top_face() const
  {
    return m_top_face;
  }
};

} //namespace CGAL

#endif // CGAL_ARR_TOROIDAL_BATCHED_PL_HELPER_H
