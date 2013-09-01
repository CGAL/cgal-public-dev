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
// $URL$
// $Id$
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>

/*! \file
 * Definition of the internal _Arr_default_overlay_traits_base class template.
 */

#ifndef CGAL_ARR_RED_BLUE_OVERLAY_TRAITS_BASE_H
#define CGAL_ARR_RED_BLUE_OVERLAY_TRAITS_BASE_H

namespace CGAL {

/*!
 * \class
 * An overlay-traits class for computing the overlay of two arrangement that
 * are templated with the default DCEL classes, namely they store extra
 * data for red-blue-overlay with 
 * their DCEL features. T
 */
template <class ArrangementA, class ArrangementB, class ArrangementR>
class Arr_red_blue_overlay_traits_base {
public:

    // TODO
    // ArrangementR must be of type
#if 0
    CGAL::Arrangement_2< ArrangementA::Geo_traits, 
    CGAL::Arr_extended_dcel< ArrangementA::Geo_traits, 
                             std::pair< bool, bool >, 
                             std::pair< bool, bool >, 
                             std::pair< bool, bool > >
    // the bool bool is important!
#endif

    typedef typename ArrangementA::Vertex_const_handle    Vertex_handle_A;
    typedef typename ArrangementA::Halfedge_const_handle  Halfedge_handle_A;
    typedef typename ArrangementA::Face_const_handle      Face_handle_A;
    
    typedef typename ArrangementB::Vertex_const_handle    Vertex_handle_B;
    typedef typename ArrangementB::Halfedge_const_handle  Halfedge_handle_B;
    typedef typename ArrangementB::Face_const_handle      Face_handle_B;
    
    typedef typename ArrangementR::Vertex_handle          Vertex_handle_R;
    typedef typename ArrangementR::Halfedge_handle        Halfedge_handle_R;
    typedef typename ArrangementR::Face_handle            Face_handle_R;
    
    /*! Destructor. */
    virtual ~Arr_red_blue_overlay_traits_base ()
    {}
    
    /*!
     * Create a vertex v that corresponds to the coinciding vertices v1 and v2.
     */
    virtual void create_vertex (Vertex_handle_A v1,
                                Vertex_handle_B v2,
                                Vertex_handle_R v) const
    {
        v->data().first = true;
        v->data().second = true;
    }
    
    /*!
     * Create a vertex v that mathces v1, which lies of the edge e2.
     */
    virtual void create_vertex (Vertex_handle_A v1,
                                Halfedge_handle_B e2,
                                Vertex_handle_R v) const
    {
        v->data().first = true;
        v->data().second = true;
    }
    
    /*!
     * Create a vertex v that mathces v1, contained in the face f2.
     */
    virtual void create_vertex (Vertex_handle_A v1,
                                Face_handle_B f2,
                                Vertex_handle_R v) const
    {
        v->data().first = true;
        v->data().second = false;
    }
    
    /*!
     * Create a vertex v that mathces v2, which lies of the edge e1.
     */
    virtual void create_vertex (Halfedge_handle_A e1,
                                Vertex_handle_B v2,
                                Vertex_handle_R v) const
    {
        v->data().first = true;
        v->data().second = true;
    }
    
    /*!
     * Create a vertex v that mathces v2, contained in the face f1.
     */
    virtual void create_vertex (Face_handle_A f1,
                                Vertex_handle_B v2,
                                Vertex_handle_R v) const
    {
        v->data().first = false;
        v->data().second = true;
    }
    
    /*!
     * Create a vertex v that mathces the intersection of the edges e1 and e2.
     */
    virtual void create_vertex (Halfedge_handle_A e1,
                                Halfedge_handle_B e2,
                                Vertex_handle_R v) const
    {
        v->data().first = true;
        v->data().second = true;
    }
    
    /*!
     * Create an edge e that matches the overlap between e1 and e2.
     */
    virtual void create_edge (Halfedge_handle_A e1,
                              Halfedge_handle_B e2,
                              Halfedge_handle_R e) const
    {
        e->data().first = true;
        e->data().second = true;
        e->twin()->data() = e->data();
    }
    
    /*!
     * Create an edge e that matches the edge e1, contained in the face f2.
     */
    virtual void create_edge (Halfedge_handle_A e1,
                              Face_handle_B f2,
                              Halfedge_handle_R e) const
    {
        e->data().first = true;
        e->data().second = false;
        e->twin()->data() = e->data();
    }
    
    /*!
     * Create an edge e that matches the edge e2, contained in the face f1.
     */
    virtual void create_edge (Face_handle_A f1,
                              Halfedge_handle_B e2,
                              Halfedge_handle_R e) const
    {
        e->data().first = false;
        e->data().second = true;
        e->twin()->data() = e->data();
    }
    
    /*!
     * Create a face f that matches the overlapping region between f1 and f2.
     */
    virtual void create_face (Face_handle_A f1,
                              Face_handle_B f2,
                              Face_handle_R f) const
    {
        f->data().first = false;
        f->data().second = false;
    }
};

} //namespace CGAL

#endif
