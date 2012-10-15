// Copyright (c) 2008  Max-Planck-Institute fuer Informatik (Germany), 
// Tel-Aviv University (Israel).
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
// $Date$
// 
//
// Author(s)     : Eric Berberich <ericb@post.tau.ac.il>

#ifndef CGAL_ARR_GEOMETRY_OF_DCEL_RECORD_2_H
#define CGAL_ARR_GEOMETRY_OF_DCEL_RECORD_2_H

/*! \file
 * Definition of a functor object returning geometric objects of dcel-records
 */

#include <CGAL/config.h>

namespace CGAL {

/*! \class
 * Functor object that returns geometric objects belonging to different 
 * dcel-records
 */
template <class Arrangement_2_>
class Arr_geometry_of_dcel_record_2 {

public:

    //! this instance's template parameter
    typedef Arrangement_2_ Arrangement_2;
    
    //! type of x-monotone curve
    typedef typename Arrangement_2::X_monotone_curve_2 X_monotone_curve_2;
    //! type of geometric point
    typedef typename Arrangement_2::Point_2 Point_2;
    
    typedef typename Arrangement_2::Vertex_const_handle Vertex_const_handle;
    typedef typename Arrangement_2::Halfedge_const_handle 
    Halfedge_const_handle;
    typedef typename Arrangement_2::Face_const_handle Face_const_handle;
    
    typedef typename Arrangement_2::Outer_ccb_const_iterator 
    Ccb_const_iterator;

#if 0
    // TODO check that both are identical
    typedef typename Arrangement_2::Inner_ccb_const_iterator 
    Inner_ccb_const_iterator;
#endif

    typedef typename Arrangement_2::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
    
    //typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
    //Halfedge_around_vertex_const_circulator;


    //!\name Operators
    //!@{
    
    template < class OutputIterator >
    OutputIterator operator()(Vertex_const_handle vh, 
                              OutputIterator oi) const {
        
        // value-type is CGAL::Object
        if (!vh->is_fictitious()) {
            *oi++ = CGAL::make_object(vh->point()); 
        }
        
        return oi;
    }


    template < class OutputIterator >
    OutputIterator operator()(Halfedge_const_handle heh,
                              OutputIterator oi) const {
        
        // value-type is CGAL::Object
        if (!heh->is_fictitious()) {
            *oi++ = CGAL::make_object(heh->curve()); 
        }
        
        return oi;
    }


#if 0
    template < class OutputIterator >
    OutputIterator operator()(Outer_ccb_const_iterator occb, 
                              OutputIterator oi) const {
        
        // value-type is CGAL::Object
        
        
        return oi;
    }
#endif

    template < class OutputIterator >
    OutputIterator operator()(Ccb_const_iterator iccb, 
                              OutputIterator oi) const {
        
        // value-type is CGAL::Object

        std::list< Halfedge_const_handle > seen;
        
        Ccb_halfedge_const_circulator hec = iccb;
        Ccb_halfedge_const_circulator curr = hec;
        do {
            Halfedge_const_handle he = curr;
            if (!he->is_fictitious()) {
                if (seen.find(he->twin()) == seen.end()) {
                    seen.push_front(he);
                    // twin does not already exists
                    *oi++ = CGAL::make_object(he->curve());
                }
            }
            ++curr;
        } while(curr != hec);

        return oi;
    }

    template < class OutputIterator >
    OutputIterator operator()(Face_const_handle fh, 
                              OutputIterator oi,
                              bool just_outer_ccb = false) const {
        
        // value-type is CGAL::Object
        
        for (Ccb_const_iterator occb = fh->outer_ccbs_begin();
             occb != fh->outer_ccbs_end(); occb++) {
            this->operator()(*occb, oi);
        }
        if (!just_outer_ccb) {
            for (Ccb_const_iterator iccb = fh->inner_ccbs_begin();
                 iccb != fh->inner_ccbs_end(); iccb++) {
                this->operator()(*iccb, oi);
            }
            for (Vertex_const_handle vh = fh->isolated_vertices_begin();
                 vh != fh->isolated_vertices_end(); vh++) {
                CGAL_assertion(!vh->is_fictitious());
                *oi++ = CGAL::make_object(vh->point());
            }
        }
        
        return oi;
    }

#if 0
    // Question: Write further operators?
    template < class OutputIterator >
    OutputIterator operator()(OutputIterator oi) const {
        
        // value-type is CGAL::Object
        
        
        return oi;
    }
#endif
    
    //!@}
};

} //namespace CGAL

#endif // CGAL_ARR_GEOMETRY_OF_DCEL_RECORD_2_H

