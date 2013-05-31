// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/webxalci/include/rasterizer.h
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:11 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*! \file rasterizer.h
 *  \brief Defines class \c Rasterizer 
 *  
 *  Algebraic curve analysis and rasterization routines for \c Xalci_server
 */

#ifndef XALCI_RASTERIZER_H
#define XALCI_RASTERIZER_H

// maximal number of arcs, vertices and faces to print respectively
#define MAX_ARCS_PRINT_OUT          1000
#define MAX_VERTICES_PRINT_OUT      500
#define MAX_FACES_PRINT_OUT         500

#define MAX_COORDS_PRINTOUT 120*1024 // maximal size of coordinates print-out
//! \class Rasterizer
//! provides the server with an interface to EXACUS curve analysis and
//! rendering
class Rasterizer 
{
public:
    //!\name public methods
    //!@{    

    //! default constructor
    Rasterizer(); 
    
    
    //! sets up rasterizer parameters
    void setup(int width_, int height_);
            
    //! analyses an algebraic curve and returns a set of arcs in
    //! \c arcs parameter
    void analyse_curve(Analysis_entry& analysis,
       const Poly_int_vector& poly_vec);
       
    void draw_topology_graph(const Analysis_entry& analysis,
        const CGAL::Bbox_2& box, void *dst);
    
    //! renders a set of curve arcs onto the bitmap \c plot, with
    //! \c box defining drawing window and \c indices - a set of arc
    //! indices, \c mode = 0: default rasterizer,
    //! \c mode = 1: rasterizer complete curve in one-color
    void plot_arcs(const Analysis_entry& analysis, uint *indices,
            uint n_indices, const CGAL::Bbox_2& box, Rasterize_mode mode,
            void *dst);

    void plot_points(const Analysis_entry& analysis, uint *indices,
            uint n_indices, const CGAL::Bbox_2& box, Rasterize_mode mode,
            QPixmap *plot);

    void plot_faces(const Analysis_entry& analysis, uint *indices,
        uint n_indices, const CGAL::Bbox_2& box, Rasterize_mode mode,
        QPixmap *plot);

    //! plots a curve equation using 1D space subdivision
    void plot_subdivision(const Poly_int_vector& poly_vec,
            const CGAL::Bbox_2& box, QPixmap *plot);
            
    //! draw coordinate axis
    void draw_axis(const CGAL::Bbox_2& box, QPixmap *plot);

    //! prints out a set of curve arcs to \c out
    void print_out(Analysis_entry& analysis);

    void print_endpoint(const Arc_2& arc, CGAL::Arr_curve_end end,
        std::ostream& os);

    //! issues a point location query on given point
    bool locate_point(const Analysis_entry& analysis,
        const CGAL::Bbox_2& box, const Point_coords& coords,
            SHM_Point_query_reply *reply);
    
    //!@}
public:
    //!\name private data members
    //!@{

    //! dimensions of the drawing window
    static int width;
    static int height;
    
    //!@}
};

#endif // XALCI_RASTERIZER_H

