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
// Library       : QdX
// File          : demos/xsurface/include/arrangements.h
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef XSURFACE_ARRANGEMENTS_H
#define XSURFACE_ARRANGEMENTS_H

/*!\file arrangements.h 
 * defines interface to arrangments computation
 */

#include "includes_common.h"

class XSurface_arrangements {

public:

    XSurface_arrangements() : base_index(-1), valid_flag(false) {
    }

    void compute_spatial_coords(Point_vec_3f& in, 
        Point_vec_3f& verts, Point_vec_3f *normals) const;

    //! which = 0: outer circle, which = 1: tube circle
    void compute_aux_coords(Point_vec_3f& in, Point_vec_3f& verts, 
            int which) const;

    //! some value which relates to base surface dimenstions 
    double get_base_surf_extent() const;
    
    //! origin ?
    void get_base_surf_origin(Point_3d& origin) const;

    void get_pole_coords(Point_3f& pt) const;

    bool read_base_surfaces(const char *filename);
    bool read_surface_set(const char *filename, unsigned which, 
            bool clear_flag);

    //! approximate spatial points and arcs of an arrangement  
    void render(const CGAL::Bbox_2& bbox, int res_w, int res_h); 

    void compute_arr_on_surface(); //! computes arrangement

    void compute_overlay_arr();    //! overlays two arrangements

    void clear_approximations() {
        xarcs_approx.clear();
        points_approx.clear();
    }

    const std::vector< Point_vec_3d >& get_arcs_approx() const
    { return xarcs_approx; }

    const Point_vec_3d& get_points_approx() const 
    { return points_approx; }

    const std::vector< unsigned >& get_color_index() const 
    { return xarcs_color_index; }

    const std::vector< unsigned >& get_points_color_index() const 
    { return points_color_index; }

    bool arrangement_valid() const 
    { return valid_flag; }    

    bool base_surface_valid() const;
    
protected:

    //! index of the base surface to compute arrangements on
    unsigned base_index; 
    bool valid_flag; //! indicates whether arrangement is ready
   
    std::vector< Point_vec_3d > xarcs_approx;
    std::vector< unsigned > xarcs_color_index;
    std::vector< unsigned > points_color_index;
    
    Point_vec_3d points_approx;
};

#endif // XSURFACE_ARRANGEMENTS_H
