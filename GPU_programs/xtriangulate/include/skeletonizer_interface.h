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
// File          : demos/xtriangulate/include/skeletonizer.h
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.2 $
// Revision_date : $Date: 2009-07-24 13:21:30 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef XSKELETONIZER_INTERFACE_H
#define XSKELETONIZER_INTERFACE_H

/*!\file skeletonizer.h 
 * produces skeletons from surface bodies (undead factory)
 */

#include "includes_common.h"

class XSkeletonizer {

public:

    XSkeletonizer();

    void reset_parameters(unsigned _n_slices_x = 1, unsigned _n_slices_y = 1,
            double _en_left = 0.1, double _en_right = 0.1,
            double _en_btm = 0.1, double _en_top = 0.1,
            double _z_below = -2.0, double _z_above = 2.0);

    bool triangulation_valid();

    bool surface_valid();    

    bool load_surface(const char *filename);
    bool overlay_surface(const char *filename);

    //bool load_raw_data(const char *);

    void skeletonize(bool use_abs_bounds = false);

    void compute_normals(const Point_vec_3d& in, Point_vec_3d& normals,
        Point_3d& centroid, bool is_inversed = false);

    const std::vector< unsigned >& get_zlift_ranges();
    const Triangles_vector& get_triangles();
    const Triangles_vector& get_sil_triangles();

    const Point_vec_3d& get_vertices();
};

#endif // XSKELETONIZER_H
