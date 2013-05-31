// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRcomANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : QdX
// File          : demos/xtriangulate/zy_skeletonizer.C
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.6 $
// Revision_date : $Date: 2009-09-19 01:39:16 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

//! startup:    xtriangulate -sx=100 -sy=100 -btm=-1.3 -top=2 -below=-2.3 -above=2.3 11

// #define NDEBUG 1

//! \file zy_skeletonizer.C
//! produces skeletons from surface bodies (undead factory)

// #define CGAL_DEBUG_GL_LISTS 1

#include "includes_common.h"
#include "include/surface_analysis_interface.h"
#include "include/skeletonizer_interface.h"
#include "include/Surface_skeletonizer_impl.h"

CGAL::Timer tm_external_timer;

// Kernel_2 kernel_2_inst;

// NOTE NOTE: make sure both surf analysis and skeletonizer use the same
// kernel_2 instance: i.e. rely on get_static_instance() or smth..

static XSurface_analysis surf_engine_interface;

static CGAL::Surface_skeletonizer_impl< Kernel_2, XSurface_analysis >
    skeletonizer_impl(0);
//     &kernel_2_inst);

XSkeletonizer::XSkeletonizer() {
    skeletonizer_impl.set_surface_engine(&surf_engine_interface);
}

void XSkeletonizer::reset_parameters(unsigned _n_slices_x,
           unsigned _n_slices_y, double _en_left, double _en_right,
            double _en_btm, double _en_top,
            double _z_below, double _z_above) {
    
    skeletonizer_impl.reset_parameters(_n_slices_x, _n_slices_y, _en_left,
           _en_right, _en_btm, _en_top, _z_below, _z_above);
}

bool XSkeletonizer::triangulation_valid() {
    return skeletonizer_impl.triangulation_valid();
}

bool XSkeletonizer::load_surface(const char *filename) {

    return skeletonizer_impl.load_surface(filename);
}

bool XSkeletonizer::overlay_surface(const char *filename) {
    return skeletonizer_impl.overlay_surface(filename);
}

const Triangles_vector& XSkeletonizer::get_triangles() {
    return skeletonizer_impl.get_triangles();
}

const std::vector< unsigned >& XSkeletonizer::get_zlift_ranges() {
    return skeletonizer_impl.get_zlift_ranges();
}

const Triangles_vector& XSkeletonizer::get_sil_triangles() {
    return skeletonizer_impl.get_sil_triangles();
}

const Point_vec_3d& XSkeletonizer::get_vertices() {
    return skeletonizer_impl.get_vertices();
}


void XSkeletonizer::compute_normals(const Point_vec_3d& in, 
        Point_vec_3d& normals, Point_3d& centroid, bool is_inversed) {

    skeletonizer_impl.compute_normals(in, normals, centroid, is_inversed);
}

// too many 'skeletonize' ? ))
void XSkeletonizer::skeletonize(bool use_abs_bounds) {
    skeletonizer_impl.skeletonize(use_abs_bounds);
}

