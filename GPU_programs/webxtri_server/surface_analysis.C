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
// File          : demos/webxalci/rasterizer.C
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:04 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#define NDEBUG 1

#include "include/IPC.h"
#include "include/CGAL_includes.h"
#include "include/skeletonizer.h"

#include "include/Surface_analysis_impl.h"
#include "include/SFA_interface.h"

typedef CGAL::Surface_analysis_impl< Kernel_2 > Surface_analysis_3;

typedef std::map< MD5_digest, Surface_analysis_3, MD5_compare >
    Surface_analysis_cache;

static Surface_analysis_cache SFA_cache;
    
static Surface_analysis_3 *g_SFA_inst = 0;

//static CGAL::Surface_analysis_impl< Kernel_2 > surf_engine(0);

bool SFA_interface::analyse_surface(const MD5_digest& surface_ID,
        const Polynomial_3& in) {

    std::pair< Surface_analysis_cache::iterator, bool > res =
        SFA_cache.insert(std::make_pair(surface_ID, Surface_analysis_3()));

    g_SFA_inst = &res.first->second; // set the correct SFA instance
    if(res.second) {// if insertion took place, load surface
        g_SFA_inst->load_surface(in);
        return g_SFA_inst->surface_valid();
    }
    return true;
}

bool SFA_interface::load_surface(const MD5_digest& surface_ID) {

    Surface_analysis_cache::iterator i = SFA_cache.find(surface_ID);
    if(i == SFA_cache.end())
        return false;
    g_SFA_inst = &i->second; // set the current instance in use
    return true;
}

bool SFA_interface::surface_valid() {
    return g_SFA_inst->surface_valid();
}

const Kernel_2& SFA_interface::kernel_2_inst() {
    return g_SFA_inst->kernel_2_inst();
}

Curve_analysis_2 SFA_interface::silhouette_curve() {
    return g_SFA_inst->silhouette_curve();
}

Polynomial_2 SFA_interface::poly_at_rational_x(const Rational& x) {
    return g_SFA_inst->poly_at_rational_x(x);
}

unsigned SFA_interface::locate_pt_at_sil(const X_coordinate_1& x,
        int arcno, Dcel_feature_vec& vec) {
    return g_SFA_inst->locate_pt_at_sil(x, arcno, vec);
}

unsigned SFA_interface::locate_pt_in_face(const X_coordinate_1& x,
        const Rational& y, Dcel_feature_vec& vec) {

    return g_SFA_inst->locate_pt_in_face(x, y, vec);
}

void SFA_interface::approximate_zs_at_sil(const X_coordinate_1& x,
    int arcno, const Dcel_feature_data& dcel, const double& threshold,
                              Double_vector& zs) {

    g_SFA_inst->approximate_zs_at_sil(x, arcno, dcel, threshold, zs);
}


void SFA_interface::approximate_zs_in_face(const Polynomial_2& poly_at_x,
        const Rational& y, const Rational& threshold,
        const Rational& z_below, const Rational& z_above,
         Double_vector& zs, unsigned& n_skips_below,
          unsigned& n_skips_above, bool use_z_bounds) {

    return g_SFA_inst->approximate_zs_in_face(poly_at_x, y, threshold,
        z_below, z_above, zs, n_skips_below, n_skips_above, use_z_bounds);
}

int SFA_interface::univariate_roots(const Polynomial_1& poly, const
        Rational& threshold, X_coord_vector& zs, bool sq_free) {

    return g_SFA_inst->univariate_roots(poly, threshold, zs, sq_free);
}

#if 0
    void exact_zs_in_face(const Xy_coordinate_2& xy,
        unsigned idcel, const Rational& threshold, const Rational& z_below,
        const Rational& z_above, Double_vector& zs) const;
#endif

void SFA_interface::roots_at_y_boundary(const Rational& y_bnd,
       const Rational& threshold, X_coord_vector& zs,
        bool is_sq_free, Polynomial_2 *ppoly) {

    g_SFA_inst->roots_at_y_boundary(y_bnd, threshold, zs,
           is_sq_free, ppoly);
}

bool SFA_interface::adjacency_info(const Dcel_feature_data& dcel1,
             const Dcel_feature_data& dcel2, Adjacency_vector& adj) {

    return g_SFA_inst->adjacency_info(dcel1, dcel2, adj);
}

void SFA_interface::compute_normals(const Point_vec_3d& in,
        Point_vec_3d& normals, bool is_inversed) {

    g_SFA_inst->compute_normals(in, normals, is_inversed);
}

Curve_analysis_2 SFA_interface::z_curve(Rational& z_below,
           Rational& z_above) {

    return g_SFA_inst->z_curve(z_below, z_above);
}

void SFA_interface::compute_y_range(double& y_min, double& y_max,
         bool& is_set) {

    g_SFA_inst->compute_y_range(y_min, y_max, is_set);
}

void SFA_interface::dump_arr() {
    g_SFA_inst->dump_arr();
}

