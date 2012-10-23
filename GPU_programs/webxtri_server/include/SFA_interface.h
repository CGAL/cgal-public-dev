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
// File          : demos/xtriangulate/include/surf_analysis.h
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.5 $
// Revision_date : $Date: 2009-07-24 13:21:30 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef XTRI_SFA_INTERFACE_H
#define XTRI_SFA_INTERFACE_H


//! minimal interface to surface analysis required by Skeletonizer 
class SFA_interface {

public:
    SFA_interface() {
    }

    bool surface_valid();

    void clear();

    bool analyse_surface(const MD5_digest& surface_ID,
        const Polynomial_3& poly);
        
    bool load_surface(const MD5_digest& surface_ID);

    //! returns an instance of ACK_2
    const Kernel_2& kernel_2_inst();

    //! returns silhouette curve object
    Curve_analysis_2 silhouette_curve();

    Polynomial_2 poly_at_rational_x(const Rational& x);

    unsigned locate_pt_at_sil(const X_coordinate_1& x,
        int arcno, Dcel_feature_vec& vec);

    unsigned locate_pt_in_face(const X_coordinate_1& x,
        const Rational& y, Dcel_feature_vec& vec);

    //! computes approximations of z-coordinates at a given point lying on
    //! silhouette, returns an arrangement feature assigned to this point
    void approximate_zs_at_sil(const X_coordinate_1& x, int arcno,
       const Dcel_feature_data& dcel, const double& threshold,
                  Double_vector& zs);
    
    //! computes approximations of z-coordinates of rational point within face 
    //! returns an arrangement feature (face handle) and # of roots skipped
    //! on bottom boundary \c z_below
    void approximate_zs_in_face(const Polynomial_2& poly_at_x,
        const Rational& y, const Rational& threshold,
        const Rational& z_below, const Rational& z_above,
         Double_vector& zs, unsigned& n_skips_below,
          unsigned& n_skips_above, bool use_z_bounds = true);
    
    int univariate_roots(const Polynomial_1& poly, const Rational& threshold,
             X_coord_vector& zs, bool sq_free = true);

#if 0
    void exact_zs_in_face(const Xy_coordinate_2& xy,
        unsigned idcel, const Rational& threshold, const Rational& z_below, 
        const Rational& z_above, Double_vector& zs) const;
#endif    

    //! computes approximation of real roots of biv polynomial at rational y
    //! if \c ppoly == NULL uses silhouette curve
    void roots_at_y_boundary(const Rational& y_bnd, 
       const Rational& threshold, X_coord_vector& zs,
        bool is_sq_free, Polynomial_2 *ppoly = NULL);

    //! returns adjacency info of two points given by index coordinates \c i1
    //! and \c i2 , returns \c true if two points belong to features of the
    //! same time
    bool adjacency_info(const Dcel_feature_data& dcel1,
               const Dcel_feature_data& dcel2, Adjacency_vector& adj);

    void compute_normals(const Point_vec_3d& in, Point_vec_3d& normals,
        bool is_inversed = false);

    Curve_analysis_2 z_curve(Rational& z_below, Rational& z_above);

    void compute_y_range(double& y_min, double& y_max, bool& is_set);

// DEBUG only
    void dump_arr();
};

#endif // XTRI_SFA_INTERFACE_H
