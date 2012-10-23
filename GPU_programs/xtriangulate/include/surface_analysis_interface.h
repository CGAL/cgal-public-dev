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

#ifndef XTRIANGULATE_SURF_ANALYSIS_H
#define XTRIANGULATE_SURF_ANALYSIS_H

#define CGAL_BISOLVE_ENABLE_ARCAVOID 1 // default TODO?
#define CGAL_ACK_BITSTREAM_USES_E08_TREE 1 // do not change
#define CGAL_BISOLVE_USE_ADJUSTABLE_PRECISION 1 // 1 is default

#define CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER 1
#define CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE 1

#define CGAL_BISOLVE_DEBUG 0
#define CGAL_BISOLVE_VERBOSE 0

#define CGAL_BISOLVE_USE_RS_AK 0
#define CGAL_BISOLVE_USE_RS_ISOLATOR 1 // 1 is default

#define CGAL_ACK_DEBUG_FLAG 0
#define CGAL_ACK_DEBUG_PRINT std::cout

#define CGAL_BISOLVE_USE_RESULTANT_COFACTORS 1
#define CGAL_BISOLVE_ARRANGEMENTS 1 // 0 is default

#define CGAL_BISOLVE_USE_GMP 1
#define CGAL_BISOLVE_USE_CORE 0

//#define CGAL_USE_M_K_DESCARTES 1
// #define CGAL_AK3_USE_NEW_ISOLATORS 1

#define CGAL_AK_USE_OLD_BITSTREAM_DESCARTES 0

#define CGAL_MODULAR_FILTER_OFF

// #define CGAL_WEBXACLI_MULTITHREADED

#define STILL_ALIVE std::cout << __FILE__ << "; " << __LINE__ << "\n";

#include "includes_common.h"

#include <CGAL/Arithmetic_kernel.h>

#define CGAL_WEBXALCI_DEFMUTEX(mutex) \
    extern pthread_mutex_t mutex;

#define CGAL_WEBXALCI_MUTEXLOCK(mutex) do { \
    printf("try lock: %x\n", pthread_self()); \
    pthread_mutex_lock(&mutex); \
    printf("lock acquired: %x\n", pthread_self()); \
    } while(0);

#define CGAL_WEBXALCI_MUTEXUNLOCK(mutex) do { \
    printf("unlock: %x\n", pthread_self()); \
    pthread_mutex_unlock(&mutex); \
    } while(0);

#if CGAL_BISOLVE_USE_GMP
#include <CGAL/GMP_arithmetic_kernel.h>
typedef CGAL::GMP_arithmetic_kernel AT;
#endif

#if CGAL_BISOLVE_USE_CORE
#include <CGAL/CORE_arithmetic_kernel.h>
typedef CGAL::CORE_arithmetic_kernel AT;
#endif

#define CGAL_XTri_USE_TIMERS

/** **************************************************************************/

#define CGAL_BISOLVE_USE_BIGCD 1
#define CGAL_BIGCD_USE_SHIFT 0
#define CGAL_BIGCD_CHECK_SANITY 1

#define CGAL_BISOLVE_USE_GPU_RESULTANTS 1 // default?
#define CGAL_BISOLVE_CHECK_GPU_RESULTANTS_SANITY 0 // default 0

#define CGAL_BISOLVE_USE_GPU_GCDS 1  // default?
#define CGAL_BISOLVE_CHECK_GPU_GCDS_SANITY 0 // default 1

#define CGAL_BISOLVE_USE_NTL  1 // default 1 ??

#include <CGAL/symbolic_standalone.h>

/** **************************************************************************/

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Algebraic_kernel_2/Rounding_ak_d_1.h>

#include <CGAL/Algebraic_kernel_d/Generic_isolator.h>

#if CGAL_BISOLVE_USE_RS_AK
#include <CGAL/Algebraic_kernel_d/Float_traits.h>
#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>
#else
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d_1_generator.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#endif

#include <CGAL/Algebraic_kernel_d_2.h>
#include <CGAL/Arrangement_2l/Adjacencies_3.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_enums.h>

#include <CGAL/Algebraic_kernel_1/Arcavoid_root_isolator.h>

typedef AT::Integer Integer;
typedef AT::Rational Rational;

#if CGAL_BISOLVE_USE_RS_AK
  typedef Algebraic_kernel_rs_gmpz_d_1 Algebraic_kernel_d_1;
#else

#if CGAL_BISOLVE_USE_RS_ISOLATOR
#warning using RS root isolator !

   typedef CGAL::Algebraic_kernel_d_1_generator< Integer, Rational >
    ::Algebraic_kernel_with_qir_and_rs_1 Internal_ak_1;

//         typedef CGAL::Algebraic_kernel_d_1
//     < Integer, Rational,
//       CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi
//            < Integer, Rational >,
//       CGAL::internal::Real_solve< CGAL::Polynomial< Integer >, Rational >
//     > Algebraic_kernel_with_qir_and_rs_1;

    typedef CGAL::internal::Rounding_ak_d_1< Internal_ak_1 > AK_1;

#else
#warning using bitstream root isolator !
// TODO use Algebraic_kernel_d_1_generator?
    typedef CGAL::Algebraic_kernel_d_1
    < Integer, Rational,
      CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi
           < Integer, Rational >,
      CGAL::internal::Bitstream_descartes
           < CGAL::internal::Bitstream_coefficient_kernel< Integer > >
    > AK_1;
#endif
#endif // !CGAL_BISOLVE_USE_RS_AK

typedef CGAL::Algebraic_curve_kernel_2< AK_1 > Kernel_2;

typedef CGAL::Polynomial< Integer > Polynomial_1;
typedef CGAL::Polynomial< Polynomial_1 > Polynomial_2;
typedef CGAL::Polynomial< Polynomial_2 > Polynomial_3;

typedef Kernel_2::Curve_analysis_2 Curve_analysis_2;
typedef Curve_analysis_2::Status_line_1 Status_line_1;
typedef Kernel_2::Algebraic_real_1 X_coordinate_1;
typedef Kernel_2::Algebraic_real_2 Xy_coordinate_2;

typedef std::vector< double > Double_vector;
typedef std::vector< X_coordinate_1 > X_coord_vector;

typedef CGAL::Adjacencies_3::Adjacency_vector Adjacency_vector;
typedef CGAL::Adjacencies_3::Adjacency_pair Adj_pair;

typedef std::pair < CGAL::Dcel_feature, const void * > Dcel_feature_data;
typedef std::vector < Dcel_feature_data > Dcel_feature_vec;

class XSurface_analysis {

public:
    XSurface_analysis() {
    }

    bool surface_valid();

    void clear();

    //! loads and analyses the surface
    bool load_surface(const char *filename);

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

#endif // XTRIANGULATE_SURF_ANALYSIS_H
