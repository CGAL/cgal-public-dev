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
// Library       : ??
// File          : demos/webxalci/include/CGAL_includes.h
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.14 $
// Revision_date : $Date: 2007/02/06 12:03:11 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*! \file CGAL_includes.h
 *  \brief Includes CGAL code
 *  
 *  Algebraic curve analysis and rasterization routines for \c xalci_server
 */

#ifndef XTRI_CGAL_INCLUDES_H
#define XTRI_CGAL_INCLUDES_H

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

// #define CGAL_AK3_USE_NEW_ISOLATORS 1

#define CGAL_AK_USE_OLD_BITSTREAM_DESCARTES 0
#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1

#define CGAL_BISOLVE_USE_GMP 1
#define CGAL_BISOLVE_USE_CORE 0

#define CGAL_MODULAR_FILTER_OFF

// #define CGAL_WEBXACLI_MULTITHREADED

#define STILL_ALIVE std::cout << __FILE__ << "; " << __LINE__ << "\n";

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

/** **************************************************************************/

#define CGAL_BISOLVE_USE_BIGCD 1
#define CGAL_BIGCD_USE_SHIFT 0
#define CGAL_BIGCD_CHECK_SANITY 0

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

#if CGAL_BISOLVE_USE_RS_AK
#include <CGAL/Algebraic_kernel_d/Float_traits.h>
#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>
#else
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d_1_generator.h>
#endif

#include <CGAL/Algebraic_kernel_d_2.h>
#include <CGAL/Arrangement_2l/Adjacencies_3.h>
#include <CGAL/Arrangement_2l/Restricted_cad_3_enums.h>

#include <CGAL/Arcavoid_root_isolator.h>


typedef AT::Integer Integer;
typedef AT::Rational Rational;

#if CGAL_BISOLVE_USE_RS_AK
  typedef Algebraic_kernel_rs_gmpz_d_1 Algebraic_kernel_d_1;
#else

#if CGAL_BISOLVE_USE_RS_ISOLATOR
#warning using RS root isolator !
// #include <include/Unreal_solve.h>

   typedef CGAL::Algebraic_kernel_d_1_generator< Integer, Rational >
    ::Algebraic_kernel_with_qir_and_rs_1 Internal_ak_1;

    typedef CGAL::internal::Rounding_ak_d_1< Internal_ak_1 > AK_1;

#else
#warning using bitstream root isolator !
    typedef CGAL::Algebraic_kernel_d_1
    < Integer, Rational,
      CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi
           < Integer, Rational >,
      CGAL::internal::Bitstream_descartes
        < CGAL::internal::Bitstream_descartes_rndl_tree_traits
            < CGAL::internal::Bitstream_coefficient_kernel< Integer > >
        >
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

#include <boost/array.hpp>

template <class NT, std::size_t N>
std::ostream& operator <<(std::ostream& os,
        const boost::array<NT, N>& v) {

    CGAL::output_range(os, v.begin(), v.end(), ", ", "[", "]");
    return os;
}

typedef boost::array< float, 3 >  Point_3f;
typedef boost::array< double, 3 > Point_3d;
typedef boost::array< double, 2 > Point_2d;

typedef std::vector< Point_3f > Point_vec_3f;
typedef std::vector< Point_3d > Point_vec_3d;

typedef boost::array< unsigned, 3 > Tri_index;
typedef std::vector< Tri_index > Triangles_vector;

//! surface analysis cache entry
// struct Analysis_entry {
//     Analysis_entry() { }
//     
//     MD5_digest surface_ID;   // unique surface identifier
//     CGAL_Arrangement_2 arr;  // stores complete arrangement
// 
// };

#endif // XTRI_CGAL_INCLUDES_H

