// Copyright (c) 1999,2000,2001  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Menelaos Karavelas <mkaravel@cse.nd.edu>

// This file is automatically generated by
// scripts/filtered_predicates_generator.pl

// MK: January 19, 2004
// This file was originally automatically generated by
// scripts/filtered_predicates_generator.pl
// Modifications have been made on top of it.

#ifndef CGAL_ARITHMETIC_FILTER_APOLLONIUS_GRAPH_FTC2_H
#define CGAL_ARITHMETIC_FILTER_APOLLONIUS_GRAPH_FTC2_H

#include <CGAL/Profile_counter.h>

CGAL_BEGIN_NAMESPACE

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/* inline */
bool
ad_is_hidden_test_ring_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA ad_is_hidden_test_alg_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_is_hidden_test_ring_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter failures("IA ad_is_hidden_test_alg_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_is_hidden_test_ring_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact());
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/* inline */
bool
ad_is_hidden_test_sqrtf_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA ad_is_hidden_test_sqrtf_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_is_hidden_test_sqrtf_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter failures("IA ad_is_hidden_test_sqrtf_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_is_hidden_test_sqrtf_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact());
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/* CGAL_MEDIUM_INLINE */
Comparison_result
compare_ad_distances_test_ring_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA compare_ad_distances_test_ring_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return compare_ad_distances_test_ring_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x.interval(),
		y.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter
      failures("IA compare_ad_distances_test_ring_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return compare_ad_distances_test_ring_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x.exact(),
		y.exact());
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/* CGAL_MEDIUM_INLINE */
Comparison_result
compare_ad_distances_test_sqrtf_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      calls("IA compare_ad_distances_test_sqrtf_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return compare_ad_distances_test_sqrtf_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x.interval(),
		y.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      failures("IA compare_ad_distances_test_sqrtf_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return compare_ad_distances_test_sqrtf_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x.exact(),
		y.exact());
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
Sign
ad_incircle_test_sqrtf_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA ad_incircle_test_sqrtf_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_incircle_test_sqrtf_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter failures("IA ad_incircle_test_sqrtf_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_incircle_test_sqrtf_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact());
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
Sign
ad_incircle_test_ring_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA ad_incircle_test_ring_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_incircle_test_ring_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter failures("IA ad_incircle_test_ring_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_incircle_test_ring_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact());
  }
}


template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
Sign
ad_incircle_test_sqrtf_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA ad_incircle_test_sqrtf_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_incircle_test_sqrtf_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x3.interval(),
		y3.interval(),
		w3.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter failures("IA ad_incircle_test_sqrtf_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_incircle_test_sqrtf_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x3.exact(),
		y3.exact(),
		w3.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact());
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
Sign
ad_incircle_test_ring_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA ad_incircle_test_ring_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_incircle_test_ring_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x3.interval(),
		y3.interval(),
		w3.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter failures("IA ad_incircle_test_ring_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_incircle_test_ring_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x3.exact(),
		y3.exact(),
		w3.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact());
  }
}


template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
bool
ad_finite_edge_test_sqrtf_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw,
    bool b)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA ad_finite_edge_test_sqrtf_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_finite_edge_test_sqrtf_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x3.interval(),
		y3.interval(),
		w3.interval(),
		x4.interval(),
		y4.interval(),
		w4.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval(),
		b);
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      failures("IA ad_finite_edge_test_sqrtf_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_finite_edge_test_sqrtf_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x3.exact(),
		y3.exact(),
		w3.exact(),
		x4.exact(),
		y4.exact(),
		w4.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact(),
		b);
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
bool
ad_finite_edge_test_ring_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw,
    bool b)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA ad_finite_edge_test_ring_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_finite_edge_test_ring_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x3.interval(),
		y3.interval(),
		w3.interval(),
		x4.interval(),
		y4.interval(),
		w4.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval(),
		b);
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter failures("IA ad_finite_edge_test_ring_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_finite_edge_test_ring_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x3.exact(),
		y3.exact(),
		w3.exact(),
		x4.exact(),
		y4.exact(),
		w4.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact(),
		b);
  }
}


template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
bool
ad_finite_edge_test_degenerated_sqrtf_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw,
    bool b)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      calls("IA ad_finite_edge_test_degenerated_sqrtf_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_finite_edge_test_degenerated_sqrtf_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval(),
		b);
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      failures("IA ad_finite_edge_test_degenerated_sqrtf_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_finite_edge_test_degenerated_sqrtf_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact(),
		b);
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
bool
ad_finite_edge_test_degenerated_ring_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw,
    bool b)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      calls("IA ad_finite_edge_test_degenerated_ring_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_finite_edge_test_degenerated_ring_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval(),
		b);
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      failures("IA ad_finite_edge_test_degenerated_ring_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_finite_edge_test_degenerated_ring_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact(),
		b);
  }
}


template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
bool
ad_finite_edge_test_degenerated_sqrtf_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw,
    bool b)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      calls("IA ad_finite_edge_test_degenerated_sqrtf_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_finite_edge_test_degenerated_sqrtf_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x3.interval(),
		y3.interval(),
		w3.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval(),
		b);
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      failures("IA ad_finite_edge_test_degenerated_sqrtf_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_finite_edge_test_degenerated_sqrtf_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x3.exact(),
		y3.exact(),
		w3.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact(),
		b);
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
bool
ad_finite_edge_test_degenerated_ring_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw,
    bool b)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      calls("IA ad_finite_edge_test_degenerated_ring_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_finite_edge_test_degenerated_ring_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x3.interval(),
		y3.interval(),
		w3.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval(),
		b);
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      failures("IA ad_finite_edge_test_degenerated_ring_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_finite_edge_test_degenerated_ring_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x3.exact(),
		y3.exact(),
		w3.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact(),
		b);
  }
}


template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
bool
ad_infinite_edge_test_sqrtf_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw,
    bool b)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA ad_infinite_edge_test_sqrtf_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_infinite_edge_test_sqrtf_C2(
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x3.interval(),
		y3.interval(),
		w3.interval(),
		x4.interval(),
		y4.interval(),
		w4.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval(),
		b);
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      failures("IA ad_infinite_edge_test_sqrtf_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_infinite_edge_test_sqrtf_C2(
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x3.exact(),
		y3.exact(),
		w3.exact(),
		x4.exact(),
		y4.exact(),
		w4.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact(),
		b);
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
bool
ad_infinite_edge_test_ring_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &qw,
    bool b)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter calls("IA ad_infinite_edge_test_ring_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_infinite_edge_test_ring_C2(
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x3.interval(),
		y3.interval(),
		w3.interval(),
		x4.interval(),
		y4.interval(),
		w4.interval(),
		qx.interval(),
		qy.interval(),
		qw.interval(),
		b);
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      failures("IA ad_infinite_edge_test_ring_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_infinite_edge_test_ring_C2(
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x3.exact(),
		y3.exact(),
		w3.exact(),
		x4.exact(),
		y4.exact(),
		w4.exact(),
		qx.exact(),
		qy.exact(),
		qw.exact(),
		b);
  }
}


template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
bool
ad_is_degenerate_edge_test_sqrtf_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w4)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      calls("IA ad_is_degenerate_edge_test_sqrtf_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_is_degenerate_edge_test_sqrtf_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x3.interval(),
		y3.interval(),
		w3.interval(),
		x4.interval(),
		y4.interval(),
		w4.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      failures("IA ad_is_degenerate_edge_test_sqrtf_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_is_degenerate_edge_test_sqrtf_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x3.exact(),
		y3.exact(),
		w3.exact(),
		x4.exact(),
		y4.exact(),
		w4.exact());
  }
}

template < class CGAL_IA_CT, class CGAL_IA_ET, bool CGAL_IA_PROTECTED,
           class CGAL_IA_CACHE >
/*  */
bool
ad_is_degenerate_edge_test_ring_C2(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w1,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic,
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w2,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w3,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &x4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &y4,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, 
    CGAL_IA_PROTECTED, CGAL_IA_CACHE> &w4)
{
  try
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      calls("IA ad_is_degenerate_edge_test_ring_C2 calls");
    ++calls;
#endif
    Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
    return ad_is_degenerate_edge_test_ring_C2(
		x1.interval(),
		y1.interval(),
		w1.interval(),
		x2.interval(),
		y2.interval(),
		w2.interval(),
		x3.interval(),
		y3.interval(),
		w3.interval(),
		x4.interval(),
		y4.interval(),
		w4.interval());
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
#ifdef CGAL_PROFILE
    static Profile_counter 
      failures("IA ad_is_degenerate_edge_test_ring_C2 failures");
    ++failures;
#endif
    Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
    return ad_is_degenerate_edge_test_ring_C2(
		x1.exact(),
		y1.exact(),
		w1.exact(),
		x2.exact(),
		y2.exact(),
		w2.exact(),
		x3.exact(),
		y3.exact(),
		w3.exact(),
		x4.exact(),
		y4.exact(),
		w4.exact());
  }
}


CGAL_END_NAMESPACE

#endif // CGAL_ARITHMETIC_FILTER_APOLLONIUS_GRAPH_FTC2_H
