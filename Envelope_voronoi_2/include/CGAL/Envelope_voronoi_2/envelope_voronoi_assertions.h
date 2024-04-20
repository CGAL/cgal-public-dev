// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Scripts/developer_scripts/create_assertions.sh $
// $Id: create_assertions.sh 44139 2008-07-14 11:16:41Z spion $
// 
//
// Author(s)     : Geert-Jan Giezeman, Sven Sch�nherr

// Generated from script create_assertions.sh

// macro definitions
// =================
// assertions
// ----------


#if defined(CGAL_ENVELOPE_VORONOI_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_assertion(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_assertion_code(CODE)
#else
#  define CGAL_envelope_voronoi_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_assertion_code(CODE) CODE
#  define CGAL_envelope_voronoi_assertions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_ASSERTIONS

#if defined(CGAL_ENVELOPE_VORONOI_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_exactness_assertion_code(CODE)
#else
#  define CGAL_envelope_voronoi_exactness_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_exactness_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_exactness_assertion_code(CODE) CODE
#  define CGAL_envelope_voronoi_exactness_assertions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_ASSERTIONS

#if defined(CGAL_ENVELOPE_VORONOI_NO_ASSERTIONS) \
  || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_expensive_assertion(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_assertion_code(CODE)
#else
#  define CGAL_envelope_voronoi_expensive_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_expensive_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_expensive_assertion_code(CODE) CODE
#  define CGAL_envelope_voronoi_expensive_assertions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_ASSERTIONS

#if defined(CGAL_ENVELOPE_VORONOI_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_expensive_exactness_assertion(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_exactness_assertion_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_exactness_assertion_code(CODE)
#else
#  define CGAL_envelope_voronoi_expensive_exactness_assertion(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_expensive_exactness_assertion_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::assertion_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_expensive_exactness_assertion_code(CODE) CODE
#  define CGAL_envelope_voronoi_expensive_exactness_assertions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_ASSERTIONS


// preconditions
// -------------

#if defined(CGAL_ENVELOPE_VORONOI_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_precondition(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_precondition_code(CODE)
#else
#  define CGAL_envelope_voronoi_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_precondition_code(CODE) CODE
#  define CGAL_envelope_voronoi_preconditions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_PRECONDITIONS

#if defined(CGAL_ENVELOPE_VORONOI_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_exactness_precondition_code(CODE)
#else
#  define CGAL_envelope_voronoi_exactness_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_exactness_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_exactness_precondition_code(CODE) CODE
#  define CGAL_envelope_voronoi_exactness_preconditions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_PRECONDITIONS

#if defined(CGAL_ENVELOPE_VORONOI_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_expensive_precondition(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_precondition_code(CODE)
#else
#  define CGAL_envelope_voronoi_expensive_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_expensive_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_expensive_precondition_code(CODE) CODE
#  define CGAL_envelope_voronoi_expensive_preconditions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_PRECONDITIONS

#if defined(CGAL_ENVELOPE_VORONOI_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_expensive_exactness_precondition(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_exactness_precondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_exactness_precondition_code(CODE)
#else
#  define CGAL_envelope_voronoi_expensive_exactness_precondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_expensive_exactness_precondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::precondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_expensive_exactness_precondition_code(CODE) CODE
#  define CGAL_envelope_voronoi_expensive_exactness_preconditions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_PRECONDITIONS


// postconditions
// --------------

#if defined(CGAL_ENVELOPE_VORONOI_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_postcondition(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_postcondition_code(CODE)
#else
#  define CGAL_envelope_voronoi_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_postcondition_code(CODE) CODE
#  define CGAL_envelope_voronoi_postconditions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_POSTCONDITIONS

#if defined(CGAL_ENVELOPE_VORONOI_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_exactness_postcondition_code(CODE)
#else
#  define CGAL_envelope_voronoi_exactness_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_exactness_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_exactness_postcondition_code(CODE) CODE
#  define CGAL_envelope_voronoi_exactness_postconditions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_POSTCONDITIONS

#if defined(CGAL_ENVELOPE_VORONOI_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_expensive_postcondition(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_postcondition_code(CODE)
#else
#  define CGAL_envelope_voronoi_expensive_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_expensive_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_expensive_postcondition_code(CODE) CODE
#  define CGAL_envelope_voronoi_expensive_postconditions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_POSTCONDITIONS

#if defined(CGAL_ENVELOPE_VORONOI_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_expensive_exactness_postcondition(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_exactness_postcondition_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_exactness_postcondition_code(CODE)
#else
#  define CGAL_envelope_voronoi_expensive_exactness_postcondition(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_expensive_exactness_postcondition_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::postcondition_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_expensive_exactness_postcondition_code(CODE) CODE
#  define CGAL_envelope_voronoi_expensive_exactness_postconditions 1
#endif // CGAL_ENVELOPE_VORONOI_NO_POSTCONDITIONS


// warnings
// --------

#if defined(CGAL_ENVELOPE_VORONOI_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_warning(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_warning_code(CODE)
#else
#  define CGAL_envelope_voronoi_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_warning_code(CODE) CODE
#  define CGAL_envelope_voronoi_warnings 1
#endif // CGAL_ENVELOPE_VORONOI_NO_WARNINGS

#if defined(CGAL_ENVELOPE_VORONOI_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_exactness_warning_code(CODE)
#else
#  define CGAL_envelope_voronoi_exactness_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_exactness_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_exactness_warning_code(CODE) CODE
#  define CGAL_envelope_voronoi_exactness_warnings 1
#endif // CGAL_ENVELOPE_VORONOI_NO_WARNINGS

#if defined(CGAL_ENVELOPE_VORONOI_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_expensive_warning(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_warning_code(CODE)
#else
#  define CGAL_envelope_voronoi_expensive_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_expensive_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_expensive_warning_code(CODE) CODE
#  define CGAL_envelope_voronoi_expensive_warnings 1
#endif // CGAL_ENVELOPE_VORONOI_NO_WARNINGS

#if defined(CGAL_ENVELOPE_VORONOI_NO_WARNINGS) || defined(CGAL_NO_WARNINGS) \
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXACTNESS) && !defined(CGAL_CHECK_EXACTNESS))\
  || (!defined(CGAL_ENVELOPE_VORONOI_CHECK_EXPENSIVE) && !defined(CGAL_CHECK_EXPENSIVE)) \
  || defined(NDEBUG)
#  define CGAL_envelope_voronoi_expensive_exactness_warning(EX) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_exactness_warning_msg(EX,MSG) (static_cast<void>(0))
#  define CGAL_envelope_voronoi_expensive_exactness_warning_code(CODE)
#else
#  define CGAL_envelope_voronoi_expensive_exactness_warning(EX) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__))
#  define CGAL_envelope_voronoi_expensive_exactness_warning_msg(EX,MSG) \
   (CGAL::possibly(EX)?(static_cast<void>(0)): ::CGAL::warning_fail( # EX , __FILE__, __LINE__, MSG))
#  define CGAL_envelope_voronoi_expensive_exactness_warning_code(CODE) CODE
#  define CGAL_envelope_voronoi_expensive_exactness_warnings 1
#endif // CGAL_ENVELOPE_VORONOI_NO_WARNINGS


