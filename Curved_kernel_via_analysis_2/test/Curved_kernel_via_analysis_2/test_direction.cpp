/*! \file test_direction.cpp
 * Some tests for directions of arcs
 */

#include <CGAL/Algebraic_kernel_d/flags.h>

#include <CGAL/basic.h>

#include <iostream>

#ifndef CGAL_USE_CORE

int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return (0);
}

#else

// Allows to use the Filtered_curve_kernel_via_analysis_2
#ifndef CGAL_ACK_USE_FILTERED_CKvA_2
#define CGAL_ACK_USE_FILTERED_CKvA_2 0
#endif

// What is the coefficient type of the input?
#ifndef CGAL_ACK_COEFFICIENT
#define CGAL_ACK_COEFFICIENT CGAL::CORE_arithmetic_kernel::Integer
#define CGAL_ACK_COEFFICIENT_IS_INTEGER 1
#endif

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_curve_kernel_2_generator.h>

#if CGAL_ACK_USE_FILTERED_CKvA_2
#include <CGAL/Filtered_algebraic_curve_kernel_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Filtered_curved_kernel_via_analysis_2_impl.h>
#else
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>
#endif

typedef CGAL_ACK_COEFFICIENT Coefficient;

#if !CGAL_ACK_USE_FILTERED_CKvA_2
typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>
::Algebraic_curve_kernel_with_qir_and_bitstream_2
Algebraic_curve_kernel_2;
#else
typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>
::Filtered_algebraic_curve_kernel_with_qir_and_bitstream_2
Algebraic_curve_kernel_2;
#endif

#if !CGAL_ACK_USE_FILTERED_CKvA_2
typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
Curved_kernel_2; 
#else
typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
Exact_curved_kernel_2; 
typedef CGAL::Filtered_curved_kernel_via_analysis_2<Exact_curved_kernel_2>
Curved_kernel_2; 
#endif

typedef Curved_kernel_2                                 Traits_2; 
typedef Traits_2::Curve_2                               Bezier_curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;

// The main program.
int main (int argc, char **argv)
{

  X_monotone_curve_2 arc;

  Traits_2 traits;
  X_monotone_curve_2 arc2 = traits.construct_opposite_2_object()(arc);
  assert(arc.is_left_to_right());
  assert(arc2.is_left_to_right());
  assert(traits.compare_endpoints_xy_2_object()(arc) == CGAL::SMALLER);
  assert(traits.compare_endpoints_xy_2_object()(arcf) == CGAL::LARGER);
  
  return 0;
}

#endif
