
// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

// This file defines several instances of Algebraic_curve_kernel_2 
// specializations that are used in tests and demos. 

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_2_GENERATOR_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_2_GENERATOR_H 1

#include <CGAL/basic.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1_generator.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h>
#include <CGAL/Filtered_algebraic_curve_kernel_2.h>


namespace CGAL {

template<typename Coefficient, 
         typename Bound = typename Get_arithmetic_kernel< Coefficient >::Arithmetic_kernel::Rational>
struct Algebraic_curve_kernel_2_generator {

    // Unfiltered kernels

    typedef CGAL::Algebraic_curve_kernel_2 
    < typename CGAL::Algebraic_kernel_d_1_generator<Coefficient,Bound>
        ::Default_algebraic_kernel_1 >
        Default_algebraic_curve_kernel_2;

    typedef CGAL::Algebraic_curve_kernel_2 
    < typename CGAL::Algebraic_kernel_d_1_generator<Coefficient,Bound>
        ::Algebraic_kernel_with_bisection_and_descartes_1 >
        Algebraic_curve_kernel_with_bisection_and_descartes_2;

    typedef CGAL::Algebraic_curve_kernel_2 
    < typename CGAL::Algebraic_kernel_d_1_generator<Coefficient,Bound>
        ::Algebraic_kernel_with_qir_and_bitstream_1 >
         Algebraic_curve_kernel_with_qir_and_bitstream_2;

     typedef CGAL::Algebraic_curve_kernel_2 
    < typename CGAL::Algebraic_kernel_d_1_generator<Coefficient,Bound>
        ::Algebraic_kernel_with_qir_and_descartes_1 >
         Algebraic_curve_kernel_with_qir_and_descartes_2;

    // Filtered kernels

    typedef CGAL::Filtered_algebraic_curve_kernel_2 
    < typename CGAL::Algebraic_kernel_d_1_generator<Coefficient,Bound>
        ::Default_algebraic_kernel_1 >
        Default_filtered_algebraic_curve_kernel_2;

    typedef CGAL::Filtered_algebraic_curve_kernel_2 
    < typename CGAL::Algebraic_kernel_d_1_generator<Coefficient,Bound>
        ::Algebraic_kernel_with_bisection_and_descartes_1 >
        Filtered_algebraic_curve_kernel_with_bisection_and_descartes_2;

    typedef CGAL::Filtered_algebraic_curve_kernel_2 
    < typename CGAL::Algebraic_kernel_d_1_generator<Coefficient,Bound>
        ::Algebraic_kernel_with_qir_and_bitstream_1 >
         Filtered_algebraic_curve_kernel_with_qir_and_bitstream_2;

};

} //namespace CGAL


#endif
