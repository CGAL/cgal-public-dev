// Copyright (c) 2011 INRIA Nancy-Grand Est (France).
// Copyright (c) 2011 National and Kapodistrian University of Athens (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Authors: Yacine Bouzidi <yacine.bouzidi@inria.fr>
//          Luis Peñaranda <luis.penaranda@gmx.com>
//          Marc Pouget <marc.pouget@inria.fr>
//          Fabrice Rouillier <fabrice.rouillier@inria.fr>

#ifndef CGAL_RS_ALGEBRAIC_KERNEL_RS_2
#define CGAL_RS_ALGEBRAIC_KERNEL_RS_2

// TODO: replace the old AK1 by the new one
#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/RS/rs_calls_1.h>
#include <CGAL/RS/functors_2.h>

namespace CGAL{

// _Polynomial: the type of bivariate polynomials forming the systems that
// the kernel will be able to solve.
// _Bound: the type of the numbers representing the bounds of the boxes of
// the solutions of the system.
// _AK1: the univariate algebraic kernel that this bivariate kernel
// refines. It defaults to the RS univariate algebraic kernel.
// _Ptraits: the polynomial traits containing the operations on
// polynomials. It defaults to the CGAL polynomial traits.
template <class _Polynomial,
          class _Bound,
          class _AK1=CGAL::Algebraic_kernel_rs_gmpz_d_1,
          class _Ptraits=CGAL::Polynomial_traits_d<_Polynomial> >
struct Algebraic_kernel_rs_2:
public _AK1{

        typedef _Polynomial                             Polynomial_2;
        typedef _Bound                                  Bound;
        typedef _Ptraits                                Ptraits;
        typedef _AK1                                    AK1;
        typedef typename AK1::Polynomial_1              Polynomial_1;
        typedef int                                     size_type;
        typedef typename AK1::Algebraic_real_1          Algebraic_real_1;
        typedef CGAL::RS3::Algebraic_2<Polynomial_2>    Algebraic_real_2;
        typedef RS3::Multiplicity                       Multiplicity_type;
        typedef typename Ptraits::Innermost_coefficient_type
                                                        Coefficient;

        // TODO: write new functions in RS3 namespace to init/reset
        Algebraic_kernel_rs_2(){
                CGAL::init_solver();
        };

        ~Algebraic_kernel_rs_2(){
                CGAL::reset_solver();
        };

        typedef RS3::Construct_alg_2<Polynomial_2,
                                     Bound,
                                     Coefficient,
                                     Algebraic_real_1>
                                                Construct_algebraic_real_2;
        typedef RS3::Compute_polynomial_x_2<Algebraic_real_2,Polynomial_1>
                                                Compute_polynomial_x_2;
        typedef RS3::Compute_polynomial_y_2<Algebraic_real_2,Polynomial_1>
                                                Compute_polynomial_y_2;
        typedef RS3::Isolate_2<Polynomial_2>    Isolate_2;
        typedef RS3::Isolate_x_2<Polynomial_2>  Isolate_x_2;
        typedef RS3::Isolate_y_2<Polynomial_2>  Isolate_y_2;
        typedef RS3::Is_square_free_2<Ptraits>  Is_square_free_2;
        typedef Ptraits::Make_square_free       Make_square_free_2;
        typedef Ptraits::Square_free_factorize_up_to_constant_factor
                                                Square_free_factorize_2;
        typedef RS3::Is_coprime_2<Ptraits>      Is_coprime_2;
        typedef RS3::Make_coprime_2<Ptraits>    Make_coprime_2;
        typedef RS3::Solve_2<Polynomial_2,Bound>
                                                Solve_2;
        /*
        typedef RS3::Number_of_solutions_2<Polynomial_2>
                                                Number_of_solutions_2;
        typedef RS3::Sign_at_2<Polynomial_2>    Sign_at_2;
        typedef RS3::Compare_x_2                Compare_x_2;
        typedef RS3::Compare_y_2                Compare_y_2;
        typedef RS3::Compare_xy_2               Compare_xy_2;
        typedef RS3::Bound_between_x_2          Bound_between_x_2;
        typedef RS3::Bound_between_y_2          Bound_between_y_2;
        typedef RS3::Approximate_absolute_x_2
                                                Approximate_absolute_x_2;
        typedef RS3::Approximate_absolute_y_2
                                                Approximate_absolute_y_2;
        typedef RS3::Approximate_relative_x_2
                                                Approximate_relative_x_2;
        typedef RS3::Approximate_relative_y_2
                                                Approximate_relative_y_2;
        */

#define CGALRS_CREATE_FUNCTION_OBJECT(T,N) \
        T N##_object()const{return T();}
        CGALRS_CREATE_FUNCTION_OBJECT(Construct_algebraic_real_2,
                                      construct_algebraic_real_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Compute_polynomial_x_2,
                                      compute_polynomial_x_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Compute_polynomial_y_2,
                                      compute_polynomial_y_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Isolate_2,isolate_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Isolate_x_2,isolate_x_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Isolate_y_2,isolate_y_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Is_square_free_2,is_square_free_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Make_square_free_2,make_square_free_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Square_free_factorize_2,
                                      square_free_factorize_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Is_coprime_2,is_coprime_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Make_coprime_2,make_coprime_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Solve_2,solve_2)
        /*
        CGALRS_CREATE_FUNCTION_OBJECT(Number_of_solutions_2,
                                      number_of_solutions_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Sign_at_2,sign_at_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Compare_x_2,compare_x_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Compare_y_2,compare_y_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Compare_xy_2,compare_xy_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Bound_between_x_2,bound_between_x_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Bound_between_y_2,bound_between_y_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Approximate_absolute_x_2,
                                      approximate_absolute_x_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Approximate_absolute_y_2,
                                      approximate_absolute_y_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Approximate_relative_x_2,
                                      approximate_relative_x_2)
        CGALRS_CREATE_FUNCTION_OBJECT(Approximate_relative_y_2,
                                      approximate_relative_y_2)
        */
#undef CGALRS_CREATE_FUNCTION_OBJECT
};  // Algebraic_kernel_rs_2

} // namespace CGAL

#endif  // CGAL_RS_ALGEBRAIC_KERNEL_RS_2
