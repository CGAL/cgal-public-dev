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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_SIMPLE_AK_1_H
#define CGAL_RS_SIMPLE_AK_1_H

#include <cstddef> // included only to define size_t
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/RS/simple_algebraic_1.h>
#include <CGAL/RS/simple_comparator_1.h>
#include <CGAL/RS/simple_signat_1.h>

namespace CGAL{
namespace SimpleAK1{

template <class Polynomial_,
          class Bound_,
          class Isolator_,
          class Refiner_,
          class Ptraits_=CGAL::Polynomial_traits_d<Polynomial_> >
class Simple_algebraic_kernel_1{
        public:
        typedef Polynomial_                             Polynomial_1;
        typedef typename Polynomial_1::NT               Coefficient;
        typedef Bound_                                  Bound;
        private:
        typedef Isolator_                               Isolator;
        typedef Refiner_                                Refiner;
        typedef Ptraits_                                Ptraits;
        typedef CGAL::SimpleAK1::Simple_signat_1<Polynomial_1,
                                                 Bound>
                                                        Signat;
        typedef CGAL::SimpleAK1::Simple_comparator_1<Polynomial_1,
                                                     Bound,
                                                     Refiner,
                                                     Signat,
                                                     Ptraits>
                                                        Comparator;
        public:
        typedef CGAL::SimpleAK1::Simple_algebraic_1<Polynomial_1,
                                                    Bound,
                                                    Refiner,
                                                    Comparator,
                                                    Ptraits>
                                                        Algebraic_real_1;
        typedef size_t                                  size_type;
        typedef unsigned                                Multiplicity_type;

        // default constructor and destructor
        public:
        Simple_algebraic_kernel_1(){};
        ~Simple_algebraic_kernel_1(){};

        // functors from the CGAL concept
        public:
        // TODO: Construct_algebraic_real_1
        typedef CGAL::SimpleAK1::Compute_polynomial_1<Polynomial_1,
                                                      Algebraic_real_1>
                                                        Compute_polynomial_1;
        typedef CGAL::SimpleAK1::Isolate_1<Polynomial_1,
                                           Bound,
                                           Algebraic_real_1,
                                           Isolator,
                                           Comparator,
                                           Signat,
                                           Ptraits>
                                                        Isolate_1;
        typedef typename Ptraits::Is_square_free        Is_square_free_1;
        typedef typename Ptraits::Make_square_free      Make_square_free_1;
        typedef typename Ptraits::Square_free_factorize Square_free_factorize_1;
        typedef CGAL::SimpleAK1::Is_coprime_1<Polynomial_1,Ptraits>
                                                        Is_coprime_1;
        typedef CGAL::SimpleAK1::Make_coprime_1<Polynomial_1,Ptraits>
                                                        Make_coprime_1;
        typedef CGAL::SimpleAK1::Solve_1<Polynomial_1,
                                         Bound,
                                         Algebraic_real_1,
                                         Isolator,
                                         Signat,
                                         Ptraits>
                                                        Solve_1;
        typedef CGAL::SimpleAK1::Number_of_solutions_1<Polynomial_1,Isolator>
                                                        Number_of_solutions_1;

        // TODO: Sign_at_1
        typedef CGAL::SimpleAK1::Compare_1<Algebraic_real_1,
                                           Bound,
                                           Comparator>
                                                        Compare_1;
        // TODO: Bound_between_1
        // TODO: Approximate_absolute_1
        // TODO: Approximate_relative_1

#define CREATE_SIMPLE_FUNCTION_OBJECT(T,N) \
        T N##_object()const{return T();}
        /*
        CREATE_SIMPLE_FUNCTION_OBJECT(Construct_algebraic_real_1,
                                      construct_algebraic_real_1)
        */
        CREATE_SIMPLE_FUNCTION_OBJECT(Compute_polynomial_1,
                                      compute_polynomial_1)
        CREATE_SIMPLE_FUNCTION_OBJECT(Isolate_1,
                                      isolate_1)
        CREATE_SIMPLE_FUNCTION_OBJECT(Is_square_free_1,
                                      is_square_free_1)
        CREATE_SIMPLE_FUNCTION_OBJECT(Make_square_free_1,
                                      make_square_free_1)
        CREATE_SIMPLE_FUNCTION_OBJECT(Square_free_factorize_1,
                                      square_free_factorize_1)
        CREATE_SIMPLE_FUNCTION_OBJECT(Is_coprime_1,
                                      is_coprime_1)
        CREATE_SIMPLE_FUNCTION_OBJECT(Make_coprime_1,
                                      make_coprime_1)
        CREATE_SIMPLE_FUNCTION_OBJECT(Solve_1,
                                      solve_1)
        CREATE_SIMPLE_FUNCTION_OBJECT(Number_of_solutions_1,
                                      number_of_solutions_1)
        /*
        CREATE_SIMPLE_FUNCTION_OBJECT(Sign_at_1,
                                      sign_at_1)
        */
        CREATE_SIMPLE_FUNCTION_OBJECT(Compare_1,
                                      compare_1)
        /*
        CREATE_SIMPLE_FUNCTION_OBJECT(Bound_between_1,
                                      bound_between_1)
        CREATE_SIMPLE_FUNCTION_OBJECT(Approximate_absolute_1,
                                      approximate_absolute_1)
        CREATE_SIMPLE_FUNCTION_OBJECT(Approximate_relative_1,
                                      approximate_relative_1)
        */
#undef CREATE_SIMPLE_FUNCTION_OBJECT

}; // class Simple_algebraic_kernel_1

} // namespace SimpleAK1
} // namespace CGAL

#endif // CGAL_RS_SIMPLE_AK_1_H
