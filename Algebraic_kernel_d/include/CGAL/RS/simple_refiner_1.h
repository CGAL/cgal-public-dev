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

#ifndef CGAL_RS_SIMPLE_REFINER_1_H
#define CGAL_RS_SIMPLE_REFINER_1_H

#include <CGAL/Polynomial_traits_d.h>
#include <rs_exports.h> // only for rs_init_rs()
#include <rs3_fncts.h>

namespace CGAL{
namespace RS3{

template <class Polynomial_,class Bound_>
struct Simple_refiner{
        void operator()(const Polynomial_&,Bound_&,Bound_&,int);
}; // class Simple_refiner

template <class Polynomial_,class Bound_>
void
Simple_refiner<Polynomial_,Bound_>::
operator()(const Polynomial_&,Bound_&,Bound_&,int){
        CGAL_error_msg("RS3 refiner not implemented for these types");
        return;
};

template<>
void
Simple_refiner<Polynomial<Gmpz>,Gmpfr>::
operator()(const Polynomial<Gmpz> &pol,Gmpfr &left,Gmpfr &right,int prec){
        typedef Polynomial<Gmpz>                        Polynomial;
        typedef Polynomial_traits_d<Polynomial>         Ptraits;
        typedef Ptraits::Degree                         Degree;
        CGAL_precondition(left<=right);
        int deg=Degree()(pol);
        mpz_t* coefficients=(mpz_t*)malloc((deg+1)*sizeof(mpz_t));
        __mpfi_struct interval;
        if(!left.is_unique()){
                mpfr_t new_left_mpfr;
                mpfr_init2(new_left_mpfr,left.get_precision());
                mpfr_set(new_left_mpfr,left.fr(),GMP_RNDN);
                CGAL_assertion(mpfr_equal_p(new_left_mpfr,left.fr()));
                left=Gmpfr(new_left_mpfr);
        }
        if(!right.is_unique()){
                mpfr_t new_right_mpfr;
                mpfr_init2(new_right_mpfr,right.get_precision());
                mpfr_set(new_right_mpfr,right.fr(),GMP_RNDN);
                CGAL_assertion(mpfr_equal_p(new_right_mpfr,right.fr()));
                right=Gmpfr(new_right_mpfr);
        }
        interval.left=*(left.fr());
        interval.right=*(right.fr());
        for(int i=0;i<=deg;++i)
                coefficients[i][0]=*(pol[i].mpz());
        // TODO: avoid calling this rs_init_rs() many times
        rs_init_rs();
        rs3_refine_u_root(&interval,coefficients,deg,prec,0,0);
        mpfr_custom_init_set(left.fr(),
                             mpfr_custom_get_kind(&interval.left),
                             mpfr_custom_get_exp(&interval.left),
                             mpfr_get_prec(&interval.left),
                             mpfr_custom_get_mantissa(&interval.left));
        mpfr_custom_init_set(right.fr(),
                             mpfr_custom_get_kind(&interval.right),
                             mpfr_custom_get_exp(&interval.right),
                             mpfr_get_prec(&interval.right),
                             mpfr_custom_get_mantissa(&interval.right));
        CGAL_postcondition(left<=right);
        return;
};

} // namespace RS3
} // namespace CGAL

#endif // CGAL_RS_SIMPLE_REFINER_1_H
