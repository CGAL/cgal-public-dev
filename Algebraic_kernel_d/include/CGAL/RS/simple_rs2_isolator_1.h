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

#ifndef CGAL_RS_SIMPLE_RS2_ISOLATOR_1_H
#define CGAL_RS_SIMPLE_RS2_ISOLATOR_1_H

#include <CGAL/RS/simple_rs2_calls.h>
#include <CGAL/Gmpfi.h>
#include <vector>

namespace CGAL{
namespace RS2{

template <class Polynomial_,class Bound_>
class Simple_rs2_isolator_1{
        public:
        typedef Polynomial_                             Polynomial;
        typedef Bound_                                  Bound;
        private:
        typedef Gmpfi                                   Interval;
        public:
        Simple_rs2_isolator_1(const Polynomial&);
        Polynomial polynomial()const;
        int number_of_real_roots()const;
        bool is_exact_root(int i)const;
        Bound left_bound(int i)const;
        Bound right_bound(int i)const;
        private:
        Polynomial _polynomial;
        std::vector<Interval> _real_roots;
};

template <class Polynomial_,class Bound_>
Simple_rs2_isolator_1<Polynomial_,Bound_>::
Simple_rs2_isolator_1(const Polynomial_ &p){
        CGAL_error_msg("not implemented for these polynomial/bound types");
};

template <>
Simple_rs2_isolator_1<CGAL::Polynomial<CGAL::Gmpz>,Gmpfr>::
Simple_rs2_isolator_1(const CGAL::Polynomial<CGAL::Gmpz> &p):_polynomial(p){
        int numsols;
        unsigned int degree=p.degree();
        mpz_t *coeffs=(mpz_t*)malloc((degree+1)*sizeof(mpz_t));
        mpfi_ptr *intervals_mpfi=(mpfi_ptr*)malloc(degree*sizeof(mpfi_ptr));
        for(int i=0;i<=degree;++i)
                coeffs[i][0]=*(p[i].mpz());
        RS2::Simple_rs2_calls::init_solver();
        RS2::Simple_rs2_calls::create_rs_upoly(coeffs,
                                               degree,
                                               rs_get_default_up());
        free(coeffs);
        set_rs_precisol(0);
        set_rs_verbose(0);
        rs_run_algo((char*)"UISOLE");
        numsols=RS2::Simple_rs2_calls::affiche_sols_eqs(intervals_mpfi);
        for(int j=0;j<numsols;++j)
                _real_roots.push_back(Gmpfi(intervals_mpfi[j]));
        free(intervals_mpfi);
};

template <class Polynomial_,class Bound_>
Polynomial_
Simple_rs2_isolator_1<Polynomial_,Bound_>::polynomial()const{
        return _polynomial;
};

template <class Polynomial_,class Bound_>
int
Simple_rs2_isolator_1<Polynomial_,Bound_>::number_of_real_roots()const{
        return _real_roots.size();
};

template <class Polynomial_,class Bound_>
bool
Simple_rs2_isolator_1<Polynomial_,Bound_>::is_exact_root(int i)const{
        return _real_roots[i].inf()==_real_roots[i].sup();
};

template <class Polynomial_,class Bound_>
Bound_
Simple_rs2_isolator_1<Polynomial_,Bound_>::left_bound(int i)const{
        return _real_roots[i].inf();
};

template <class Polynomial_,class Bound_>
Bound_
Simple_rs2_isolator_1<Polynomial_,Bound_>::right_bound(int i)const{
        return _real_roots[i].sup();
};

} // namespace RS2
} // namespace CGAL

#endif // CGAL_RS_SIMPLE_RS2_ISOLATOR_1_H
