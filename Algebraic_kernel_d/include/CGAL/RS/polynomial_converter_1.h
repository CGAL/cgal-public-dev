// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.

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

#ifndef CGAL_RS_POLYNOMIAL_CONVERTER_1_H
#define CGAL_RS_POLYNOMIAL_CONVERTER_1_H

namespace CGAL{
namespace RS_AK1{

template <class InputPolynomial_,class OutputPolynomial_>
struct Polynomial_converter_1:
public std::unary_function<InputPolynomial_,OutputPolynomial_>{
        typedef InputPolynomial_                        InpPolynomial_1;
        typedef OutputPolynomial_                       OutPolynomial_1;
        OutPolynomial_1 operator()(const InpPolynomial_1&)const;
}; // class Polynomial_converter_1

template <>
Polynomial<Gmpz>
Polynomial_converter_1<Polynomial<Gmpq>,Polynomial<Gmpz> >::operator()(
                const Polynomial<Gmpq> &p)const{
        std::vector<Gmpz> outcoeffs;
        unsigned degree=p.degree();
        mpz_t lcm;
        mpz_init(lcm);
        mpz_lcm(lcm,mpq_denref(p[0].mpq()),mpq_denref(p[degree].mpq()));
        for(unsigned i=1;i<degree;++i)
                mpz_lcm(lcm,lcm,mpq_denref(p[i].mpq()));
        for(unsigned i=0;i<=degree;++i){
                Gmpz c;
                mpz_divexact(c.mpz(),lcm,mpq_denref(p[i].mpq()));
                mpz_mul(c.mpz(),c.mpz(),mpq_numref(p[i].mpq()));
                outcoeffs.push_back(c);
        }
        mpz_clear(lcm);
        return Polynomial<Gmpz>(outcoeffs.begin(),outcoeffs.end());
}

template <>
Polynomial<Gmpz>
Polynomial_converter_1<Polynomial<Gmpfr>,Polynomial<Gmpz> >::operator()(
                const Polynomial<Gmpfr> &p)const{
        std::vector<Gmpz> outcoeffs;
        unsigned degree=p.degree();
        mpfr_exp_t min_e=mpfr_zero_p(p[degree].fr())?
                          0:
                          mpfr_get_exp(p[degree].fr());
        Gmpfr::Precision_type max_prec=p[degree].get_precision();
        for(unsigned i=0;i<degree;++i){
                mpfr_exp_t current_e=mpfr_zero_p(p[i].fr())?
                                     0:
                                     mpfr_get_exp(p[i].fr());
                if(min_e>current_e)
                        min_e=current_e;
                Gmpfr::Precision_type current_prec=p[i].get_precision();
                if(max_prec<current_prec)
                        max_prec=current_prec;
        }
        // We have now the maximum and minimum exponents of the
        // coefficients.
        if(min_e<0){
                // Multiply the coefficients by 2^(max_prec-min_e-1). The
                // precision must be added because we can have a mantissa
                // of the form .000...01. Other possibility consists in
                // checking the mantissa and adding the number of trailing
                // zeroes to min_e.
                unsigned long mulby2=max_prec-min_e-1;
                for(unsigned i=0;i<=degree;++i){
                        Gmpz cz;
                        Gmpfr czfr(p[i]);
                        CGAL_assertion_code(int round=)
                        mpfr_mul_2ui(czfr.fr(),p[i].fr(),mulby2,GMP_RNDD);
                        CGAL_assertion(!round);
                        mpfr_get_z(cz.mpz(),czfr.fr(),GMP_RNDD);
                        CGAL_assertion(!mpfr_cmp_z(czfr.fr(),cz.mpz()));
                        outcoeffs.push_back(cz);
                }
        }else{
                // Otherwise, we have only integer coefficients. We may
                // reduce the size of the integer coefficients when min_e>0
                // by dividing by 2^(min_e), but this would incur in
                // another copy of the coefficient.
                for(unsigned i=0;i<=degree;++i){
                        Gmpz cz;
                        CGAL_assertion_code(int round=)
                        mpfr_get_z(cz.mpz(),p[i].fr(),GMP_RNDD);
                        CGAL_assertion(!round);
                        outcoeffs.push_back(cz);
                }
        }
        return Polynomial<Gmpz>(outcoeffs.begin(),outcoeffs.end());
}

template <>
Polynomial<Gmpz>
Polynomial_converter_1<Polynomial<Gmpfr>,Polynomial<Gmpz> >::operator()(
                const Polynomial<Gmpfr> &p)const{
        std::vector<Gmpz> outcoeffs;
        unsigned degree=p.degree();
        mpfr_exp_t min_e=mpfr_zero_p(p[degree].fr())?
                          0:
                          mpfr_get_exp(p[degree].fr());
        Gmpfr::Precision_type max_prec=p[degree].get_precision();
        for(unsigned i=0;i<degree;++i){
                mpfr_exp_t current_e=mpfr_zero_p(p[i].fr())?
                                     0:
                                     mpfr_get_exp(p[i].fr());
                if(min_e>current_e)
                        min_e=current_e;
                Gmpfr::Precision_type current_prec=p[i].get_precision();
                if(max_prec<current_prec)
                        max_prec=current_prec;
        }
        // We have now the maximum and minimum exponents of the
        // coefficients.
        if(min_e<0){
                // Multiply the coefficients by 2^(max_prec-min_e-1). The
                // precision must be added because we can have a mantissa
                // of the form .000...01. Other possibility consists in
                // checking the mantissa and adding the number of trailing
                // zeroes to min_e.
                unsigned long mulby2=max_prec-min_e-1;
                for(unsigned i=0;i<=degree;++i){
                        Gmpz cz;
                        Gmpfr czfr(p[i]);
                        CGAL_assertion_code(int round=)
                        mpfr_mul_2ui(czfr.fr(),p[i].fr(),mulby2,GMP_RNDD);
                        CGAL_assertion(!round);
                        mpfr_get_z(cz.mpz(),czfr.fr(),GMP_RNDD);
                        CGAL_assertion(!mpfr_cmp_z(czfr.fr(),cz.mpz()));
                        outcoeffs.push_back(cz);
                }
        }else{
                // Otherwise, we have only integer coefficients. We may
                // reduce the size of the integer coefficients when min_e>0
                // by dividing by 2^(min_e), but this would incur in
                // another copy of the coefficient.
                for(unsigned i=0;i<=degree;++i){
                        Gmpz cz;
                        CGAL_assertion_code(int round=)
                        mpfr_get_z(cz.mpz(),p[i].fr(),GMP_RNDD);
                        CGAL_assertion(!round);
                        outcoeffs.push_back(cz);
                }
        }
        return Polynomial<Gmpz>(outcoeffs.begin(),outcoeffs.end());
}

} // namespace RS_AK1
} // namespace CGAL

#endif // CGAL_RS_POLYNOMIAL_CONVERTER_1_H
