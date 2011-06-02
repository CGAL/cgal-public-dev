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

#ifndef CGAL_RS_SIMPLE_SIGNAT_1_H
#define CGAL_RS_SIMPLE_SIGNAT_1_H

#include <CGAL/Gmpfi.h>
//#include <boost/mpl/assert.hpp>

namespace CGAL{
namespace SimpleAK1{

template <class Polynomial_,class Bound_/*,class Ptraits_*/>
struct Simple_signat_1{
        typedef Polynomial_                                     Polynomial;
        typedef Bound_                                          Bound;
        //typedef Ptraits_                                        Ptraits;
        typedef Polynomial_traits_d<Polynomial>                 Ptraits;
        Polynomial pol;
        Simple_signat_1(const Polynomial &p):pol(p){};
        CGAL::Sign operator()(const Bound&);
}; //struct Simple_signat_1{

template <class Polynomial_,class Bound_/*,class Ptraits_*/>
inline CGAL::Sign
Simple_signat_1<Polynomial_,Bound_>::operator()(const Bound_ &x){
        typedef Polynomial_                                     Polynomial;
        typedef Bound_                                          Bound;
        //typedef Ptraits_                                        Ptraits;
        typedef Polynomial_traits_d<Polynomial>                 Ptraits;
        typedef typename Ptraits::Degree                        Degree;
        typedef Real_embeddable_traits<Bound>                   REtraits;
        typedef typename REtraits::Sgn                          BSign;
        typedef Algebraic_structure_traits<Bound>               AStraits;
        // This generic signat works only when Bound_ is an exact type. For
        // non-exact types, an implementation must be provided.
        //BOOST_MPL_ASSERT((boost::is_same<AStraits::Is_exact,Tag_true>));
        int d=Degree()(pol);
        Bound h(pol[d]);
        for(int i=1;i<=d;++i)
                h=h*x+pol[d-i];
        return BSign()(h);
};

template <>
inline CGAL::Sign
Simple_signat_1<Polynomial<Gmpz>,Gmpfr>::operator()(const Gmpfr &x){
        typedef Ptraits::Degree                                 Degree;
        typedef Ptraits::Substitute                             Substitute;
        typedef Simple_signat_1<Polynomial,Gmpq>                Exact_sign;
        int d=Degree()(pol);
        if(d==0)
                return pol[0].sign();
        Gmpfi h(pol[d],x.get_precision()+2*d);
        if(h.sign()!=Uncertain<Sign>::indeterminate())
                return Exact_sign(pol)(x);
        for(int i=1;i<=d;++i){
                h*=x;
                if(h.sign()!=Uncertain<Sign>::indeterminate())
                        return Exact_sign(pol)(x);
                h+=pol[d-i];
                if(h.sign()!=Uncertain<Sign>::indeterminate())
                        return Exact_sign(pol)(x);
        }
        CGAL_assertion(h.sign()!=Uncertain<Sign>::indeterminate());
        return h.sign();
};

} // namespace SimpleAK1
} // namespace CGAL

#endif // CGAL_RS_SIMPLE_SIGNAT_1_H
