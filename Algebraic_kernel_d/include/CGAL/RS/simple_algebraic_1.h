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

#ifndef CGAL_RS_SIMPLE_ALGEBRAIC_1_H
#define CGAL_RS_SIMPLE_ALGEBRAIC_1_H

#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Polynomial_traits_d.h>

namespace CGAL{
namespace SimpleAK1{

// This class represents the simplest algebraic number one can think about.
// One algebraic number is represented by the polynomial of which it is
// root and the two endpoints of an interval that contains it, and no other
// root. Polynomial type and bound type are the first two template
// parameters.
//
// The third template parameter is a refiner, a function object that
// receives the polynomial and the bounds defining an algebraic number and
// an integer p, and modifies the two bounds until the difference between
// the two bounds is less than x*2^(-p), where x is the value of the
// represented algebraic number. The signature of a refiner must be:
//      void
//      Refiner_()(const Polynomial_&,Bound_&,Bound_&,int p);
//
// The fourth template argument is a comparator, a function object that
// receives the polynomials and bounds defining two algebraic numbres and
// just compares them, returning a CGAL::Comparison_result. The signature
// of a comparator must be:
//      CGAL::Comparison_result
//      Comparator_()(const Polynomial_&,Bound_&,Bound_&,
//                    const Polynomial_&,Bound_&,Bound_&);
// The comparator can modify the bounds, with the condition that the
// algebraic numbers remain consistent (one and only one root on each
// interval).

template <class Polynomial_,
          class Bound_,
          class Refiner_,
          class Comparator_/*,
          class Ptraits_*/>
class Simple_algebraic_1{
        private:
        typedef Polynomial_                             Polynomial;
        typedef Bound_                                  Bound;
        typedef Refiner_                                Refiner;
        typedef Comparator_                             Comparator;
        //typedef Ptraits_                                Ptraits;
        typedef Polynomial_traits_d<Polynomial>         Ptraits;
        typedef typename Ptraits::Coefficient_type      Coefficient;
        typedef typename Ptraits::Scale                 Scale;
        typedef Simple_algebraic_1<Polynomial,Bound,Refiner,Comparator>
                                                        Algebraic;
        Polynomial pol;
        mutable Bound left,right;
        public:
        Simple_algebraic_1(){};
        Simple_algebraic_1(const Polynomial &p,
                           const Bound &l,
                           const Bound &r):pol(p),left(l),right(r){
                CGAL_assertion(l<=r);
        }
        Simple_algebraic_1(const Algebraic &a):
        pol(a.pol),left(a.left),right(a.right){}
        // TODO: constructors from types such as int, unsigned and double
        ~Simple_algebraic_1(){}

        Simple_algebraic_1& operator=(const Algebraic &a){
                pol=a.get_pol();
                left=a.get_left();
                right=a.get_right();
        }

        Polynomial get_pol()const{return pol;}
        Bound& get_left()const{return left;}
        Bound& get_right()const{return right;}

        Algebraic& operator-()const{
                return Algebraic(Scale(get_pol(),Coefficient(-1)),
                                 -right,
                                 -left);
        }

#define CGAL_RS_COMPARE_THIS_WITH(_a) \
        (Comparator()(get_pol(),get_left(),get_right(), \
                      (_a).get_pol(),(_a).get_left(),(_a).get_right()))
        bool operator==(const Algebraic &a)const
                {return CGAL_RS_COMPARE_THIS_WITH(a)==CGAL::EQUAL;}
        bool operator!=(const Algebraic &a)const
                {return CGAL_RS_COMPARE_THIS_WITH(a)!=CGAL::EQUAL;}
        bool operator<(const Algebraic &a)const
                {return CGAL_RS_COMPARE_THIS_WITH(a)==CGAL::SMALLER;}
        bool operator<=(const Algebraic &a)const
                {return CGAL_RS_COMPARE_THIS_WITH(a)!=CGAL::LARGER;}
        bool operator>(const Algebraic &a)const
                {return CGAL_RS_COMPARE_THIS_WITH(a)==CGAL::LARGER;}
        bool operator>=(const Algebraic &a)const
                {return CGAL_RS_COMPARE_THIS_WITH(a)!=CGAL::SMALLER;}
#undef CGAL_RS_COMPARE_THIS_WITH

#ifdef IEEE_DBL_MANT_DIG
#define CGAL_RS_DBL_PREC        IEEE_DBL_MANT_DIG
#else
#define CGAL_RS_DBL_PREC        53
#endif
        double to_double(){
                typedef Real_embeddable_traits<Bound>                   RT;
                typedef typename RT::To_double                          TD;
                Refiner()(get_pol(),get_left(),get_right(),CGAL_RS_DBL_PREC);
                CGAL_assertion(TD()(get_left())==TD()(get_right()));
                return TD()(get_left());
        }
        std::pair<double,double> to_interval(){
                typedef Real_embeddable_traits<Bound>                   RT;
                typedef typename RT::To_interval                        TI;
                return std::make_pair(TI()(get_left().first),
                                      TI()(get_right().second));
        }
#undef CGAL_RS_DBL_PREC

}; // class Simple_algebraic_1

} // namespace SimpleAK1

// We define Simple_algebraic_1 as real-embeddable
template <class Polynomial_,class Bound_,class Refiner_,class Comparator_>
class Real_embeddable_traits<SimpleAK1::Simple_algebraic_1<Polynomial_,
                                                           Bound_,
                                                           Refiner_,
                                                           Comparator_> >{
        typedef Polynomial_                             P;
        typedef Bound_                                  B;
        typedef Refiner_                                R;
        typedef Comparator_                             C;

        public:

        typedef SimpleAK1::Simple_algebraic_1<P,B,R,C>  Type;
        typedef CGAL::Tag_true                          Is_real_embeddable;
        typedef bool                                    Boolean;
        typedef CGAL::Sign                              Sign;
        typedef CGAL::Comparison_result                 Comparison_result;

        typedef INTERN_RET::Real_embeddable_traits_base<Type,CGAL::Tag_true>
                                                        Base;
        typedef typename Base::Compare                  Compare;

        class Sgn:public std::unary_function<Type,CGAL::Sign>{
                public:
                CGAL::Sign operator()(const Type &a)const{
                        return Compare()(a,Type(0));
                }
        };

        class To_double:public std::unary_function<Type,double>{
                public:
                double operator()(const Type &a)const{return a.to_double();}
        };

        class To_interval:
        public std::unary_function<Type,std::pair<double,double> >{
                public:
                std::pair<double,double> operator()(const Type &a)const{
                        return a.to_interval();
                }
        };

        class Is_zero:public std::unary_function<Type,Boolean>{
                public:
                bool operator()(const Type &a)const{
                        return Sgn()(a)==CGAL::ZERO;
                }
        };

        class Is_finite:public std::unary_function<Type,Boolean>{
                public:
                bool operator()(const Type&)const{return true;}
        };

        class Abs:public std::unary_function<Type,Type>{
                public:
                Type operator()(const Type &a)const{
                        return Sgn()(a)==CGAL::NEGATIVE?-a:a;
                }
        };
};

} // namespace CGAL

#endif // CGAL_RS_SIMPLE_ALGEBRAIC_1_H
