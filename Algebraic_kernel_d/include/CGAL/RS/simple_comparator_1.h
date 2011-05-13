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

#ifndef CGAL_RS_SIMPLE_COMPARATOR_1_H
#define CGAL_RS_SIMPLE_COMPARATOR_1_H

#include <CGAL/Polynomial_traits_d.h>

namespace CGAL{
namespace SimpleAK1{

template <class Polynomial_,
          class Bound_,
          class Refiner_,
          class Signat_/*,
          class Ptraits_*/>
struct Simple_comparator_1{
        typedef Polynomial_                                     Polynomial;
        typedef Bound_                                          Bound;
        typedef Refiner_                                        Refiner;
        typedef Signat_                                         Signat;
        //typedef Ptraits_                                        Ptraits;
        typedef CGAL::Polynomial_traits_d<Polynomial>           Ptraits;
        typedef typename Ptraits::Gcd_up_to_constant_factor     Gcd;
        typedef typename Ptraits::Degree                        Degree;

        CGAL::Comparison_result
        operator()(const Polynomial &p1,Bound &l1,Bound &r1,
                   const Polynomial &p2,Bound &l2,Bound &r2)const{
                if(l1<=l2){
                        if(r1<l2)
                                return SMALLER;
                }else{
                        if(r2<l1)
                                return LARGER;
                }
                Polynomial G=Gcd()(p1,p2);
                if(Degree()(G)==1)
                        return compare_unequal(p1,l1,r1,p2,l2,r2);
                Signat sg(G);
                CGAL::Sign sleft=sg(l1>l2?l1:l2);
                if(sleft==ZERO)
                        return EQUAL;
                CGAL::Sign sright=sg(r1<r2?r1:r2);
                if(sleft!=sright)
                        return EQUAL;
                else
                        return compare_unequal(p1,l1,r1,p2,l2,r2);
        }

        // This function compares two algebraic numbers, assuming that they
        // are not equal.
        CGAL::Comparison_result
        compare_unequal(const Polynomial &p1,Bound &l1,Bound &r1,
                        const Polynomial &p2,Bound &l2,Bound &r2)const{
                do{
                        Refiner()(p1,l1,r1,100);
                        Refiner()(p2,l2,r2,100);
                }while(l1<=l2?r1>=l2:r2>=l1);
                return (r1<l2?SMALLER:LARGER);
        }

}; // struct Simple_comparator_1

} // namespace SimpleAK1
} // namespace CGAL

#endif // CGAL_RS_SIMPLE_COMPARATOR_1_H
