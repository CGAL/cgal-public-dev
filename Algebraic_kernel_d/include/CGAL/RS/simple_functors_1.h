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

#ifndef CGAL_RS_SIMPLE_FUNCTORS_1_H
#define CGAL_RS_SIMPLE_FUNCTORS_1_H

#include <CGAL/Polynomial_traits_d.h>
#include <vector>

namespace CGAL{
namespace SimpleAK1{

template <class Polynomial_,
          class Bound_,
          class Algebraic_,
          class Isolator_,
          class Signat_>
struct Solve_1{
        typedef Polynomial_                                     Polynomial_1;
        typedef Bound_                                          Bound;
        typedef Algebraic_                                      Algebraic;
        typedef Isolator_                                       Isolator;
        typedef Signat_                                         Signat;
        typedef Polynomial_traits_d<Polynomial_1>               Ptraits;
        typedef typename Ptraits::Gcd_up_to_constant_factor     Gcd;
        typedef typename Ptraits::Square_free_factorize_up_to_constant_factor
                                                                Sqfr;
        typedef typename Ptraits::Degree                        Degree;
        typedef typename Ptraits::Make_square_free              Sfpart;

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_1 &p,
                                  OutputIterator res)const{
                typedef std::pair<Polynomial_1,int>     polmult;
                typedef std::vector<polmult>            sqvec;

                Polynomial_1 sfp=Sfpart()(p);
                sqvec sfv;
                Sqfr()(p,std::back_inserter(sfv));
                Isolator isol(sfp);
                int *m=(int*)calloc(isol.number_of_real_roots(),sizeof(int));
                for(typename sqvec::iterator i=sfv.begin();i!=sfv.end();++i){
                        int k=Degree()(i->first);
                        Signat signof(i->first);
                        for(int j=0;k&&j<isol.number_of_real_roots();++j){
                                if(!m[j]){
                                        CGAL::Sign sg_l=
                                                signof(isol.left_bound(j));
                                        CGAL::Sign sg_r=
                                                signof(isol.right_bound(j));
                                        if((sg_l!=sg_r)||
                                           ((sg_l==CGAL::ZERO)&&
                                            (sg_r==CGAL::ZERO))){
                                                m[j]=i->second;
                                                --k;
                                        }
                                }
                        }
                }
                for(int l=0;l<isol.number_of_real_roots();++l)
                        *res++=std::make_pair(Algebraic(p,
                                                        isol.left_bound(l),
                                                        isol.right_bound(l)),
                                              m[l]);
                free(m);
                return res;
        }

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_1 &p,
                                  bool known_to_be_square_free,
                                  OutputIterator res)const{
                // TODO
                CGAL_error_msg("not implemented");
                return res;
        }

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_1 &p,
                                  const Bound &l,
                                  const Bound &u,
                                  OutputIterator res)const{
                // TODO 
                CGAL_error_msg("not implemented");
                return res;
        }

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_1 &p,
                                  bool known_to_be_square_free,
                                  const Bound &l,
                                  const Bound &u,
                                  OutputIterator res)const{
                // TODO 
                CGAL_error_msg("not implemented");
                return res;
        }

};  // Solve_1

} // namespace SimpleAK1
} // namespace CGAL

#endif // CGAL_RS_SIMPLE_FUNCTORS_1_H
