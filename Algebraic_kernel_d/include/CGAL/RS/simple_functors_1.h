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
                Isolator isol(p);
                for(int l=0;l<isol.number_of_real_roots();++l)
                        *res++=Algebraic(p,
                                         isol.left_bound(l),
                                         isol.right_bound(l));
                return res;
        }

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_1 &p,
                                  const Bound &l,
                                  const Bound &u,
                                  OutputIterator res)const{
                typedef std::vector<Algebraic>                  RV;
                typedef typename RV::iterator                   RVI;
                typedef std::pair<Polynomial_1,int>             PM;
                typedef std::vector<PM>                         PMV;
                typedef typename PMV::iterator                  PMVI;
                CGAL_precondition_msg(l<=u,
                                      "left bound must be <= right bound");
                RV roots; // all roots of the polynomial
                this->operator()(p,false,std::back_inserter(roots));
                size_t nb_roots=roots.size();
                // indices of the first and last roots to be reported:
                int index_l=0,index_u;
                while(index_l<nb_roots&&roots[index_l]<l)
                        ++index_l;
                CGAL_assertion(index_l<=nb_roots);
                if(index_l==nb_roots)
                        return res;
                index_u=index_l;
                while(index_u<nb_roots&&roots[index_u]<u)
                        ++index_u;
                CGAL_assertion(index_u<=nb_roots);
                if(index_u==index_l)
                        return res;
                // now, we have to return roots in [index_l,index_u)
                PMV sfv;
                Sqfr()(p,std::back_inserter(sfv)); // square-free fact. of p
                // array to store the multiplicities
                int *m=(int*)calloc(nb_roots,sizeof(int));
                // we iterate over all the pairs <root,mult> and match the
                // roots in the interval [index_l,index_u)
                for(PMVI i=sfv.begin();i!=sfv.end();++i){
                        int k=Degree()(i->first);
                        Signat signof(i->first);
                        for(int j=index_l;k&&j<index_u;++j){
                                if(!m[j]){
                                        CGAL::Sign sg_l=
                                                signof(roots[j].get_left());
                                        CGAL::Sign sg_r=
                                                signof(roots[j].get_right());
                                        if((sg_l!=sg_r)||
                                           ((sg_l==CGAL::ZERO)&&
                                            (sg_r==CGAL::ZERO))){
                                                m[j]=i->second;
                                                --k;
                                        }
                                }
                        }
                }
                for(int l=index_l;l<index_u;++l)
                        *res++=std::make_pair(roots[l],m[l]);
                free(m);
                return res;
        }

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_1 &p,
                                  bool known_to_be_square_free,
                                  const Bound &l,
                                  const Bound &u,
                                  OutputIterator res)const{
                typedef std::vector<Algebraic>                  RV;
                typedef typename RV::iterator                   RVI;
                CGAL_precondition_msg(l<=u,
                                      "left bound must be <= right bound");
                RV roots;
                this->operator()(p,
                                 known_to_be_square_free,
                                 std::back_inserter(roots));
                for(RVI it=roots.begin();it!=roots.end();it++)
                        if(*it>=l&&*it<=u)
                                *res++=*it;
                return res;
        }

};  // Solve_1

} // namespace SimpleAK1
} // namespace CGAL

#endif // CGAL_RS_SIMPLE_FUNCTORS_1_H
