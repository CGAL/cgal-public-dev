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

#ifndef CGAL_RS_FUNCTORS_2_H
#define CGAL_RS_FUNCTORS_2_H

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/RS/algebraic_2.h>
#include <CGAL/Polynomial_type_generator.h>

namespace CGAL{
namespace RS3{

typedef CGAL::Polynomial_type_generator<Gmpz,2>::Type   RS_polynomial_2;
typedef int                                             Multiplicity;

template <class _Polynomial,class _Bound,class _Coefficient,class _AR1>
struct Construct_alg_2{
        typedef _Polynomial                             Polynomial;
        typedef _Bound                                  Bound;
        typedef _Coefficient                            Coefficient;
        typedef _AR1                                    Algebraic_real_1;
        typedef CGAL::RS3::Algebraic_2<Polynomial>      Algebraic_real_2;

        Algebraic_real_2 operator()(int x,int y)const{
                return Algebraic_real_2(x,y);
        };
        Algebraic_real_2 operator()(Bound x,Bound y)const{
                return Algebraic_real_2(x,y);
        };
        Algebraic_real_2 operator()(Coefficient x,Coefficient y)const{
                return Algebraic_real_2(x,y);
        };
        Algebraic_real_2 operator()(Algebraic_real_1 x,
                                    Algebraic_real_1 y)const{
                // TODO
                return Algebraic_real_2();
        };
        Algebraic_real_2 operator()(const Polynomial&,
                                    const Polynomial&,
                                    int)const{
                // TODO: do this in terms of Solve_2?
                return Algebraic_real_2();
        };
        Algebraic_real_2 operator()(const Polynomial&,
                                    const Polynomial&,
                                    const Bound&,const Bound&,
                                    const Bound&,const Bound&)const{
                // TODO: do this in terms of Solve_2?
                return Algebraic_real_2();
        };
}; // struct Construct_alg_2<_P,_B,_C,_AR1>

 template <class _AR2,class _Polynomial_1, class _Ptraits>
   struct Compute_polynomial_x_2{
     typedef _Ptraits                                Ptraits;
     typedef _AR2                                    Algebraic_real_2;
     typedef _Polynomial_1                           Polynomial_1;
     typedef typename Ptraits::Resultant Resultant;
     
     Polynomial_1 operator()( Algebraic_real_2 &a) const{
	  if (a.is_pol_x()){
	    return a.get_pol_x();
	  }
	  else{
	    Polynomial_1 res;
	    res = Resultant()(a.get_f(),a.get_g());
	    a.set_pol_x(res);
	    return res;
	  }
	  
        }
};

template <class _AR2,class _Polynomial_1, class _Ptraits>
struct Compute_polynomial_y_2{
  typedef _Ptraits                                Ptraits;
  typedef _AR2                                    Algebraic_real_2;
  typedef _Polynomial_1                           Polynomial_1;
  typedef typename Ptraits::Resultant Resultant;
  
  Polynomial_1 operator()(Algebraic_real_2 &a)const{
	  if (a.is_pol_y()){
	    return a.get_pol_y();
	  }
	  else{
	    Polynomial_1 res;
	    res = Resultant()(a.get_f(),a.get_g());
	    a.set_pol_y(res);
	    return res;
	    
	  }
        }
};

/* template <class _AR1,class _AR2,class _Polynomial_2,class _Bound> */
/* struct Isolate_2{ */
/*         typedef _AR1                                    Algebraic_real_1; */
/*         typedef _AR2                                    Algebraic_real_2; */
/*         typedef _Polynomial_2                           Polynomial_2; */
/*         typedef _Bound                                  Bound; */

/*         CGAL::cpp0x::array<AlgebraicKernel_d_1::Bound,4> */
/*         operator()(const Algebraic_real_2 &a,const Polynomial_2 &f)const{ */
/*                 // TODO */
/*         }; */

/*         CGAL::cpp0x::array<AlgebraicKernel_d_1::Bound,4> */
/*         operator()(const Algebraic_real_2 &a, */
/*                    const Polynomial_2 &f */
/*                    const Polynomial_2 &g)const{ */
/*                 // TODO */
/*         }; */
/* }; // struct Isolate_2 */

/* template <class _AR2,class _Polynomial_1,class _Bound> */
/* struct Isolate_x_2{ */
/*         typedef _AR2                                    Algebraic_real_2; */
/*         typedef _Polynomial_1                           Polynomial_1; */
/*         typedef _Bound                                  Bound; */

/*         std::pair<Bound,Bound> */
/*         operator()(const Algebraic_real_2 &a,const Polynomial_1 &f)const{ */
/*                 // TODO */
/*         }; */
/* }; // struct Isolate_x_2 */

/* template <class _AR2,class _Polynomial_1,class _Bound> */
/* struct Isolate_y_2{ */
/*         typedef _AR2                                    Algebraic_real_2; */
/*         typedef _Polynomial_1                           Polynomial_1; */
/*         typedef _Bound                                  Bound; */

/*         std::pair<Bound,Bound> */
/*         operator()(const Algebraic_real_2 &a,const Polynomial_1 &f)const{ */
/*                 // TODO */
/*         }; */
/* }; // struct Isolate_y_2 */

/* template <class _Ptraits> */
/* struct Is_square_free_2{ */
/*         typedef _Ptraits                                Ptraits; */
/*         typedef typename Ptraits::Polynomial_d          Polynomial; */

/*         bool operator()(const Polynomial &p)const{ */
/* 	  // TODO */
/*         }; */
/* }; // struct Is_square_free_2 */

template <class _Ptraits>
struct Is_coprime_2{
        typedef _Ptraits                                Ptraits;
        typedef typename Ptraits::Polynomial_d          Polynomial;
        typedef typename Ptraits::Gcd_up_to_constant_factor
                                                        Gcd;
        typedef typename Ptraits::Degree                Degree;

        bool operator()(const Polynomial &p1,const Polynomial &p2)const{
                Polynomial g=Gcd()(p1,p2);
                return (Degree()(g,0)==0&&Degree()(g,0)==0);
        };
}; // struct Is_coprime_2

/* TODO: for RS3, this functor may need to be specialized, in order to use */
/* fast functions for exact division */
template <class _Ptraits>
struct Make_coprime_2{
  typedef _Ptraits                                Ptraits;
  typedef typename Ptraits::Polynomial_d          Polynomial;
  typedef typename Ptraits::Integral_division_up_to_constant_factor
	IDiv;
  typedef typename Ptraits::Gcd_up_to_constant_factor          Gcd;
  typedef typename Ptraits::Degree          Degree;
  
        bool operator()(const Polynomial &p1,
                        const Polynomial &p2,
                        Polynomial &g,
                        Polynomial &q1,
                        Polynomial &q2)const{
	  g=Gcd()(p1,p2);
	  q1=IDiv()(p1,g);
	  q2=IDiv()(p2,g);
	  return (Degree()(g,0)==0&&Degree()(g,0)==0);
        };
}; // struct Make_coprime_2

template <class _Polynomial,class _Bound>
struct Solve_2{
        typedef _Polynomial                             Polynomial;
        typedef _Bound                                  Bound;

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial&,const Polynomial&,
                                  OutputIterator res)const{
                CGAL_error_msg("not implemented for this polynomial/bound");
                return res;
        }

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial&,const Polynomial&,
                                  const Bound&,const Bound&,
                                  const Bound&,const Bound&,
                                  OutputIterator res)const{
                CGAL_error_msg("not implemented for this polynomial/bound");
                return res;
        }
}; // struct Solve_2<_Polynomial,_Bound>

template <>
struct Solve_2<RS_polynomial_2,CGAL::Gmpfr>{
  typedef  RS_polynomial_2::NT             Polynomial;
        typedef CGAL::Gmpfr                             Bound;
        typedef CGAL::RS3::Algebraic_2<Polynomial>      Algebraic_real_2;
	typedef CGAL::RS3::Rur_2<Polynomial> rur_2;
	
        template <class OutputIterator>
	  OutputIterator operator()(const RS_polynomial_2 &f,const RS_polynomial_2 &g,
                                  OutputIterator res)const{
	  // TODO: solve the system {f=0,g=0}
	  // 1. call RS
	  // 2. store solutions in res
	  // (for the moment, we always return two hardcoded roots)
	  std::vector< rur_2 > rurs;
	  CGAL::RS3::decomposition_in_rurs_2(f,g, back_inserter(rurs));
	  int mult_sys = -1;
	  for (int i=0;i<rurs.size();i++)
	    {
	      std::vector< std::pair<CGAL::Gmpfi,CGAL::Gmpfi> > boxes;
	      CGAL::RS3::isolate_rurs_2(rurs[i], back_inserter(boxes));
	      for (int j=0;j<boxes.size();j++)
	      	{
		  Algebraic_real_2 alg = Algebraic_real_2(f, g, rurs[i],boxes[j].first,boxes[j].second);
	      	  *res++ = std::make_pair(alg,mult_sys);
	      	}
		
	    }
                return res;
        }

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial &f,const Polynomial &g,
                                  const Bound &xl,const Bound &xu,
                                  const Bound &yl,const Bound &yu,
                                  OutputIterator res)const{
                // TODO: solve inside the box [xl,xu]*[yl,yu]
                for(int i=6;i<8;++i)
                         *res++=std::make_pair(Algebraic_real_2(i,0),i+1);
                return res;
        }
	
 }; // struct Solve_2<RS_polynomial_2,CGAL::Gmpfr>
 
} // namespace RS3
} // namespace CGAL

#endif  // CGAL_RS_FUNCTORS_2_H
