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

#include <stdlib.h>


#ifndef CGAL_RS_RUR_2_H
#define CGAL_RS_RUR_2_H

namespace CGAL{

namespace RS3{

#include <CGAL/Handle_for.h>

template <class Polynomial_>
class Rur_2_rep{
        public:
  typedef Polynomial_                             Polynomial_1;
  // the four polynomials that form the RUR of a bivariate system
  Polynomial_1 _f,_g,_g_1,_g_2;
  // the separating element (a,b) representing the polynomial a*(first variable)+b(second variable)
  std::pair<int,int> _sep_elm;
  // the multiplicity of the rur which is the multiplicity of all the solutions since the decomposition is performed with respect to the multiplcities 
  int _multiplicity;
  
  Rur_2_rep(){};
  ~Rur_2_rep(){};
 private:
  Rur_2_rep(const Rur_2_rep&);
  Rur_2_rep& operator=(const Rur_2_rep);
 };
 
template <class Polynomial_>
  class Rur_2:Handle_for<Rur_2_rep<Polynomial_> >{
        private:
        typedef Polynomial_                             Polynomial_1;
        typedef Rur_2_rep<Polynomial_1>                 Base;

        public:

        // default constructor (gives an inconsistent RUR)
        // TODO: remove this constructor?
        Rur_2(){};

        // constructor from the four polynomials the separating element and the multiplicity
        Rur_2(const Polynomial_1 &f,const Polynomial_1 &g,
              const Polynomial_1 &g_1,const Polynomial_1 &g_2, std::pair<int,int> &sep_elm, int multiplicity){
	 
	  this->ptr()->_f = f;
	 
	  this->ptr()->_g=g;
	 
	  this->ptr()->_g_1=g_1;
	 
	  this->ptr()->_g_2=g_2;
	 
	  this->ptr()->_sep_elm=sep_elm;
	 
	  this->ptr()->_multiplicity=multiplicity;
	 
        };

        // get members as const
        const Polynomial_1& get_f()const{return (this->ptr())->_f;};
        const Polynomial_1& get_g()const{return (this->ptr())->_g;};
        const Polynomial_1& get_g1()const{return (this->ptr())->_g_1;};
        const Polynomial_1& get_g2()const{return (this->ptr())->_g_2;};
	const std::pair<int,int> get_sep_elm()const{return (this->ptr())->_sep_elm;};
	const Polynomial_1& get_multiplicity()const{return (this->ptr())->_multiplicity;};
	
	
 }; // class Rur_2
 
// write a RUR to a stream
// TODO: binary mode
/* template <class Polynomial_> */
/* inline std::ostream& operator<<(std::ostream &o,const Rur_2<Polynomial_> &r){ */
/*         if(is_pretty(o)){ */
/*                 o<<"Rur_2("<<r.get_f()<<','<<r.get_g()<<','<< */
/*                         r.get_g1()<<','<<r.get_g2()<<')'; */
/*         }else{ */
/*                 o<<'['<<r.get_f()<<','<<r.get_g()<<','<< */
/*                         r.get_g1()<<','<<r.get_g2()<<']'; */
/*         } */
/*         return o; */
/* } */

} // namespace RS3
} // namespace CGAL

#endif // CGAL_RS_RUR_2_H
