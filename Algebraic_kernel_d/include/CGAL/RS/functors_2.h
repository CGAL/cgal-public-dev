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

#ifndef CGAL_RS_FUNCTORS_2_H
#define CGAL_RS_FUNCTORS_2_H

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/RS/algebraic_2.h>

namespace CGAL{
namespace RS3{

typedef CGAL::Gmpz                                      Innermost_coefficient;
typedef CGAL::Polynomial<Innermost_coefficient>         Polynomial_1;
typedef CGAL::Polynomial_type_generator<Innermost_coefficient,2>::Type
                                                        Polynomial_2;
typedef CGAL::RS3::Algebraic_2                          Algebraic_real_2;
typedef CGAL::Gmpfr                                     Bound;
typedef int                                             Multiplicity;

template <class Polynomial_>
struct Solve_2{
        typedef Polynomial_                                 Polynomial;

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial &f,const Polynomial &g,
                                  OutputIterator res)const{
                // TODO: solve the system {f=0,g=0}
                // 1. call RS
                // 2. store solutions in res
                // (for the moment, we always return two hardcoded roots)
                for(int i=0;i<2;++i)
                         *res++=std::make_pair(Algebraic_real_2(i,0),i+1);
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

}; // Solve_2

} // namespace RS3
} // namespace CGAL

#endif  // CGAL_RS_FUNCTORS_2_H
