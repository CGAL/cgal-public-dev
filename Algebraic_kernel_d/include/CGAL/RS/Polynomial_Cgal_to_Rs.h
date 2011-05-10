// Copyright (c) 2006-2009 Inria Lorraine (France). All rights reserved.
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
// $URL: svn+ssh://algerbya@scm.gforge.inria.fr/svn/cgal/branches/features/Algebraic_kernel_d-RS_bivariate-nancy/Algebraic_kernel_d/include/CGAL/RS/solve_1.h $
// $Id: solve_1.h 63060 2011-04-20 12:13:38Z penarand $
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>
// Author: Yacine Bouzidi <Bouzidi.yacine@gmail.com>

#ifndef CGAL_RS_POLY_FROM_H
#define CGAL_RS_POLY_FROM_H


#include <CGAL/basic.h>
#include <CGAL/RS/basic_1.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Polynomial/Polynomial_type.h>
#include <stdlib.h>
#include <gmp.h>

namespace CGAL {  
  
  
  
  typedef CGAL::Polynomial<CGAL::Gmpz> Polynomial_1;
  typedef CGAL::Polynomial_type_generator<CGAL::Gmpz,2>::Type Polynomial_2;
  
  // return the degree respectively in y and x of a bivariate polynomialx
  template<>
  std::pair<int,int> Degrees_Rs(const Polynomial_2& p) {
    
    int deg1 = p.degree();
    int deg2 = 0;
    for (int j=0;j<deg1;j++)
      {
	if (deg2 < p[j].degree())
	  deg2 = p[j].degree();
	
      }
    return std::make_pair(deg1,deg2);
  }

    // return the coefficients of a bivariate polynomial in the RS format 

  template<>
    MP_INT** Get_Coeff_Bi_For_Rs2(const Polynomial_2& p ) {
    std::pair<int,int> deg = Degrees_Rs(p);
    int Size_y = deg.first;
    int Size_x = deg.second;
    std::vector<Polynomial_1> v = p.coeffs();
      
      MP_INT** res = (MP_INT **)malloc((Size_y+1)*sizeof(MP_INT*));
      
      for(int i=0;i<=Size_y;i++)
	{
	  res[i] = (MP_INT *)malloc((Size_x+1)*sizeof(MP_INT));
	  for(int j=0;j<=Size_x;j++)
	    {
	      mpz_init(&res[i][j]);
	      if (j <= v[i].degree())
		mpz_set(&res[i][j],v[i][j].mpz());
	      else 
		mpz_set_ui(&res[i][j],0);
	    }
	}
      return res;
    }

  // return the coefficients of a univariate polynomial in the RS format (ie MP_INT)
  
  template<>
    MP_INT* Get_Coeff_For_Rs1(const Polynomial_1& p) {
    int Size = (p.coeffs()).size();
    std::vector<CGAL::Gmpz> v = p.coeffs();
    
    MP_INT *res = (MP_INT *)malloc(Size*sizeof(MP_INT));
    
    for(int j=0;j<Size;j++)
      {
	
	mpz_init(&res[j]);
	mpz_set(&res[j],(v[j].mpz()));
	
      }
    
    return res;
  }
    
  
  
} // namespace CGAL


#endif  //CGAL_RS_POLY_FROM_H
