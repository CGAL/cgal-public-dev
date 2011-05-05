// Copyright (c) 2006-2008 Inria Lorraine (France). All rights reserved.
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
// $URL: svn+ssh://algerbya@scm.gforge.inria.fr/svn/cgal/branches/features/Algebraic_kernel_d-RS_bivariate-nancy/Algebraic_kernel_d/include/CGAL/RS/compare_1.h $
// $Id: compare_1.h 63060 2011-04-20 12:13:38Z penarand $
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>
// Author: Yacine Bouzidi <Bouzidi.yacine@gmail.com>

#ifndef CGAL_RS_COMPARE_1_H
#define CGAL_RS_COMPARE_1_H

#include <mpfr.h>
#include <CGAL/RS/New_refine_rs_1.h>


// default refinement and sign functions
//#define CGALRS_REFINE_N(A,N)        CGAL::NewAlg::refine_1(A,N)
#define CGALRS_REFSTEPS             4
//#define CGALRS_SIGNAT(P,M)          RSSign::signat(P,M)

namespace CGAL{
  
  namespace NewAlg{
    
  namespace RS_COMPARE{
    
    
    //typedef CGAL::NewAlg::New_Algebraic_1 New_Algebraic_1;
    // compare two algebraic numbers, knowing they are not equal
    template < class Coefficient_ >
      Comparison_result 
      compare_1_unequal(const CGAL::NewAlg::New_Algebraic_1<Coefficient_> &r1,const CGAL::NewAlg::New_Algebraic_1<Coefficient_> &r2){
      /*typedef _Gcd_policy     Gcd;
        mp_prec_t prec1=r1.get_prec();
        mp_prec_t prec2=r2.get_prec();*/
      CGAL::NewAlg::refine_1<Coefficient_>(r1,CGALRS_REFSTEPS);
      CGAL::NewAlg::refine_1<Coefficient_>(r2,CGALRS_REFSTEPS);
        while(r1.overlaps(r2)){
                /*if(prec1<prec2 || r2.lefteval()==ZERO){
                        CGALRS_REFINE_N(r1,CGALRS_REFSTEPS);
                        prec1=r1.get_prec();
                }else{
                        CGALRS_REFINE_N(r2,CGALRS_REFSTEPS);
                        prec2=r2.get_prec();
                }
                */
                CGAL::NewAlg::refine_1<Coefficient_>(r1,CGALRS_REFSTEPS);
                CGAL::NewAlg::refine_1<Coefficient_>(r2,CGALRS_REFSTEPS);
        }
        return(mpfr_less_p(r1.right(),r2.left())?SMALLER:LARGER);
    }
 
    template <class Coefficient_>
      Comparison_result
      compare_1(const CGAL::NewAlg::New_Algebraic_1<Coefficient_> &r1,const CGAL::NewAlg::New_Algebraic_1<Coefficient_> &r2){
      // typedef _Gcd_policy     Gcd;
      /*if(r1.pol()==r2.pol())
	return(r1.nr()!=r2.nr()?(r1.nr()<r2.nr()?SMALLER:LARGER):EQUAL);
      */
      
      typedef typename CGAL::Polynomial_type_generator<Coefficient_,1>::Type Polynomial_1;
      typedef CGAL::Polynomial_traits_d< Polynomial_1 >
	Polynomial_traits_1;
      
      if(mpfr_lessequal_p(r1.left(),r2.left())){
	if(mpfr_less_p(r1.right(),r2.left()))
	  return SMALLER;
      }else{
	if(mpfr_less_p(r2.right(),r1.left()))
		  return LARGER;
        }
      typename Polynomial_traits_1::Gcd_up_to_constant_factor gcd_utcf;
      Polynomial_1 gcd=gcd_utcf(r1.pol(),r2.pol());
      if(!gcd.degree())
	return compare_1_unequal<Coefficient_>(r1,r2);

        Sign sleft,sright;
        if(mpfr_greater_p(r1.left(),r2.left()))
	  sleft=gcd.sign_at(Gmpfr(r1.left()));
	
        else
	  sleft=gcd.sign_at(Gmpfr(r2.left()));
	
        if(sleft==ZERO)
                return EQUAL;
        if(mpfr_less_p(r1.right(),r2.right()))
	  sright=gcd.sign_at(Gmpfr(r1.right()));
                
        else
	  sright=gcd.sign_at(Gmpfr(r2.right()));
                
        if(sleft!=sright)
                return EQUAL;
        else
	  return RS_COMPARE::compare_1_unequal<Coefficient_>(r1,r2);
}


    
  } // namespace RS_COMPARE
  } // namespace NewAlg
} // namespace CGAL

//#undef CGALRS_REFINE_N
#undef CGALRS_REFSTEPS


#endif  // CGAL_RS_COMPARE_1_H
