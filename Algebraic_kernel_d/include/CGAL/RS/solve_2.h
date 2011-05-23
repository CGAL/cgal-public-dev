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
// $URL: svn+ssh://algerbya@scm.gforge.inria.fr/svn/cgal/branches/features/Algebraic_kernel_d-RS_bivariate-nancy/Algebraic_kernel_d/include/CGAL/RS/decomposition_rur_2.h $
// $Id: solve_1.h 63060 2011-04-20 12:13:38Z penarand $
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>
// Author: Yacine Bouzidi <Bouzidi.yacine@gmail.com>

#ifndef CGAL_RS_DECOMPOSE_2_H
#define CGAL_RS_DECOMPOSE_2_H

#include <CGAL/basic.h>
#include <CGAL/RS/basic_1.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/RS/rs3_calls.h>
#include <CGAL/RS/rur_2.h>




namespace CGAL {  
  
  namespace RS3 {
    
    
    template<class OutputIterator, class Polynomial_>
      OutputIterator decomposition_in_rurs_2( const CGAL::Polynomial<Polynomial_> &p1, const CGAL::Polynomial<Polynomial_> &p2, OutputIterator res ,unsigned int prec=CGAL_RS_DEF_PREC){
      
      typedef  Polynomial_ Polynomial_1;
      typedef  CGAL::Polynomial<Polynomial_1> Polynomial_2;
      typedef  CGAL::Polynomial_traits_d<Polynomial_1> Polynomial_traits_1;
      typedef  Rur_2<Polynomial_1> rur_2;
      typename  Polynomial_traits_1::Construct_polynomial Construct_polynomial;
      typename  Polynomial_traits_1::Differentiate Differentiate; 
      typename  Polynomial_traits_1::Shift Shift; 
      
      std::vector< std::vector< std::vector<CGAL::Gmpz> > > vect;
      RS3::init_solver();
      RS3::create_rs_bisys(p1,p2);
      set_rs_precisol(prec);
      set_rs_verbose(3);
      rs_run_algo(CGALRS_CSTR("RURBIV"));
      vect = RS3::Rurs_sys_list();
      
    
    // construct the set of the Rurs by extracting the corresponding polynomials 
    for (int k=0; k<vect.size();k=k+2)
      {
	
	
    	// construct the extension polynomial of the RUR
    	Polynomial_1 _f = Construct_polynomial(vect[k][0].begin(),vect[k][0].end());
	
	// construct the denominator of the X,Y coordinates (square free of _f)
	Polynomial_1 _g = Differentiate(_f);
	
	// construct the numerator of the X coordinate
	Polynomial_1 _gy;
    	if (vect[k+1][0].size() != 0)
    	  _gy = Construct_polynomial(vect[k+1][0].begin(),vect[k+1][0].end());
    	else
    	  _gy = Shift(_g, 1);
    	
	// construct the numerator of the y coordinate
	Polynomial_1 _gx;
    	if (vect[k+1][1].size() != 0)
    	  _gx = Construct_polynomial(vect[k+1][1].begin(),vect[k+1][1].end());
    	else
    	  _gx = Shift(_g, 1);
    	
	// construct the separating element
	std::pair<int,int> sep_elm;
	if (vect[k+1][2].size() == 1)
	  sep_elm = std::make_pair(0,mpz_get_si(vect[k+1][2][0].mpz()));
	else
    	  sep_elm = std::make_pair(mpz_get_si(vect[k+1][2][1].mpz()), mpz_get_si(vect[k+1][2][0].mpz()));
	
	
    	// the multiplicity of the current RUR
	int multiplicity = mpz_get_ui(vect[k+1][3][0].mpz());
    	
	// construct the RUR and add it to the list.
	*res++ = rur_2(_f,_g,_gx,_gy,sep_elm,multiplicity);
	
      }
    
    
    return res;
    
  }



  //isolate_2 : return disjoint boxes namely pair<mpfi,mpfi> representing the solutions of the rurs
  // TODO : construct for each pair of mpfi the corresponding algebraic_2 
    template< class OutputIterator, class Polynomial_ >
      OutputIterator isolate_rurs_2(CGAL::RS3::Rur_2<Polynomial_>& rur,OutputIterator res, unsigned int prec=CGAL_RS_DEF_PREC){
    
      typedef CGAL::RS3::Rur_2<Polynomial_> rur_2;
    // the solutions are expressed as a vector of pair<Gmpfi, Gmpfi>
    //std::vector< std::pair<CGAL::Gmpfi, CGAL::Gmpfi> > solutions;
    
    // loop on the set of the rurs
    // for(InputIterator it = begin; it != end; it++){
      rs_reset_all();
      set_rs_verbose(3);
      set_rs_precisol(prec);
      // create the Rs object which represente the CGAL rur
      create_rs_rur(rur);
      // isolate the current rur stored inside Rs
      rs_run_algo(CGALRS_CSTR("RISOLE"));
      // add the solutions of the current rur to the vector of the solutions
      affiche_sols_eqs(res);
      
      
  
    }
  
  } // namespace RS3  

} // namespace CGAL

#endif //CGAL_RS_DECOMPOSE_2_H
