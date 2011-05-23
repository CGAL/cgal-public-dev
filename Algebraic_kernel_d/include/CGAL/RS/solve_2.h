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
#include <CGAL/Gmpfi.h>



namespace CGAL {  
  
  namespace RS3 {
    
  typedef  CGAL::Polynomial<CGAL::Gmpz> Polynomial_1;
  typedef  CGAL::Polynomial<Polynomial_1> Polynomial_2;
  typedef  CGAL::Polynomial_traits_d<Polynomial_1> Polynomial_traits_1;
  typedef  CGAL::RS3::Rur_2<Polynomial_1> rur_2;

  template<class OutputIterator>
    OutputIterator decomposition_in_rurs_2( const Polynomial_2 &p1, const Polynomial_2 &p2, OutputIterator res ,unsigned int prec=CGAL_RS_DEF_PREC){
    
    
    std::vector< std::vector< std::vector<CGAL::Gmpz> > > vect;
    RS3::init_solver();
    RS3::create_rs_bisys<Polynomial_2>(p1,p2);
    set_rs_precisol(prec);
    set_rs_verbose(4);
    rs_run_algo(CGALRS_CSTR("RURBIV"));
    vect = RS3::Rurs_sys_list();
    affiche_bivariate_decomp_raw();
    
    // construct the set of the Rurs by extracting the corresponding polynomials 
    for (int k=0; k<vect.size();k=k+2)
      {
	
	
    	// construct the extension polynomial of the RUR
    	Polynomial_1 _f = Polynomial_traits_1::Construct_polynomial()(vect[k][0].begin(),vect[k][0].end());
	fprintf(stderr,"extension \n");
	for (int i=0; i<= _f.degree();i++){
	  
	  mpz_out_str(stderr,10,_f[i].mpz());
	  fprintf(stderr," ");
	}
	// construct the denominator of the X,Y coordinates (square free of _f)
    	printf("\n");
	Polynomial_1 _g = Polynomial_traits_1::Differentiate()(_f);
	fprintf(stderr,"denominateur \n");
    	for (int i=0; i<= _g.degree();i++){
	  
	  mpz_out_str(stderr,10,_g[i].mpz());
	   fprintf(stderr," ");
	}
	printf("fin \n");
	// construct the numerator of the X coordinate
	Polynomial_1 _gy;
    	if (vect[k+1][0].size() != 0)
    	  _gy = Polynomial_traits_1::Construct_polynomial()(vect[k+1][0].begin(),vect[k+1][0].end());
    	else
    	  _gy = Polynomial_traits_1::Shift()(_g, 1);
    	
	// construct the numerator of the y coordinate
	Polynomial_1 _gx;
    	if (vect[k+1][1].size() != 0)
    	  _gx = Polynomial_traits_1::Construct_polynomial()(vect[k+1][1].begin(),vect[k+1][1].end());
    	else
    	  _gx = Polynomial_traits_1::Shift()(_g, 1);
    	
	// construct the separating element
	std::pair<int,int> sep_elm;
	if (vect[k+1][2].size() == 1)
	  sep_elm = std::make_pair(0,mpz_get_si(vect[k+1][2][0].mpz()));
	else
    	  sep_elm = std::make_pair(mpz_get_si(vect[k+1][2][1].mpz()), mpz_get_si(vect[k+1][2][0].mpz()));
	
	
    	// the multiplicity of the current RUR
	int multiplicity = mpz_get_ui(vect[k+1][3][0].mpz());
    	
	// construct the RUR and add it to the list.
	*res++ = CGAL::RS3::rur_2(_f,_g,_gx,_gy,sep_elm,multiplicity);
	
      }
    
    
    return res;
    
  }



  //isolate_2 : return disjoint boxes namely pair<mpfi,mpfi> representing the solutions of the rurs
  // TODO : construct for each pair of mpfi the corresponding algebraic_2 
  template< class InputIterator >
    int isolate_rurs_2(InputIterator begin, InputIterator end, unsigned int prec=CGAL_RS_DEF_PREC){
    
    // the solutions are expressed as a vector of pair<Gmpfi, Gmpfi>
    std::vector< std::pair<CGAL::Gmpfi, CGAL::Gmpfi> > solutions;
    
    // loop on the set of the rurs
    for(InputIterator it = begin; it != end; it++){
      rs_reset_all();
      set_rs_precisol(prec);
      // create the Rs object which represente the CGAL rur
      create_rs_rur<rur_2>(*it);
      // isolate the current rur stored inside Rs
      rs_run_algo(CGALRS_CSTR("RISOLE"));
      // add the solutions of the current rur to the vector of the solutions
      affiche_sols_eqs(solutions);
      
    }
    
  }
  
  } // namespace RS3  
  
} // namespace CGAL

#endif //CGAL_RS_DECOMPOSE_2_H
