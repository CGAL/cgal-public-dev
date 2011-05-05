// Copyright (c) 2009-2010 Inria Lorraine (France). All rights reserved.
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
// $URL: svn+ssh://algerbya@scm.gforge.inria.fr/svn/cgal/branches/features/Algebraic_kernel_d-RS_bivariate-nancy/Algebraic_kernel_d/include/CGAL/RS/sign_rs_1.h $
// $Id: sign_rs_1.h 63060 2011-04-20 12:13:38Z penarand $
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>
// Author: Yacine Bouzidi <Bouzidi.yacine@gmail.com>

#ifndef CGAL_RS_SIGN_1_RS_H
#define CGAL_RS_SIGN_1_RS_H
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>

#include <CGAL/assertions.h>

namespace RS3{
  
  CGAL::Sign sign_1(const CGAL::NewAlg::New_Algebraic_1<CGAL::Gmpz>::Polynomial_1 &p, const CGAL::NewAlg::New_Algebraic_1<CGAL::Gmpz> &a){
    
    mpz_t* list_constr[]={p.Get_Coeff_For_Rs()};
    int list_degs[]={p.degree()};
    mpfi_t* sols_u=(mpfi_t*)malloc(sizeof(mpfi_t));
    *(sols_u[0])=*(a.mpfi());
    mpfi_t **vals_constr=(mpfi_t**)malloc(sizeof(mpfi_t*));
    vals_constr[0]=(mpfi_t*)malloc(sizeof(mpfi_t));
    mpfi_init(vals_constr[0][0]);
    rs3_refine_eval_u(a.pol().Get_Coeff_For_Rs(),a.pol().degree(),NULL,0,
		      (const mpz_t **)list_constr,list_degs,
		      1,sols_u,1,
		      vals_constr,NULL,MPFR_PREC_MIN,1,1,1);
    if(mpfr_zero_p(&vals_constr[0][0]->left)!=0 &&
       mpfr_zero_p(&vals_constr[0][0]->right)!=0){
      return CGAL::ZERO;
    }
    CGAL_assertion(mpfi_has_zero(vals_constr[0][0])<=0);
    if(mpfi_is_pos(vals_constr[0][0])){
      return CGAL::POSITIVE;
    }
        return CGAL::NEGATIVE;
  }
  
} // namespace RS3

#endif  // CGAL_RS_SIGN_1_RS_H
