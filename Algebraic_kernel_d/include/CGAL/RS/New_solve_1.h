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

#ifndef CGAL_RS_SOLVE_1_H
#define CGAL_RS_SOLVE_1_H

#include <CGAL/basic.h>
#include <CGAL/RS/basic_1.h>
#include <CGAL/RS/dyadic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/RS/New_algebraic_1.h>
#include <CGAL/RS/rs_calls_1.h>
#include <CGAL/Gmpfi.h>
#include <vector>

namespace CGAL{

  //class RS_polynomial_1;
  //class Algebraic_1;

  // solve given the precision, returns the number of roots
  typedef CGAL::NewAlg::New_Algebraic_1<Gmpz> New_Algebraic_1;
  typedef New_Algebraic_1::Polynomial_1 Polynomial_1;


    inline int solve_1(mpfi_ptr *x,
		       const Polynomial_1 &p1,
		       unsigned int prec=CGAL_RS_DEF_PREC){
      
      rs_init_rs();
      rs_reset_all();
      create_rs_upoly2(p1.Get_Coeff_For_Rs1(),p1.degree(),rs_get_default_up());
      set_rs_precisol(prec);
      set_rs_verbose(3);
      rs_run_algo(CGALRS_CSTR("UISOLE"));
      return affiche_sols_eqs(x);
    
    }

// calculate the sign of a polynomial evaluated at the root of another
    inline Sign sign_rs_1(const Polynomial_1 &p1,
			  const New_Algebraic_1 &a,
			  unsigned int prec=CGAL_RS_MIN_PREC){
      mpz_t **constr;
      int *degs;
      CGAL_assertion(a.is_consistent());
      rs_reset_all ();
        // tell RS to find the roots of this polynomial
      create_rs_upoly (a.pol().Get_Coeff_For_Rs(), a.pol().degree(),
                        rs_get_default_up());
        // the constraint vector will have just one element
        constr = (mpz_t**)malloc(sizeof(mpz_t*));
        *constr = p1.Get_Coeff_For_Rs();
        degs = (int*)malloc(sizeof(int));
        *degs = p1.degree ();
        create_rs_uconstr (constr, degs, rs_get_default_ineqs_u ());
        set_rs_precisol (prec);
        set_rs_verbose (CGAL_RS_VERB);
        rs_run_algo(CGALRS_CSTR("UISOLES"));
        return affiche_signs_constr(a.nr());
}

} // namespace CGAL

#endif  // CGAL_RS_SOLVE_1_H
