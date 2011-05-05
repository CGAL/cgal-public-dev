// Copyright (c) 2009 Inria Lorraine (France). All rights reserved.
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
// $URL: svn+ssh://algerbya@scm.gforge.inria.fr/svn/cgal/branches/features/Algebraic_kernel_d-RS_bivariate-nancy/Algebraic_kernel_d/include/CGAL/RS/refine_rs_1.h $
// $Id: refine_rs_1.h 63060 2011-04-20 12:13:38Z penarand $

// Author: Luis Peñaranda <luis.penaranda@gmx.com>
// Author: Yacine Bouzidi <Yacine.Bouzidi@gmail.com>

#ifndef CGAL_RS_REFINE_1_RS_H
#define CGAL_RS_REFINE_1_RS_H

#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <CGAL/Polynomial/Polynomial_type.h>
#include <CGAL/assertions.h>
#include <CGAL/RS/New_algebraic_1.h>
#include <rs3_fncts.h>

namespace CGAL{

namespace NewAlg{
  
  
  
  template< class Coefficient_> inline
    void refine_1(const New_Algebraic_1<Coefficient_>& a,unsigned int s=10000){
    const mpz_t* coef = (a.pol()).Get_Coeff_For_Rs(); 
    CGAL_precondition(a.inf()<=a.sup());
    rs3_refine_u_root((mpfi_ptr)a.mpfi(),
		      coef,
		      a.pol().degree(),
		      mpfi_get_prec(a.mpfi())+s,
		      0,
		      0);
    CGAL_assertion(a.inf()<=a.sup());
  }
  
  
} // namespace NewAlg

} // namespace CGAL
#endif  // CGAL_RS_REFINE_1_RS_H
