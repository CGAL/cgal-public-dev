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
// $URL: svn+ssh://algerbya@scm.gforge.inria.fr/svn/cgal/branches/features/Algebraic_kernel_d-RS_bivariate-nancy/Algebraic_kernel_d/include/CGAL/RS/algebraic_1_operators.h $
// $Id: algebraic_1_operators.h 61907 2011-03-22 10:11:01Z penarand $
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>
// Author: Yacine Bouzidi <Bouzidi.yacine@gmail.com>

#ifndef CGAL_RS_ALGEBRAIC_1_OPERATORS_H
#define CGAL_RS_ALGEBRAIC_1_OPERATORS_H

namespace CGAL{
  
  namespace NewAlg{
    
    template<class Coefficient_> inline
      New_Algebraic_1<Coefficient_> New_Algebraic_1<Coefficient_>::operator+()const{
      return *this;
    }

    template<class Coefficient_> inline
      New_Algebraic_1<Coefficient_> New_Algebraic_1<Coefficient_>::operator-()const{
      mpfi_t inv;
      mpfi_init2(inv,mpfi_get_prec(mpfi()));
      mpfi_neg(inv,mpfi());
      New_Algebraic_1 *inverse=new New_Algebraic_1(inv,
						   pol(),// dont know if an inverse polynomial function existe for Polynomial_1
						   nr(),
						   mult(),
						   NULL,
						   NULL,
						   -lefteval());
      return *inverse;
    }
 
  } // namespace NewAlg
} // namespace CGAL

#endif  // CGAL_RS_ALGEBRAIC_1_OPERATORS_H
