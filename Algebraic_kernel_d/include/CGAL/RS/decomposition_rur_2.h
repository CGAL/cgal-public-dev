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

//#ifndef CGAL_RS_SOLVE_1_H
//#define CGAL_RS_SOLVE_1_H
//#include <stdlib.h>
#include <CGAL/basic.h>
#include <CGAL/RS/basic_1.h>
//#include <CGAL/RS/dyadic.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Polynomial.h>
#include <CGAL/RS/rs_calls_1.h>
#include <CGAL/Gmpfi.h>
#include <stdlib.h>

namespace CGAL {  
  
  
 
  
  typedef  CGAL::Polynomial_type_generator<CGAL::Gmpz,2>::Type Polynomial_2;
  
    std::vector< std::vector< std::vector<MP_INT> > > decomposition_in_rurs_2( const Polynomial_2 &p1, const Polynomial_2 &p2,unsigned int prec=CGAL_RS_DEF_PREC){
    
    rs_init_rs();
    rs_reset_all();
    create_rs_bisys(p1.Get_Coeff_Bi_For_Rs2(),p1.Degrees_Rs().first,p1.Degrees_Rs().second, p2.Get_Coeff_Bi_For_Rs2(),p2.Degrees_Rs().first, p2.Degrees_Rs().second);
    set_rs_precisol(prec);
    set_rs_verbose(0);
    rs_run_algo(CGALRS_CSTR("RURBIV"));
    return Rurs_sys_list();
    
    
  }
  
  
  
} // namespace CGAL

//#endif  // CGAL_RS_SOLVE_1_H
