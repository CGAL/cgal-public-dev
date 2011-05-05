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
// $URL: svn+ssh://algerbya@scm.gforge.inria.fr/svn/cgal/branches/features/Algebraic_kernel_d-RS_bivariate-nancy/Algebraic_kernel_d/include/CGAL/RS/algebraic_1_comparisons.h $
// $Id: algebraic_1_comparisons.h 61907 2011-03-22 10:11:01Z penarand $
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>
// Author: Yacine Bouzidi <Bouzidi.yacine@gmail.com>

#ifndef CGAL_RS_ALGEBRAIC_1_COMPARISONS_H
#define CGAL_RS_ALGEBRAIC_1_COMPARISONS_H

#include <boost/config.hpp>

#include <CGAL/RS/New_algebraic_1.h>
#include <CGAL/RS/New_compare_1.h>

namespace CGAL{

  namespace NewAlg{

    
    template< class Coefficient_ >inline
      bool operator<(const New_Algebraic_1<Coefficient_> &n1,const New_Algebraic_1<Coefficient_> &n2){
      
      return(CGAL::NewAlg::RS_COMPARE::compare_1<Coefficient_>(n1,n2)==SMALLER);
    }
    
    
    template< class Coefficient_ >inline
      bool operator==(const New_Algebraic_1<Coefficient_> &n1,const New_Algebraic_1<Coefficient_> &n2){
      
      return((CGAL::NewAlg::RS_COMPARE::compare_1<Coefficient_>(n1,n2))==EQUAL);
    }
    
    
    template< class Coefficient_ >inline
      New_Algebraic_1<Coefficient_> min BOOST_PREVENT_MACRO_SUBSTITUTION (const New_Algebraic_1<Coefficient_> &a,const New_Algebraic_1<Coefficient_> &b){
      return (a<b?a:b);
    }
    
    
    template< class Coefficient_ >inline
      New_Algebraic_1<Coefficient_> max BOOST_PREVENT_MACRO_SUBSTITUTION (const New_Algebraic_1<Coefficient_> &a,const New_Algebraic_1<Coefficient_> &b){
      return (a>b?a:b);
    }
    
  }// namespace NewAlg
} // namespace CGAL

#endif  // CGAL_RS_ALGEBRAIC_1_COMPARISONS_H
