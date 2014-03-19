// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Salahuddin Pasha <s9mdpash@stud.uni-saarland.de>
//
// ============================================================================
#ifndef CGAL_POLYNOMIAL_SAGE_D_H
#define CGAL_POLYNOMIAL_SAGE_D_H

#include <CGAL/basic.h>
#include <functional>
#include <list>
#include <vector>
#include <utility>

#include <CGAL/Polynomial/Sage/Polynomial_Sage.h>

#define CGAL_POLYNOMIAL_SAGE_D_BASE_TYPEDEFS

namespace CGAL {

namespace internal {

} // namespace internal


//template< class Polynomial >
class Polynomial_Sage_d {

 private:
  std::string sagePolynomial;
  
 public:
  Polynomial_Sage_d(std::string inputPolynomial) {
    sagePolynomial = inputPolynomial;
  }

  const std::string str() {
    return sagePolynomial;
  }
  
};

} //namespace CGAL

std::ostream &operator<<(std::ostream &os, CGAL::Polynomial_Sage_d &m) { 
  return os << m.str();
}



std::string degree(CGAL::Polynomial_Sage_d &polynomialSage) 
{
  std::ostringstream oStringForSage;

  oStringForSage << "R.<y> = PolynomialRing(ZZ)\nR.<x> = PolynomialRing(ZZ)\np=(" << polynomialSage << ").degree()";

  std::string stringForSage = oStringForSage.str();
  std::string dataFromSage = getDataFromSage( stringForSage );
  return dataFromSage;
}

#endif // CGAL_POLYNOMIAL_SAGE_D_H
