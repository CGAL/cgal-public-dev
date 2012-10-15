// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id:
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//
//

#ifndef CGAL_VDL3_BISECTOR_3_H
#define CGAL_VDL3_BISECTOR_3_H

#include<CGAL/compute_hash_value.h>

namespace CGAL {

namespace VDOL_3 {

//! Definition of the bisector.
/*! 
  We represent a bisector of to lines in 3D as two polynomials with 3 
  variables. The first is the representation in R^3 (quadric) and the second
  is the representation in the parameter space with respect to the base line.
  It is a surface of degree 4 (This is still of degree 2 in each variable).
  
  Basically, bisectors are the surfaces in our case. The problem is that these
  surfaces are not regular surfaces they can approach infinity at single line
  (which the envelope can not handle). That is the reason we add another 
  symbolic surface at infinity such an intersection with a bisector can be
  reported as line at infinity. Another advantage is that any feature 
  containing the symbolic surface will be labeled as at infinity.
  \todo Move this class to another file.
*/

template <typename SVCET_3> 
class Bisector
  : boost::totally_ordered< Bisector<SVCET_3> > 
{
  typedef typename SVCET_3::Poly_int_2 Poly_int_2;
  typedef typename SVCET_3::Poly_int_3 Poly_int_3;
  typedef typename SVCET_3::Poly_rat_3 Poly_rat_3;
  typedef typename SVCET_3::Poly_lazy_rat_3 Poly_lazy_rat_3;

  typedef typename SVCET_3::Line_3 Line_3;
  typedef typename SVCET_3::X_monotone_curve_2 X_monotone_curve_2;

  typedef Bisector<SVCET_3> Self;
public:
  typedef std::list< std::pair< X_monotone_curve_2, Oriented_side > >
  Boundary_container;
  
protected:
  Line_3     m_line;
  Poly_int_3 m_bisector;
  Poly_int_3 m_transformed_bisector;
  Poly_int_2 m_sqf_lcoeff; 
  Poly_int_2 m_sqf_constant_term;
  Poly_lazy_rat_3 m_lazy_rat_transformed_bisector;
  
  //! The list contains the boundary of the $xy$-monotone surface.
  // There could be a better (more compact) representation, but this 
  // representation currently enables us a lot of flexibility.
  // The current interface of Envelope_3 uses CGAL::Object, but we do not
  // need it because there are no isolated points in the boundary.
  Boundary_container m_boundary;

public:
  // We need the Bisector to be default constructable so it can be use as
  // part of the key for the trisector cache.
  // Bisector() {};
  
  Bisector(const Line_3& line, const Poly_int_3 &bisector, 
      const Poly_int_3 &transformed_bisector)
    : m_line(line), m_bisector(bisector), 
      m_transformed_bisector(transformed_bisector){
    CGAL_assertion( bisector == CGAL::canonicalize(bisector));
    CGAL_assertion( transformed_bisector 
        == CGAL::canonicalize(transformed_bisector));
    
    m_sqf_lcoeff        = CGAL::make_square_free(CGAL::leading_coefficient(m_transformed_bisector));
    m_sqf_constant_term = CGAL::make_square_free(CGAL::get_coefficient(m_transformed_bisector,0));
    
    m_lazy_rat_transformed_bisector = 
      typename CGAL::Coercion_traits<Poly_lazy_rat_3,Poly_rat_3>::Cast()(
          CGAL::canonicalize(
              typename CGAL::Coercion_traits<Poly_rat_3,Poly_int_3>::Cast()(transformed_bisector)));
   
    
#if 0 
    std::cout << " line :         " << line << std::endl;
    std::cout << " bisector :     " << bisector << std::endl;
    std::cout << " tbisector      " << transformed_bisector << std::endl;
    std::cout << " degree tbisector x: " << degree(transformed_bisector,0) << std::endl;
    std::cout << " degree tbisector y: " << degree(transformed_bisector,1) << std::endl;
    std::cout << " degree tbisector z: " << degree(transformed_bisector,2) << std::endl;
    std::cout << " lcoeff (h):    " << CGAL::leading_coefficient(transformed_bisector) << std::endl;
    std::cout << " lcoeff (h)sqf: " << make_square_free(CGAL::leading_coefficient(transformed_bisector)) << std::endl;
    std::cout << " tb(0)          " << CGAL::evaluate(transformed_bisector,Poly_int_2(0)) << std::endl;
    std::cout << " tb(0) sqf      " << make_square_free(CGAL::evaluate(transformed_bisector,Poly_int_2(0))) << std::endl;
#endif
  };
    
  const Line_3& line() const { return m_line; }
  const Poly_int_3& bisector() const { return m_bisector; }
  const Poly_int_3& transformed_bisector()       const { return m_transformed_bisector; }
  const Poly_int_2& square_free_lcoeff()         const { return m_sqf_lcoeff; }
  const Poly_int_2& square_free_constant_term()  const { return m_sqf_constant_term; }
  const Poly_lazy_rat_3& lazy_rat_transformed_bisector()  const { return m_lazy_rat_transformed_bisector; }
  

  Boundary_container& boundary() { return m_boundary; }
  const Boundary_container& boundary() const { return m_boundary; }

  mutable boost::optional<int> m_id;

  int id() const { 
    if(!m_id.is_initialized()) {
      m_id = compute_hash_value(line());
    }
    return m_id.get();
  }

  bool operator <(const Self& b) const{
    return (this->id()< b.id());
  }
};


} // namespace VDOL_3
} //namespace CGAL

#endif // CGAL_VDL3_BISECTOR_3_H
