// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/branches/features/Arrangement_on_surface_2-landmark_no_construction-ophirset/Algebraic_kernel_d/include/CGAL/Algebraic_kernel_2_1.h $
// $Id: Algebraic_kernel_2_1.h 59002 2010-10-04 11:00:27Z lrineau $
// 
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>    
//                 Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//                 Michael Kerber    <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_KERNEL_2_1_H
#define CGAL_ALGEBRAIC_KERNEL_2_1_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Root_of_traits.h>
#include <boost/foreach.hpp>
namespace CGAL {



template< typename FT  >
class Algebraic_kernel_2_1{
public:
  typedef FT                      Coefficient;
  typedef FT                      Bound;
  typedef Polynomial<FT>          Polynomial_1;
  typedef size_t                  size_type; 
  typedef size_t                  Multiplicity_type; 
private:
  typedef Root_of_traits<FT>                ROT;
  typedef Polynomial_traits_d<Polynomial_1> PT;
public:
  typedef typename ROT::Root_of_2 Algebraic_real_1; 
public:
  
  
  struct Solve_1{

    // also report multi in pair
    template < typename OutputIterator >
    OutputIterator operator()(const Polynomial_1& p, OutputIterator oit){
      std::vector<Algebraic_real_1> roots; 

      CGAL_precondition(!CGAL::is_zero(p));
      if(degree(p)==0) return oit; 
      
      CGAL_precondition( CGAL::degree(p)<=2);
      CGAL::compute_roots_of_2(
          CGAL::get_coefficient(p,2),
          CGAL::get_coefficient(p,1),
          CGAL::get_coefficient(p,0),
          std::back_inserter(roots));
      
      BOOST_FOREACH(const Algebraic_real_1& x, roots){
        *oit++=std::make_pair(x,CGAL::degree(p)-roots.size()+1); 
      }
      return oit; 
    }
  };

  Solve_1 solve_1_object() const {return Solve_1();}

  typedef typename CGAL::Real_embeddable_traits<Algebraic_real_1>::Compare 
  Compare_1; 
  
  struct Approximate_relative_1{
    std::pair<Bound,Bound> operator()(const Algebraic_real_1& x, int) const {
      // CGAL_precondition(false);  TODO
      Bound approx = CGAL::to_double(x);
      return std::make_pair(approx,approx);
    }
  };
  Approximate_relative_1 approximate_relative_1_object() const { 
    return  Approximate_relative_1();
  }
  
  typedef Approximate_relative_1 Approximate_absolute_1;
  Approximate_absolute_1 approximate_absolute_1_object() const { 
    return  Approximate_absolute_1();
  }
  
  struct Sign_at_1{    
    template<class NTX> 
    CGAL::Sign operator()(const Polynomial_1& P, const NTX& x) const {
      typename CGAL::Coercion_traits<Polynomial_1,NTX>::Cast cast; 
      return CGAL::sign_at(cast(P),&x,(&x)+1);
    }
  };
  Sign_at_1 sign_at_1_object() const {return Sign_at_1();}
  
  struct Construct_algebraic_real_1{    
    template <class NTX>
    Algebraic_real_1 operator()(const NTX& x) const {
      return Algebraic_real_1(x);
    }
  };
  Construct_algebraic_real_1 construct_algebraic_real_1_object(){ 
    return Construct_algebraic_real_1();
  }
  
};


} //namespace CGAL



#endif // CGAL_ALGEBRAIC_KERNEL_2_1_H
