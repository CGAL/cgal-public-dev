// Copyright (c) 2005-2009 Tel-Aviv University (Israel).
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

// SVCET_3 = Single_voronoi_cell_envelope_traits_3


#ifndef CGAL_VDL3_SVCET_33_COMPARE_Z_AT_XY_NEAR_ARC_3
#define CGAL_VDL3_SVCET_33_COMPARE_Z_AT_XY_NEAR_ARC_3

#include <CGAL/config.h>
#include <CGAL/VDOL_3/SVCET_3_functors/SVCET_3_functor_base.h>
#include <CGAL/VDOL_3/svcet_3_state_dependent_functions.h>
#include <CGAL/VDOL_3/rational_point_above_or_below_arc.h>

namespace CGAL { 
namespace VDOL_3 {
namespace SVCET_3_functors {

template < class SVCET_3 >
class Compare_z_at_xy_near_arc_3 : public 
Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > {
private:
  bool m_above; 
public:
  const bool& above() const {return m_above;} 
public:
  //! first template parameter
  typedef SVCET_3 Single_voronoi_cell_envelope_traits_3;
  //! the base type 
  typedef Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > Base;
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;

  typedef typename SVCET_3::Cache_container Cache_container;
  
  
  Compare_z_at_xy_near_arc_3(const SVCET_3 *traits, bool above): Base(traits), m_above(above){};
  

  
  // check which of the surfaces is closer to the envelope on the points above/below 
  // the curve cv
  // precondition: the surfaces are defined above/below cv
  //               the choise between s1 and s2 for the envelope is the same
  //               for every point in the infinitesimal region above cv
  //               the surfaces are EQUAL over the curve cv

  
  Comparison_result
  operator()(const X_monotone_curve_2& cv,
      const Xy_monotone_surface_3& s1,
      const Xy_monotone_surface_3& s2) const {

    return compare_at_arc(cv,s1,s2); 
      
//       std::cout << "new cache value " << std::endl;
//       std::cout <<" cv: "<< cv.curve().polynomial_2() << std::endl;
//       std::cout <<" cv: "<< cv << std::endl;
//       std::cout <<" id: "<< cv.id() << std::endl;
//       std::cout <<" state: "<< this->svcet_3()->state() << std::endl;
//       std::cout <<" above: "<< this->above() << std::endl;
//       std::cout <<" s1: "<< s1.bisector() << std::endl;
//       std::cout <<" s2: "<< s2.bisector() << std::endl;
//       std::cout <<" result "<< result << std::endl;
//       std::cout << std::endl;
//       this->svcet_3()->compare_at_arc_cache()
//         .insert(std::make_pair(key, result));
//     }else{
//       std::cout << "old cache value " << std::endl;
//       result = it->second;
//     }

//     if(result != compare_at_arc(cv,s1,s2)){
//       std::cout <<" cv: "<< cv.curve().polynomial_2() << std::endl;
//       std::cout <<" cv: "<< cv << std::endl;
//       std::cout <<" id: "<< cv.id() << std::endl;
//       std::cout <<" state: "<< this->svcet_3()->state() << std::endl;
//       std::cout <<" above: "<< this->above() << std::endl;
//       std::cout <<" s1: "<< s1.bisector() << std::endl;
//       std::cout <<" s2: "<< s2.bisector() << std::endl;
//       std::cout <<" result "<< result << std::endl;
//       std::cout <<" vresult "<< compare_at_arc(cv,s1,s2) << std::endl;
//     }

//     CGAL_assertion(result == compare_at_arc(cv,s1,s2) );
//     return result; 
  }
private:
  Comparison_result compare_at_arc(
      const X_monotone_curve_2& cv,
      const Xy_monotone_surface_3& s1,
      const Xy_monotone_surface_3& s2) const {
      
    // within the arrangment defined by the overlay of the curve supporting 
    // cv and the boundary curves of s1 and s2, we have to compute a rational 
    // in the cell above the arc. 
    // we therefore name the polynomials supporting the boundary curves 
    // as guards, since we are not allowed to cross them. 
    std::set<Poly_int_2> guards; 
    typename Xy_monotone_surface_3::Boundary_container::const_iterator bit; 
    for(bit = s1.boundary().begin();bit != s1.boundary().end(); bit++){
      guards.insert(bit->first.curve().polynomial_2());
    }
    for(bit = s2.boundary().begin();bit != s2.boundary().end(); bit++){
      guards.insert(bit->first.curve().polynomial_2());
    }

    std::list<std::pair<Poly_int_2,int> > sqff; 
    this->svcet_3()->projection(s1,s2,std::back_inserter(sqff));

    typedef typename std::list<std::pair<Poly_int_2,int> >::iterator SQFFIT;    
    // and add it to the guards 
    for( SQFFIT it = sqff.begin(); it != sqff.end(); it++){
      guards.insert(it->first);
    }
    
    std::pair<FT, FT> p_above = 
    rational_point_above_or_below_arc(
        this->svcet_3(), cv, guards.begin(), guards.end(), this->above()); 
    
    Comparison_result res = 
      compare_bisectors_at(this->svcet_3(),s1,s2,p_above);
        
    // in case of two parallel lines to the base line 
    // it may happen that both surfaces appear equal since both roots are negative. 
    CGAL_assertion_code(bool s1_is_parallel = (CGAL::total_degree(s1.bisector())==1));
    CGAL_assertion_code(bool s2_is_parallel = (CGAL::total_degree(s2.bisector())==1));
    CGAL_postcondition( (s1_is_parallel && s2_is_parallel) ||  res != CGAL::EQUAL);
    return res;
  }  
};

} // namespace SVCET_3_functors
} // namespace VDOL_3 
} // namespace CGAL

#endif // CGAL_VDL3_SVCET_33_COMPARE_Z_AT_XY_NEAR_ARC_3
