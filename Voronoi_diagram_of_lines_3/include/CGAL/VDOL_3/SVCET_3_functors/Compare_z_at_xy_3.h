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


#ifndef CGAL_VDL3_SVCET_33_COMPARE_Z_AT_XY_3
#define CGAL_VDL3_SVCET_33_COMPARE_Z_AT_XY_3

#include <CGAL/config.h>
#include <CGAL/VDOL_3/SVCET_3_functors/SVCET_3_functor_base.h>
#include <CGAL/VDOL_3/svcet_3_state_dependent_functions.h>

namespace CGAL { 
namespace VDOL_3 {
namespace SVCET_3_functors {

template < class SVCET_3 >
class Compare_z_at_xy_3 : public 
Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > {

public:
  //! first template parameter
  typedef SVCET_3 Single_voronoi_cell_envelope_traits_3;
  //! the base type 
  typedef Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > Base;
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
  Compare_z_at_xy_3(const SVCET_3 *traits): Base(traits) {};
  
private:
  typedef typename Xy_monotone_surface_3::Boundary_container BC;
  typedef typename BC::const_iterator BIT;

 
public:
  // check which of the surfaces is closer to the envelope at the xy
  // coordinates of point (i.e. lower if computing the lower envelope, or
  // upper if computing the upper envelope)
  // precondition: the surfaces are defined in point

  // In general we will return EQUAL if the point is_on_2 one of 
  // the reported intersection curves of the two surfaces 
  // If this is not the case we choose a sufficiantly close rational 
  // representative to perform the comparision via ray shooting
  // That is, we isolate against all boundaries and the intersection curve 
  
  // There is one exception namely the singular point in the case of 
  // intersecting lines, in general this point will be part of the intersection
  // but the surfaces are only equal if both lines intersect the base line 
  // in this point. 

  // Moreover, it may happen that the code is called with a test for an 
  // undefined surface, i.e. parallel line case. In this case  the defined 
  // surface is cosidered smaller than the undefined surface. 
  // Note that it is not enough to consider the undefined surface as 'at inifinity' 
  // since a point may be on the boundary of a defined surface which is at infinity 
  // In this case the defined surface is to be conisdered as smaller. 
  Comparison_result operator()(
      const Point_2& p,
      const Xy_monotone_surface_3& s1,
      const Xy_monotone_surface_3& s2) const {
    Comparison_result result= compare_at_point(p,s1,s2);
    //std::cout << p <<  std::endl;
    //std::cout << s1.bisector()<< std::endl;
    //std::cout << s2.bisector()<< std::endl;
    //std::cout << result<< std::endl;
    return result; 
  }

  
private:
  bool is_undefined(const Point_2& p, const Xy_monotone_surface_3& s1) const {
    //TODO: use is_undefined in svcet_3_state_dependent_functions.h 
    if(CGAL::total_degree(s1.bisector()) == 2) 
      return false; 
    
    CGAL_assertion(CGAL::total_degree(s1.bisector()) == 1 );
    CGAL_assertion(s1.boundary().size() == 1 );
    
    //if(this->svcet_3()->is_on_2_object()(p,s1.boundary().front().first))
    //  return false;
    
    typename SVCET_3::Compare_y_at_x_2
      compare_y_at_x_2(this->svcet_3()->compare_y_at_x_2_object());
    
    if( compare_y_at_x_2(p,s1.boundary().front().first) == - s1.boundary().front().second )
      return true; 
    return false;  // on boundary or inside 
  }
  
  Comparison_result compare_at_point(
      const Point_2& p,
      const Xy_monotone_surface_3& s1,
      const Xy_monotone_surface_3& s2) const {
    
    CGAL_precondition(p.is_finite());
    //std::cout << " Compare_z_at_xy_3 " << std::endl;
    
    // this can happen due to the fact that the same surface is split up
    // in several xy_monotone patches, and the envelope code may 
    // compare two boundary curves 
    if (s1.line() == s2.line()) return CGAL::EQUAL;

    // Special case of intersecting lines 
    typedef typename SVCET_3::Is_on_2 Is_on_2;
    Is_on_2 is_on_2 = this->svcet_3()->is_on_2_object();
    

    // test whether the point is on the vertical line (back()) only 
    // which means that it is at the intersection point of the line with base line 
    // but not necessarly in the direction of the singular line 
    bool is_singular_s1 = ( s1.boundary().size()==2 && 
        // is_on_2(p,s1.boundary().front().first) && 
        is_on_2(p,s1.boundary().back().first));
    
    bool is_singular_s2 = ( s2.boundary().size()==2 && 
        //is_on_2(p,s2.boundary().front().first) && 
        is_on_2(p,s2.boundary().back().first));
    
    if(is_singular_s1 && is_singular_s2) return EQUAL;
    if(is_singular_s1) return SMALLER; 
    if(is_singular_s2) return LARGER; 
   
    bool is_undefined_s1 = is_undefined(p,s1);
    bool is_undefined_s2 = is_undefined(p,s2);
    
    CGAL_assertion(!(is_undefined_s1 && is_undefined_s2));
    if(is_undefined_s1) return LARGER; 
    if(is_undefined_s2) return SMALLER; 
   
   
    typename SVCET_3::Construct_projected_intersections_2 intersect = 
      this->svcet_3()->construct_projected_intersections_2_object();
    
    
    // collect all intersections 
    // if p lies on the intersection return EQUAL 
    std::vector<X_monotone_curve_2> icurves; 
    typedef typename std::vector<X_monotone_curve_2>::iterator XVIterator;
    // get the x_mon_arcs on the current plane: 
    Obj_list arc_mult_pairs;
    intersect(s1,s2,std::back_inserter(arc_mult_pairs));
    typedef typename Obj_list::iterator OLIterator; 
    for(OLIterator it = arc_mult_pairs.begin(); it !=arc_mult_pairs.end(); it++){
      std::pair<X_monotone_curve_2, Multiplicity> arc_mult_pair;
      CGAL_assertion(CGAL::assign(arc_mult_pair,*it));
      CGAL::assign(arc_mult_pair,*it);
      icurves.push_back(arc_mult_pair.first);
      if(is_on_2(p,arc_mult_pair.first)) return EQUAL; 
    }
   
   

    // NOW comes a regular compare 
    // std::cout <<  "Compare_z_at_xy_3 regular compare" << std::endl;
    // currently we collect all curves and polynomials (guards)
    // thereafter we split all curves into xmonotone arcs 
    // find the arc the point lies on 
    // compute rational point above arc (with respect to all guards)
    // use this point from comparision 

    std::set<Curve_2, Compare_id<Curve_2> > curves;
    std::set<Poly_int_2> guards;
 
    curves.insert(p.curve());
    guards.insert(p.curve().polynomial_2());
    
    for( XVIterator it = icurves.begin(); it != icurves.end(); it++){
      curves.insert(it->curve());
      guards.insert(it->curve().polynomial_2());
    }  
    for(BIT it = s1.boundary().begin();it != s1.boundary().end(); it++){
      curves.insert(it->first.curve());
      guards.insert(it->first.curve().polynomial_2());
    } 
    for(BIT it = s2.boundary().begin();it != s2.boundary().end(); it++){
      curves.insert(it->first.curve());
      guards.insert(it->first.curve().polynomial_2());
    }
    
    std::list<X_monotone_curve_2> xcurves;
    this->svcet_3()->split_curves(curves.begin(),curves.end(),std::back_inserter(xcurves));
    
    typedef  std::pair<FT,FT> FT_pair;
    FT_pair p_above;
    CGAL_assertion_code(FT_pair p_below;);
    
    for(typename std::list<X_monotone_curve_2>::iterator it = xcurves.begin();
        it != xcurves.end(); it++){
      if(this->svcet_3()->is_on_2_object()(p,*it)){
        p_above = rational_point_above_or_below_arc(
            this->svcet_3(),*it,guards.begin(),guards.end(),true);
        CGAL_assertion_code(                                            
            p_below = rational_point_above_or_below_arc(        
                this->svcet_3(),*it,guards.begin(),guards.end(),false););
      }
    }
    Comparison_result res = 
      compare_bisectors_at(this->svcet_3(),s1,s2,p_above);

    return res; 
  }

public:
  // check which of the surfaces is closer to the envelope at the xy
  // coordinates of cv (i.e. lower if computing the lower envelope, or upper
  // if computing the upper envelope)
  // precondition: the surfaces are defined in all points of cv, and the 
  //               answer is the same for each of these points
  Comparison_result operator()(const X_monotone_curve_2& cv,
      const Xy_monotone_surface_3& s1,
      const Xy_monotone_surface_3& s2) const
  {
    // Take a representative from the interior of the point 
    // compare surfaces at that point 
    Point_2 p = this->svcet_3()->construct_interior_vertex_2_object() (cv);
    return (*this)(p, s1, s2);
  }
  // check which of the surfaces is closer to the envelope.
  // (i.e. lower if computing the lower envelope, or upper
  // if computing the upper envelope)
  // precondition: there is no intersections between the surfaces.
  Comparison_result operator()(
      const Xy_monotone_surface_3& s1,
      const Xy_monotone_surface_3& s2) const
  {
    assert(false); // should not be called ?
    return compare_bisectors_at(this->svcet_3(), s1, s2, std::make_pair(FT(0), FT(0)));
  }
};

} // namespace SVCET_3_functors
} // namespace VDOL_3 
} // namespace CGAL

#endif // CGAL_VDL3_SVCET_33_COMPARE_Z_AT_XY_3
