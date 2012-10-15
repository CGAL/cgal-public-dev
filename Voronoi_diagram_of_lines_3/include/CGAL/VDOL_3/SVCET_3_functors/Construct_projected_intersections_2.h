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

/*! \todo Make the function support the partition of the surfaces into
    several $xy$-monotone surfaces. It is OK to return the full
    trisector (without restricting it to the boundary of the 
    $xy$-monotone surface) because the envelope code inserts the curves
    only to the current faces, which is already bounded by the boudnary
    curves. It will be better to return just the correct curves
    (Less "overlaping" curves).            
    However, we must exclude curves on the wrong plane since these are 
    considered as true intersections by the envelope code.
  */

//! \todo check if we can propogate negative curve arcs information
// through singularities using the curve analysis.

#ifndef CGAL_VDL3_SVCET_33_CONSTRUCT_PROJECTED_INTERSECTIONS_2
#define CGAL_VDL3_SVCET_33_CONSTRUCT_PROJECTED_INTERSECTIONS_2

#include <CGAL/config.h>
#include <CGAL/VDOL_3/SVCET_3_functors/SVCET_3_functor_base.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Root_of_traits.h>

namespace CGAL { 
namespace VDOL_3 {
namespace SVCET_3_functors {

template < class SVCET_3 >
class Construct_projected_intersections_2 : public 
Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > {
  
public:
  //! first template parameter
  typedef SVCET_3 Single_voronoi_cell_envelope_traits_3;
  //! the base type 
  typedef Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > Base;
private:
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  typedef typename SVCET_3::Cache_container Cache_container;
  typedef typename Cache_container::Trisector_cache Trisector_cache;
  typedef typename Cache_container::Trisector_cache_key Trisector_cache_key;
  typedef typename Cache_container::Trisector_cache_data Trisector_cache_data;
  typedef std::list<X_monotone_curve_2> XCurve_list;
  typedef typename XCurve_list::iterator XCurve_iterator;
  typedef typename SVCET_3::Algebraic_kernel_d_2::Construct_curve_2 Construct_curve_2;

public: 
  Construct_projected_intersections_2(const SVCET_3 *traits) : Base(traits){}
public:
  template <class OutputIterator>
  OutputIterator
  operator()(
      const Xy_monotone_surface_3& s1,
      const Xy_monotone_surface_3& s2,
      OutputIterator o) const 
  {
    return compute_chached_intersection(s1,s2,o);    
  }  
public:
  template <class OutputIterator>
  OutputIterator
  compute_chached_intersection(
      const Xy_monotone_surface_3& s1,
      const Xy_monotone_surface_3& s2,
      OutputIterator o) const
  {
    if (s1.line() == s2.line())
      // TODO
      // we can't use cache for this, since key is currently not sufficient
      return equal_lines_case(s1,s2,o);

    
    // reoder for proper caching 
    if (s1.bisector() < s2.bisector())
      return compute_chached_intersection(s2, s1, o); 
    
    // key should be unique 
    CGAL_precondition(s1.bisector() == CGAL::canonicalize(s1.bisector()));
    CGAL_precondition(s2.bisector() == CGAL::canonicalize(s2.bisector()));
    CGAL_precondition(s1.bisector() >= s2.bisector()); 
    
    Trisector_cache_key key(this->svcet_3()->state(),s1.bisector(), s2.bisector()); 
    
    bool init_data = (this->svcet_3()->trisector_cache().count(key) == 0); 
    Trisector_cache_data& result = this->svcet_3()->trisector_cache()[key];
    CGAL_assertion(this->svcet_3()->trisector_cache().count(key) == 1 );
    if(init_data){
      compute_intersection(s1,s2,std::back_inserter(result));
    }
    
#if 0
    if(init_data){
    std::cout <<" result.size(): " << result.size() << std::endl; 
      for(typename Trisector_cache_data::iterator it = result.begin();
          it != result.end();it++){
        std::pair<X_monotone_curve_2, Multiplicity> arc_mult_pair;
        CGAL_assertion(CGAL::assign(arc_mult_pair,*it));
        CGAL::assign(arc_mult_pair,*it);
        
        std::cout << arc_mult_pair.first.curve().polynomial_2() << "\n"
                  << arc_mult_pair.first<<"\n"<<std::endl;
      }
    }
#endif 
    return std::copy(result.begin(),result.end(),o);
  }
    
  template <class OutputIterator>
  OutputIterator
  compute_intersection(
      const Xy_monotone_surface_3& s1,
      const Xy_monotone_surface_3& s2,
      OutputIterator o) const
  { 
#if 0 
    // compute projected trisector  
    Poly_int_2 projection = this->svcet_3()->projection(s1,s2);
   
#if 0
    std::cout << "intersection with line intersecting the base line  " << std::endl; 
    std::cout << "s1.bisector:" << s1.bisector() << std::endl;
    std::cout << "s2.bisector:" << s2.bisector() << std::endl;  
    std::cout << "projection degree :" << CGAL::total_degree(projection) << std::endl;
    std::cout << "projection factors:" << std::endl;
    std::list<std::pair<Poly_int_2,int> > sqff; 
    this->svcet_3()->kernel().square_free_factorize_2_object()
      (projection,std::back_inserter(sqff));
    for( typename std::list<std::pair<Poly_int_2,int> >::iterator it = sqff.begin(); 
         it != sqff.end(); it++){
      std::cout << " ( "<<it->first << " )^" << it->second << std::endl; 
    }
#endif

    // make projection square free
    projection = make_square_free(projection); 
#else
    std::list<std::pair<Poly_int_2,int> > sqff_projection; 
    Poly_int_2 projection(1);
    this->svcet_3()->projection(s1,s2,std::back_inserter(sqff_projection));
    for(typename std::list<std::pair<Poly_int_2,int> >::iterator it = sqff_projection.begin();
        it != sqff_projection.end(); it++){
      projection *= it->first; 
    }
    CGAL_postcondition(projection == CGAL::make_square_free(this->svcet_3()->projection(s1,s2))); 
#endif

    // Since the projection is square free, the lcoeff must be square free too
    // otherwise this indicates an event on the discontinuity
    if(!CGAL::is_square_free(CGAL::leading_coefficient(projection))){
      throw Event_at_discontinuity_exception("One curve event on H^*");
    }
    
    // compute common horizontal line at infinity  
#if 0    
    Poly_int_2 common_horizontal_line = CGAL::gcd(
        CGAL::make_square_free(CGAL::leading_coefficient(s1.transformed_bisector())),
        CGAL::make_square_free(CGAL::leading_coefficient(s2.transformed_bisector())));
#else
    Poly_int_2 common_horizontal_line;
    if (s1.square_free_lcoeff() == s2.square_free_lcoeff()) 
      common_horizontal_line = s1.square_free_lcoeff();
    else
      common_horizontal_line = Poly_int_2(1);
#endif
    


    CGAL_assertion(CGAL::degree(common_horizontal_line,0)==0);
    CGAL_assertion(CGAL::total_degree(common_horizontal_line) <= 1);

    // remove common horizontal line from square free projection
    if(CGAL::total_degree(common_horizontal_line)==1){
      projection = CGAL::integral_division(projection,common_horizontal_line);
    }

    // compute common vertical line, two lines intersct base line at the same point
#if 0
    Poly_int_2 common_vertical_line = CGAL::gcd(
        CGAL::make_square_free(CGAL::get_coefficient(s1.transformed_bisector(),0)),
        CGAL::make_square_free(CGAL::get_coefficient(s2.transformed_bisector(),0)));
#else
    Poly_int_2 common_vertical_line;
    if (s1.square_free_constant_term() == s2.square_free_constant_term()) 
      common_vertical_line = s1.square_free_constant_term();
    else
      common_vertical_line = Poly_int_2(1);
#endif

    CGAL_assertion(CGAL::degree(common_vertical_line,1)==0);
    
    // remove common vertical (maybe complex) lines from square free projection 
    if(CGAL::total_degree(common_vertical_line)>=1){
      projection = CGAL::integral_division(projection,common_vertical_line);
    } 
    // in this case the lines must be complex which we do not need to consider.
    if(CGAL::total_degree(common_vertical_line)==2){
      common_vertical_line = Poly_int_2(1);
    } 
    CGAL_assertion(CGAL::total_degree(common_vertical_line) <=2);

    // collect all relavant curves in order to split them into x_monoton_arcs. 
    // relevant curves are, the remaingin part of the projection
    // optional common horizontal and vertical line 
    // all curves that define a boundary of the surfaces 
    std::set<Curve_2,Compare_id<Curve_2> >  curves;

    Construct_curve_2 construct_curve_2 = 
      this->svcet_3()->kernel().construct_curve_2_object();
    curves.insert(construct_curve_2(projection));
    curves.insert(construct_curve_2(common_horizontal_line));
    curves.insert(construct_curve_2(common_vertical_line));

    typedef typename Xy_monotone_surface_3::Boundary_container BContainer;
    typedef typename BContainer::const_iterator BIterator; 
    for(BIterator it = s1.boundary().begin(); it!=s1.boundary().end(); it++){
      curves.insert(it->first.curve());
    } 
    for(BIterator it = s2.boundary().begin(); it!=s2.boundary().end(); it++){
      curves.insert(it->first.curve());
    }
  
    // split curves into non-intersection x-monotone curves 
    XCurve_list xcurves;
    this->svcet_3()->split_curves(curves.begin(),curves.end(),std::back_inserter(xcurves));

    // select the proper xcurves:
    bool report_common_horizontal_line = 
      ((CGAL::total_degree(common_horizontal_line) == 1) && 
          (s1.boundary().size() >= 1)&& (s2.boundary().size() >= 1)) ;
      
    bool report_common_vertical_line = 
      (CGAL::total_degree(common_vertical_line) == 1);
      
    for (XCurve_iterator it = xcurves.begin(); it != xcurves.end(); ++it){
        
      if(report_common_horizontal_line){
        if(it->curve().polynomial_2() == common_horizontal_line ){
          *o++ = make_object(std::make_pair(*it , Multiplicity(0)));
        }
      }
      if(report_common_vertical_line){
        if(it->curve().polynomial_2() == common_vertical_line ){
          *o++ = make_object(std::make_pair(*it , Multiplicity(0)));
        }
      }
      if(CGAL::divides(it->curve().polynomial_2(),projection)){
        if (is_odd_intersection(*it,s1,s2)){
          *o++ = make_object(std::make_pair(*it , Multiplicity(1)));
        }
      }
    }
    return o;
  }

  
private:
  template <typename OutputIterator>
  OutputIterator equal_lines_case(
      const Xy_monotone_surface_3& s1,
      const Xy_monotone_surface_3& s2,
      OutputIterator o) const {
    typedef typename Xy_monotone_surface_3::Boundary_container BContainer;
    typedef typename BContainer::const_iterator BIterator; 

    // Two different Xy_monotone_surfaces of the same surface are intersected.
    // we have to return the common boundary 
    // this also covers the case of intersecting lines .. 
   
    int count = 0; 
    for(BIterator it1 = s1.boundary().begin(); it1!=s1.boundary().end(); it1++){
      for(BIterator it2 = s2.boundary().begin(); it2!=s2.boundary().end(); it2++){
        if(it1->first == it2->first){
          // std::cout << "report:" << it1->first << std::endl;
          count ++ ;
          *o++ = make_object(std::make_pair(it1->first, Multiplicity(0)));
        }
      } 
    }
    CGAL_assertion(count == 1); // TODO report point 
    return o;
  }

  bool is_odd_intersection(
      const X_monotone_curve_2& arc, 
      const Xy_monotone_surface_3& s1, 
      const Xy_monotone_surface_3& s2) const {

#ifdef USE_IS_ODD_INTERSECTION_FILTER 

    // since a curve must at least originate from one side we can report 
    // true if the other side can be excluded by interval arithmetic. 
    
    // get interval polynomials for transformed bisectors 
    typedef CGAL::Interval_nt<true> Interval; 
    typedef CGAL::Coercion_traits<Poly_lazy_rat_3,Interval> CT_3;
    typedef typename CT_3::Type Poly_interval_3; 
    typedef CGAL::Coercion_traits<Poly_lazy_rat_1,Interval> CT_1;
    typedef typename CT_1::Type Poly_interval_1; 
    
    typename CT_3::Cast cast_3; 
    Poly_interval_3 p1_3(cast_3(s1.lazy_rat_transformed_bisector())); 
    Poly_interval_3 p2_3(cast_3(s2.lazy_rat_transformed_bisector()));
    
    //get interval approximation for point on arc 
    typename Algebraic_kernel_d_2::Approximate_relative_x_2
      approx_x = this->svcet_3()->kernel().approximate_relative_x_2_object();
    typename Algebraic_kernel_d_2::Approximate_relative_y_2
      approx_y = this->svcet_3()->kernel().approximate_relative_y_2_object();
    
    Point_2 point_on_arc = 
      this->svcet_3()->construct_interior_vertex_2_object()(arc);
    

    Interval t1_min,t1_max,t2_min,t2_max;
    int prec = 4; 
     do{
       try{
         typename Algebraic_kernel_d_2::Approximate_relative_x_2::result_type 
           xbounds(approx_x(point_on_arc.xy(),prec));
         typename Algebraic_kernel_d_2::Approximate_relative_x_2::result_type
           ybounds(approx_y(point_on_arc.xy(),prec));
         
         Interval x_interval = CGAL::hull(
             Interval(CGAL::to_interval(xbounds.first)),
             Interval(CGAL::to_interval(xbounds.second)));
         Interval y_interval = CGAL::hull(
             Interval(CGAL::to_interval(ybounds.first)),
             Interval(CGAL::to_interval(ybounds.second)));
         
         std::vector<Poly_interval_1> polys;
         polys.push_back(Poly_interval_1((x_interval)));
         polys.push_back(Poly_interval_1((y_interval)));
         polys.push_back(CGAL::shift(Poly_interval_1(Interval(1)),1)); // z^1
         
         Poly_interval_1 p1_1(CGAL::substitute(p1_3,polys.begin(),polys.end()));
         Poly_interval_1 p2_1(CGAL::substitute(p2_3,polys.begin(),polys.end()));
         
         if(CGAL::degree(p1_1)<2) throw 1; // we may add this later  
         if(CGAL::degree(p2_1)<2) throw 1; // we may add this later 
         
         t1_min = make_root_of_2(
             CGAL::get_coefficient(p1_1,2),
             CGAL::get_coefficient(p1_1,1),
             CGAL::get_coefficient(p1_1,0), true); // take the smaller one 
         t1_max = make_root_of_2(
             CGAL::get_coefficient(p1_1,2),
             CGAL::get_coefficient(p1_1,1),
             CGAL::get_coefficient(p1_1,0), false); // take the bigger one 
         
         t2_min = make_root_of_2(
             CGAL::get_coefficient(p2_1,2),
             CGAL::get_coefficient(p2_1,1),
             CGAL::get_coefficient(p2_1,0), true); // take the smaller one 
         t2_max = make_root_of_2(
             CGAL::get_coefficient(p2_1,2),
             CGAL::get_coefficient(p2_1,1),
             CGAL::get_coefficient(p2_1,0), false); // take the bigger one 
         
         // there must be at least one valid intersection on one side 
         CGAL_precondition(CGAL::overlap(t1_min,t2_min) || CGAL::overlap(t1_max,t2_max));
       }catch(...){ 
         t1_min=t2_min=t1_max=t2_max; // makes sure we continue 
       }  
       prec*=2;
     }while(CGAL::overlap(t1_min,t2_min) && CGAL::overlap(t1_max,t2_max) && prec <64);
      
        
      if (CGAL::overlap(t1_min,t2_min) && CGAL::overlap(t1_max,t2_max)){
        std::cout<< " is odd intersection filter failed" << std::endl; 
        return (
            this->svcet_3()->compare_z_at_xy_below_3_object()(arc,s1,s2) !=
            this->svcet_3()->compare_z_at_xy_above_3_object()(arc,s1,s2)); 
      }else{
        bool result; 
        if( CGAL::overlap(t1_min,t2_min) && !CGAL::overlap(t1_max,t2_max))
          result = false; 
        if(!CGAL::overlap(t1_min,t2_min) && CGAL::overlap(t1_max,t2_max))
          result = true; 
        return result; 
      } 
//       if (result != ( 
//               this->svcet_3()->compare_z_at_xy_below_3_object()(arc,s1,s2)!=
//               this->svcet_3()->compare_z_at_xy_above_3_object()(arc,s1,s2))){
//         std::cout<< " BUGG !! " << std::endl; 
//       }
#endif
    
    return (
        this->svcet_3()->compare_z_at_xy_below_3_object()(arc,s1,s2) !=
        this->svcet_3()->compare_z_at_xy_above_3_object()(arc,s1,s2)); 
  
  }  
};

} // namespace SVCET_3_functors
} // namespace VDOL_3 
} // namespace CGAL

#endif // CGAL_VDL3_SVCET_33_CONSTRUCT_PROJECTED_INTERSECTIONS_2




