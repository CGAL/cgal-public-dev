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


#ifndef CGAL_VDL3_SVCET_33_MAKE_XY_MONOTONE_3
#define CGAL_VDL3_SVCET_33_MAKE_XY_MONOTONE_3

#include <CGAL/config.h>
#include <CGAL/VDOL_3/SVCET_3_functors/SVCET_3_functor_base.h>
#include <CGAL/VDOL_3/svcet_3_state_dependent_functions.h>

namespace CGAL { 
namespace VDOL_3 {
namespace SVCET_3_functors {

template < class SVCET_3 >
class Make_xy_monotone_3 : public 
Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > {
public:
  //! first template parameter
  typedef SVCET_3 Single_voronoi_cell_envelope_traits_3;
  //! the base type 
  typedef Single_voronoi_cell_envelope_traits_3_functor_base< SVCET_3 > Base;
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  
  Make_xy_monotone_3(const SVCET_3 *traits): Base(traits) {};

private:
  template <class OutputIterator>
  OutputIterator 
  intersecting_lines_case( const Surface_3& l_s, OutputIterator o){
    CGAL_precondition(CGAL::do_intersect(l_s,this->svcet_3()->base_line()));
    
    // Construct bisector and transformed_bisector;
    Poly_int_3 bisector = construct_bisector_3 (this->svcet_3(), l_s);
    CGAL_assertion(bisector == CGAL::canonicalize(bisector));
    
    Poly_int_3 transformed_bisector = 
      construct_transformed_bisector_3(this->svcet_3(),l_s);
    CGAL_assertion( transformed_bisector 
        == CGAL::canonicalize(transformed_bisector));
    
    // now we have to construct boundary curves. 
    // there is a horizontal line which represents a boundary at infinity 
    Poly_int_2 horizontal_line = 
      CGAL::make_square_free(CGAL::leading_coefficient(transformed_bisector));
    CGAL_assertion(CGAL::total_degree(horizontal_line)<=2);
    if(CGAL::total_degree(horizontal_line) == 0){
      throw Event_at_discontinuity_exception(
          "Irregular projection direction: lcoeff is constant ");
    }
    // make sure it is horizontal, i.e. degree in x is 0, y is 1; 
    CGAL_assertion(CGAL::degree(horizontal_line,0)==0);
    CGAL_assertion(CGAL::degree(horizontal_line,1)==1);
    // there is a vertical line representing at the intersection of the lines 
    Poly_int_2 vertical_line = 
      CGAL::make_square_free(CGAL::evaluate(transformed_bisector,Poly_int_2(0)));
    CGAL_assertion(CGAL::total_degree(vertical_line)==1);
    // make sure it is vertical, i.e. degree in x is 1, y is 0; 
    CGAL_assertion(CGAL::degree(vertical_line,0)==1);
    CGAL_assertion(CGAL::degree(vertical_line,1)==0);
    
    X_monotone_curve_2 horizontal_curve, vertical_curve;
    construct_x_monotone_curves_2(this->svcet_3(),horizontal_line,&horizontal_curve);
    construct_x_monotone_curves_2(this->svcet_3(),vertical_line,&vertical_curve);
    CGAL::Object point_object; 
    this->svcet_3()->
      intersect_2_object()(horizontal_curve, vertical_curve, &point_object);
    std::pair<Point_2,Multiplicity> pm_pair;
    CGAL_assertion(CGAL::assign(pm_pair,point_object));
    CGAL::assign(pm_pair,point_object);
    X_monotone_curve_2 left_horizontal_curve, right_horizontal_curve;
    X_monotone_curve_2 lower_vertical_curve, upper_vertical_curve;
    this->svcet_3()->split_2_object()(
          horizontal_curve,pm_pair.first,
          left_horizontal_curve,right_horizontal_curve);
    this->svcet_3()->split_2_object()(
          vertical_curve,pm_pair.first,
          lower_vertical_curve,upper_vertical_curve);

    // Note that vertical curves are treated as curves with positive slope 
    // that is, left = ON_POSITIVE_SIDE
    // keep invariant that vertical curve is the second (used in Compare_z_at_xy_3)
    // the order is: 
    //  0|1
    // --+--
    //  2|3
    { // upper left 
      Bisector tmp(l_s, bisector, transformed_bisector);
      tmp.boundary().push_back(std::make_pair(left_horizontal_curve, ON_POSITIVE_SIDE));
      tmp.boundary().push_back(std::make_pair(upper_vertical_curve , ON_POSITIVE_SIDE));
      *o++ = tmp;
    }
    { // upper right 
      Bisector tmp(l_s, bisector, transformed_bisector);
      tmp.boundary().push_back(std::make_pair(right_horizontal_curve, ON_POSITIVE_SIDE));
      tmp.boundary().push_back(std::make_pair(upper_vertical_curve ,  ON_NEGATIVE_SIDE));
      *o++ = tmp;
    }
    { // lower left 
      Bisector tmp(l_s, bisector, transformed_bisector);
      tmp.boundary().push_back(std::make_pair(left_horizontal_curve, ON_NEGATIVE_SIDE));
      tmp.boundary().push_back(std::make_pair(lower_vertical_curve , ON_POSITIVE_SIDE));
      *o++ = tmp;
    }
    { // lower right 
      Bisector tmp(l_s, bisector, transformed_bisector);
      tmp.boundary().push_back(std::make_pair(right_horizontal_curve, ON_NEGATIVE_SIDE));
      tmp.boundary().push_back(std::make_pair(lower_vertical_curve  , ON_NEGATIVE_SIDE));
      *o++ = tmp;
    }
    return o;
  }
public:
  // create xy-monotone surfaces from a general surface
  // return a past-the-end iterator
  template <class OutputIterator>
  OutputIterator operator()
    (const Surface_3& l_s, bool is_lower, OutputIterator o) 
  {  
    // We allways want to compute the lower envelop 
    CGAL_precondition(is_lower == true);

    if(CGAL::do_intersect(l_s,this->svcet_3()->base_line())){
      return intersecting_lines_case(l_s,o);
    }
    
    const Linear_kernel& ker = this->svcet_3()->linear_kernel();

    CGAL_precondition_msg(
        ! CGAL::do_intersect(l_s,this->svcet_3()->base_line()),
        "Lines do intersect, this is not implemented yet !"
    );

      
    // Construct bisector and transformed_bisector;
    Poly_int_3 bisector = 
      construct_bisector_3 (this->svcet_3(), l_s);
    Poly_int_3 transformed_bisector = 
      construct_transformed_bisector_3(this->svcet_3(), l_s);
    
    // TODO: throw if bisectors are degenerate i.e. degree in z is less than in y or x 
    
    int td =  CGAL::total_degree(bisector); 
//     std::cout << this->svcet_3()->v1()<< std::endl;
//     std::cout << this->svcet_3()->v2()<< std::endl;
//     std::cout << this->svcet_3()->v3()<< std::endl;
//     std::cout << bisector << std::endl; 
//     std::cout << CGAL::total_degree(bisector) << std::endl; 
//     std::cout << CGAL::degree(bisector,0) << std::endl; 
//     std::cout << CGAL::degree(bisector,1) << std::endl; 
//     std::cout << CGAL::degree(bisector,2) << std::endl; 
    
    if( td != CGAL::degree(transformed_bisector,2) ){
      throw Event_at_discontinuity_exception("Irregular projection direction");
    }
    
//    int tdt =  CGAL::total_degree(transformed_bisector);       
//    CGAL_assertion( tdt != CGAL::degree(transformed_bisector,2) );
    
    // In the case that the lines are parallel, we get that the bisector 
    // surface is define a half-plane on each of the planes.
    
    // in this case we report the intersection of the second surface 
    // with the infinit plane. The curve is the leading coefficient 
    // with respect to r, the outer most variable of the transformed 
    // bisector.
    typename PT_int_3::Leading_coefficient lcoeff_3; 
    typename PT_int_2::Make_square_free make_square_free_2; 

    // some vectors we will need on both case
    // a point on the base line 
    Vector_3 p_base = 
      ker.construct_vector_3_object()(CGAL::ORIGIN,
          ker.construct_point_on_3_object() (this->svcet_3()->base_line(), 0));
    // a point on the other line 
    Vector_3 p_s =   
      ker.construct_vector_3_object()(CGAL::ORIGIN,
          ker.construct_point_on_3_object()(l_s, 0));
    // a vector from base line to the other line 
    Vector_3 diff  = ker.construct_difference_of_vectors_3_object()(p_s,p_base);
    
    // typename PT_int_3::Get_coefficient get_coefficient_3; 
    // if the second surface is a quadric
    if (CGAL::total_degree(bisector) == 2){
      // In the general case, the bisector is a hyperbolic paraboloid.
      // The paraboloid is not define on the whole parameter space. There
      // is one line in the parameter space (sometimes referred to as the 
      // "infinite line" where the bisector is not defined.
      // If, this line is on the current plane, the surface splits into two
      // $xy$-monotone surfaces whose boundary is that line.
      
      Poly_int_2 squared_line = lcoeff_3(transformed_bisector);
      if(CGAL::total_degree(squared_line) == 0){
        // The line equals one of the discontinuity lines
        throw Event_at_discontinuity_exception(
            "The infinite line equals discontinuity line");
      }  
      CGAL_postcondition(CGAL::total_degree(squared_line)==2);
      Poly_int_2 line = make_square_free_2(squared_line);
      CGAL_postcondition(CGAL::total_degree(line) == 1);
      
      // The line at infinity is parallel to the base line and 
      // on the opposit side of the other line l_1. 
     
      
      Vector_3 d_base = ker.construct_vector_3_object()(this->svcet_3()->base_line());
      Vector_3 d_s  = ker.construct_vector_3_object()(l_s);
      // vector perpendiculat to both lines 
      Vector_3 cross_v  = ker.construct_cross_product_vector_3_object()(d_base,d_s);
      
      
      CGAL::Sign sign_1 = 
        CGAL::sign(ker.compute_scalar_product_3_object()(cross_v,diff));
      CGAL::Sign sign_2 = 
        CGAL::sign(ker.compute_scalar_product_3_object()(cross_v,this->svcet_3()->v3()));
      
      if(sign_1 == sign_2){
        // the diff vector and the v_3 vector are on the same side 
        // with respect to the cross vector
        // that is, the line at infinity is not in direction of v_3
        // that is, it is on the negative plane. 
        // that is, we just return the full plane here. 
        *o++ = Bisector(l_s, bisector, transformed_bisector);
      }else{
        // otherwise the line at infitity cuts a bisector surface 
        // that is, we construct two xy_monotone surface, 
        // though, with the same boundary, i.e. the line. 
        Curve_2 cur =  
          this->svcet_3()->kernel().construct_curve_2_object()(line);
        Obj_list xcurves;
        this->svcet_3()->make_x_monotone_2_object() 
          (cur, std::back_inserter(xcurves));
        
        // This is a line and should be only one $x$-monotone curve.
        CGAL_assertion(std::distance(xcurves.begin(), xcurves.end()) == 1);

        X_monotone_curve_2 c;
        CGAL_assertion(CGAL::assign(c, xcurves.front()));
        CGAL::assign(c, xcurves.front());

        // We construct two sub-surfaces --- above and below the curve.
        
        Bisector above_surf(l_s, bisector, transformed_bisector);
        above_surf.boundary().push_back(std::make_pair(c, ON_POSITIVE_SIDE));
        *o++ = above_surf;

        Bisector bellow_surf(l_s, bisector, transformed_bisector);
        bellow_surf.boundary().push_back(std::make_pair(c, ON_NEGATIVE_SIDE));
        *o++ = bellow_surf;
      }
    }else{

      // In this case the bisector is a plane and the cell is
      // only defined on one half-space. The line is valid for both 
      // positive and negative planes. 
      // The only question is, whether the Xy_monotone_surface_3 is to be defined
      // above or below this horizontal line.  
      
      
      CGAL_assertion(CGAL::total_degree(bisector) == 1);
      
      Poly_int_2 line = lcoeff_3(transformed_bisector);
      if (CGAL::total_degree(line) != 1)
        throw Event_at_discontinuity_exception(
            "Intersection of the bisetor plane with "
            "the infinite plane is on discontinuity");
      
      Curve_2 cur = this->svcet_3()->kernel().construct_curve_2_object()(line);
      Obj_list xcurves;
      this->svcet_3()->make_x_monotone_2_object() 
        (cur, std::back_inserter(xcurves));
      // This is a line and should be only one $x$-monotone curve.
      CGAL_assertion(std::distance(xcurves.begin(), xcurves.end()) == 1);
      
      X_monotone_curve_2 c;
      CGAL_assertion(CGAL::assign(c, xcurves.front()));
      CGAL::assign(c, xcurves.front());
      // This was the line, now we have to determine the side on which the 
      // Xy_monotone_surface_3 is defined on. 
      
      // Sign is positive if v2 and vector to other line are in about the same direction 
      CGAL::Sign sign_1 = CGAL::sign(
          ker.compute_scalar_product_3_object()(this->svcet_3()->v2(),diff));
      CGAL_assertion(sign_1 != CGAL::ZERO); // see exception above 
      
      Bisector half_plane(l_s, bisector, transformed_bisector);
      if (sign_1 == CGAL::POSITIVE){  
        // half plane is define in direction of v2, that is above the line  
        half_plane.boundary().push_back(std::make_pair(c, ON_POSITIVE_SIDE));
      } else{
        half_plane.boundary().push_back(std::make_pair(c, ON_NEGATIVE_SIDE));
      }

#ifndef NDEBUG
      // lets  check, we get a point above the line 
      // shoot a ray, and should get a positive zero ....
      // get y coord of horizontal line 
      FT y_coord = FT(- c.curve().polynomial_2()[0][0], c.curve().polynomial_2()[1][0]); 
      y_coord += (sign_1 == CGAL::POSITIVE)?FT(1):-FT(1);
      CGAL_assertion(!is_at_infinity(this->svcet_3(),l_s,std::make_pair(FT(0),y_coord)));
#endif
      *o++ = half_plane;
    }
    return o;
  }
};

} // namespace SVCET_3_functors
} // namespace VDOL_3 
} // namespace CGAL

#endif // CGAL_VDL3_SVCET_33_MAKE_XY_MONOTONE_3
