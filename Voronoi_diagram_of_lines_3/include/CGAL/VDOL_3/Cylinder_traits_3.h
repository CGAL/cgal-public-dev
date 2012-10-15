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
// Author(s): Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//
//

#ifndef CGAL_VDL_CYLINDER_TRAITS_3_H
#define CGAL_VDL_CYLINDER_TRAITS_3_H

#include <CGAL/Transformed_output_iterator.h>

namespace CGAL {
namespace VDOL_3 {

template <class SVCET_3>
class Cylinder_traits_3 {
private: //TYPES 
  typedef Cylinder_traits_3<SVCET_3>            CT_3;
  
  typedef typename SVCET_3::Line_3              Line_3; 
  // the following types are defined in the unbounded plane hence prefix P
  typedef typename SVCET_3::Point_2             PPoint_2;
  typedef typename SVCET_3::X_monotone_curve_2  PX_monotone_curve_2;

public: // TAGS 
  typedef Tag_true                              Has_left_category;
  typedef Tag_true                              Has_merge_category;

  typedef Arr_open_side_tag                     Arr_left_side_category;
  typedef Arr_identified_side_tag               Arr_bottom_side_category;
  typedef Arr_identified_side_tag               Arr_top_side_category;
  typedef Arr_open_side_tag                     Arr_right_side_category;
  
public: // TYPES 

  typedef std::pair<Line_3,Line_3>                        Curve_2; 
  typedef unsigned int                                    Multiplicity;

  class Point_2{
    PPoint_2       m_point; 
    CGAL::Sign     m_state;  
  public:
    const PPoint_2& point() const {return m_point;}
    const Sign&     state() const {return m_state;}
    
    Point_2(){};
    Point_2(const PPoint_2& p, CGAL::Sign state):m_point(p),m_state(state){
      CGAL_assertion(point().location() == p.location());
    };
    Point_2(CGAL::Sign state, const PPoint_2& p):m_point(p),m_state(state){
      CGAL_assertion(point().location() == p.location());
    };

    template <class OutputStream>
    friend OutputStream& operator<< (OutputStream& os, Point_2 p)
    {
      os <<"Cylindrical_point_2: " << p.state() <<"; "<< p.point() << std::endl;
      return os;
    } 
  };
  
  class X_monotone_curve_2 {
    
    PX_monotone_curve_2  m_xcurve; 
    CGAL::Sign           m_state;
  public:
    const PX_monotone_curve_2& xcurve() const {return m_xcurve;}
    const Sign&                 state() const {return m_state;}
    
    X_monotone_curve_2(){};
    X_monotone_curve_2(const PX_monotone_curve_2& xcurve, const CGAL::Sign& state)
      :m_xcurve(xcurve), m_state(state){}; 
    X_monotone_curve_2(const CGAL::Sign& state,const PX_monotone_curve_2& xcurve)
      :m_xcurve(xcurve), m_state(state){};  
    
    template <class OutputStream>
    friend OutputStream& operator<< (OutputStream& os, X_monotone_curve_2 xc)
    {
      os <<"Cylindrical_xcurve_2: " << xc.state() <<"; "<< xc.xcurve() << std::endl;
      return os;
    } 
    
    struct Creator: public  std::unary_function<PX_monotone_curve_2,X_monotone_curve_2>{
      CGAL::Sign m_state; 
      Creator(const CGAL::Sign& state):m_state(state){};
      X_monotone_curve_2 operator()(const PX_monotone_curve_2& pxc) const{
        return X_monotone_curve_2(pxc,m_state);
      }
    };
    struct Object_creator: public std::unary_function<CGAL::Object,CGAL::Object>{
      CGAL::Sign m_state; 
      Object_creator(const CGAL::Sign& state):m_state(state){};
      CGAL::Object operator()(const Object& pxc_obj) const{
        PX_monotone_curve_2 pxc(object_cast<PX_monotone_curve_2>(pxc_obj));
        return CGAL::make_object(X_monotone_curve_2(pxc,m_state));
      }
    };
  };
  
private: // MEMBER 
  const SVCET_3 m_pos_svcet_3; 
  const SVCET_3 m_neg_svcet_3; 
public:
  const SVCET_3* svcet_3(const CGAL::Sign& state = POSITIVE) const { 
    if (state == POSITIVE ) return &m_pos_svcet_3; 
    else return &m_neg_svcet_3; 
  }

  int y_index(const Point_2& p) const {
    switch(p.point().location()){
    case ARR_TOP_BOUNDARY:
      if(p.state() == POSITIVE) return 3;
      else                              return 1;
    case ARR_BOTTOM_BOUNDARY:
      if(p.state() == POSITIVE) return 1;
      else                              return 3; // -1 = 3 ???
    default:
      if(p.state() == POSITIVE) return 2;
      else                              return 0;
    }
  }
  
  int y_index(const X_monotone_curve_2& xc) const {
    if(xc.state()==POSITIVE) return 2;
    else return 0; 
  }

public:  // Constructor 
  Cylinder_traits_3(const Line_3&  line_0 , int seed = 0)
    :m_pos_svcet_3(line_0,seed),m_neg_svcet_3(m_pos_svcet_3.mirror()){
  };

private: 
  class Functor_base{
  protected:
    const CT_3* m_ct_3;
    const CT_3* ct_3() const { return m_ct_3; }
    Functor_base(const CT_3* ct_3):m_ct_3(ct_3){};
  };
  
public:
  class Parameter_space_in_x_2 : public  Functor_base {
    typedef Functor_base Base; 
  public:
    Parameter_space_in_x_2(const CT_3* ct_3):Base(ct_3){};
    Arr_parameter_space operator()(const Point_2& p){
      return this->ct_3()->svcet_3()->
        parameter_space_in_x_2_object()(p.point());
    }
    Arr_parameter_space operator()(
        const X_monotone_curve_2& xc, 
        const Arr_curve_end ind){
      return this->ct_3()->svcet_3()->
        parameter_space_in_x_2_object()(xc.xcurve(),ind);
    }     
  };
  
  class Parameter_space_in_y_2 : public  Functor_base {
    typedef Functor_base Base; 
  public:
    Parameter_space_in_y_2(const CT_3* ct_3):Base(ct_3){};
    Arr_parameter_space operator()(const Point_2& p){  
      if( p.point().location() == ARR_TOP_BOUNDARY && 
          p.state()    == POSITIVE)
        return ARR_TOP_BOUNDARY;
      if( p.point().location() == ARR_BOTTOM_BOUNDARY && 
          p.state()    == NEGATIVE)
        return ARR_BOTTOM_BOUNDARY;
      return ARR_INTERIOR;
    }
    Arr_parameter_space operator()(
        const X_monotone_curve_2& xc, 
        const Arr_curve_end ind){
      Arr_parameter_space ps_y = this->ct_3()->svcet_3()->
        parameter_space_in_y_2_object()(xc.xcurve(),ind);

      CGAL_assertion(ps_y != ARR_LEFT_BOUNDARY);
      CGAL_assertion(ps_y != ARR_RIGHT_BOUNDARY);
      
      if(ps_y == ARR_TOP_BOUNDARY &&  xc.state() == POSITIVE)
        return ARR_TOP_BOUNDARY;
      if(ps_y == ARR_BOTTOM_BOUNDARY && xc.state() == NEGATIVE)
        return ARR_BOTTOM_BOUNDARY;
      return ARR_INTERIOR;
    }     
  };
      
  class Compare_x_2 : public Functor_base {
    typedef Functor_base Base; 
  public:
    Compare_x_2(const CT_3* ct_3):Base(ct_3){};
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const {
      return this->ct_3()->svcet_3()->compare_x_2_object()(p1.point(),p2.point());
    }
  };
 
  class Compare_xy_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Compare_xy_2(const CT_3* ct_3):Base(ct_3){}
    Comparison_result operator()(const Point_2 & p1, const Point_2 & p2) const {
      int y_index_1 = this->ct_3()->y_index(p1);
      int y_index_2 = this->ct_3()->y_index(p2);
      // use ordinary compare_xy_2_object() if ARR_INTERIOR of same plane 
      if(y_index_1 == y_index_2 && y_index_1 % 2 == 0)
        return this->ct_3()->svcet_3(p1.state())->
          compare_xy_2_object()(p1.point(),p2.point());
      
      switch(this->ct_3()->compare_x_2_object()(p1,p2)){
      case SMALLER: return SMALLER;
      case LARGER : return LARGER ;
      }
      // x is EQUAL 
      switch(CGAL::compare(y_index_1, y_index_2)){
      case SMALLER: return SMALLER;
      case LARGER : return LARGER; 
      }
      return EQUAL;
    } 
  };
  
  class Construct_min_vertex_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Construct_min_vertex_2(const CT_3* ct_3):Base(ct_3){}
    Point_2 operator()(const X_monotone_curve_2& c) const {
      return Point_2(
          this->ct_3()->svcet_3()->construct_min_vertex_2_object()(c.xcurve()),
          c.state());
    }
  };

  class Construct_max_vertex_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Construct_max_vertex_2(const CT_3* ct_3):Base(ct_3){}
    Point_2 operator()(const X_monotone_curve_2& c) const{
      return Point_2(
          this->ct_3()->svcet_3()->construct_max_vertex_2_object()(c.xcurve()),
          c.state());
    }
  };
  
  class Is_vertical_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Is_vertical_2(const CT_3* ct_3):Base(ct_3){}
    bool operator()(const X_monotone_curve_2& c) const {
      return this->ct_3()->svcet_3()->is_vertical_2_object()(c.xcurve());
    }
  };

  class Compare_y_at_x_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Compare_y_at_x_2(const CT_3* ct_3):Base(ct_3){}
    Comparison_result operator()(
        const Point_2 & p,
        const X_monotone_curve_2 & xc) const{
      switch(compare( this->ct_3()->y_index(p), this->ct_3()->y_index(xc))){
      case SMALLER: return SMALLER;
      case LARGER : return LARGER; 
      }
      CGAL_assertion(p.point().location() == ARR_INTERIOR);
      return this->ct_3()->svcet_3()->
        compare_y_at_x_2_object()(p.point(),xc.xcurve());
    }
  };  

  class Compare_y_at_x_left_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Compare_y_at_x_left_2(const CT_3* ct_3):Base(ct_3){}
    Comparison_result operator()(
        const X_monotone_curve_2 & xc1,
        const X_monotone_curve_2 & xc2,
        const Point_2 & p ) const{
      // should not lie on top or bottom boundary 
      CGAL_assertion( this->ct_3()->y_index(p) != 3);

      switch(compare( this->ct_3()->y_index(xc1), this->ct_3()->y_index(xc2))){
      case SMALLER: return SMALLER;
      case LARGER : return LARGER; 
      }
      
      // this should  not happen since svcet is supposed to throw
      CGAL_assertion(p.point().location() == ARR_INTERIOR);
      return this->ct_3()->svcet_3()->
        compare_y_at_x_left_2_object()(xc1.xcurve(),xc2.xcurve(),p.point());
    }
  };  

  class Compare_y_at_x_right_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Compare_y_at_x_right_2(const CT_3* ct_3):Base(ct_3){}
    Comparison_result operator()(
        const X_monotone_curve_2 & xc1,
        const X_monotone_curve_2 & xc2,
        const Point_2 & p ) const{
      // should not lie on top or bottom boundary 
      CGAL_assertion( this->ct_3()->y_index(p) != 3);

      switch(compare( this->ct_3()->y_index(xc1), this->ct_3()->y_index(xc2))){
      case SMALLER: return SMALLER;
      case LARGER : return LARGER; 
      }
      
      // this should  not happen since svcet is supposed to throw
      CGAL_assertion(p.point().location() == ARR_INTERIOR);
      return this->ct_3()->svcet_3()->
        compare_y_at_x_right_2_object()(xc1.xcurve(),xc2.xcurve(),p.point());
    }
  };  
  
  class Equal_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Equal_2(const CT_3* ct_3):Base(ct_3){};

    bool operator()(
        const X_monotone_curve_2 & xc1,
        const X_monotone_curve_2 & xc2) const 
    {
      if( this->ct_3()->y_index(xc1) != this->ct_3()->y_index(xc2))  return false; 
      
      return this->ct_3()->svcet_3()->
        equal_2_object()(xc1.xcurve(),xc2.xcurve());
    }

    bool operator()(
        const Point_2 & p1,
        const Point_2 & p2) const 
    {
      // avoids a compare in x if the two points can not be equal anyway. 
      if( this->ct_3()->y_index(p1) != this->ct_3()->y_index(p2))  
        return false; 
      return (this->ct_3()->compare_xy_2_object() (p1,p2) == EQUAL);
    }
  };

  class Compare_x_near_boundary_2: public Functor_base {    
    typedef Functor_base Base; 
  public:
    Compare_x_near_boundary_2(const CT_3* ct_3):Base(ct_3){}
    Comparison_result operator()(
        const Point_2 & p,
        const X_monotone_curve_2 & xc,
        Arr_curve_end ce) const
    {
      Parameter_space_in_x_2 ps_x 
        = this->ct_3()->parameter_space_in_x_2_object();
      Parameter_space_in_y_2 ps_y 
        = this->ct_3()->parameter_space_in_y_2_object();
      CGAL_assertion( ps_x(xc,ce)!= ARR_INTERIOR || ps_y(xc,ce)!=ARR_INTERIOR );
      CGAL_assertion( ps_y(p) == ARR_INTERIOR) ;
      return this->ct_3()->svcet_3()->
        compare_x_near_boundary_2_object()(p.point(),xc.xcurve(),ce);
    }
    Comparison_result operator()(
        const X_monotone_curve_2 & xcv1, Arr_curve_end ce1,
        const X_monotone_curve_2 & xcv2, Arr_curve_end ce2) const
    {
      CGAL_assertion(xcv1.state() == xcv2.state());
      return this->ct_3()->svcet_3()->
        compare_x_near_boundary_2_object()(xcv1.xcurve(),ce1,xcv2.xcurve(),ce2);
    }
  };
  
  class Compare_y_near_boundary_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Compare_y_near_boundary_2(const CT_3* ct_3):Base(ct_3){}
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      switch(compare(this->ct_3()->y_index(xcv1),this->ct_3()->y_index(xcv2))){
      case SMALLER: return SMALLER;
      case LARGER : return LARGER; 
      }
      return this->ct_3()->svcet_3()->
        compare_y_near_boundary_2_object()(xcv1.xcurve(),xcv2.xcurve(),ce);
    }
  };
    
  class Is_on_x_identification_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Is_on_x_identification_2(const CT_3* ct_3):Base(ct_3){}
    bool operator()(const Point_2 & p) const
    {
      retrun (this->ct_3()->y_index(p) == 3);
    }
    bool operator()(const X_monotone_curve_2 & xc) const
    {
      return false;
    }
  };
  
  class Is_on_y_identification_2 : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Is_on_y_identification_2(const CT_3* ct_3):Base(ct_3){}
    bool operator()(const Point_2 & p) const {
      return false;
    }
    bool operator()(const X_monotone_curve_2 & xc) const {
      return false;
    }
  };

  class Compare_x_on_boundary_2: public Functor_base {    
    typedef Functor_base Base; 
  public:
    Compare_x_on_boundary_2(const CT_3* ct_3):Base(ct_3){}
    Comparison_result operator()
      (const Point_2 & p1, const Point_2 & p2) const
    {
      CGAL_assertion(this->ct_3()->y_index(p1) == 3 );
      CGAL_assertion(this->ct_3()->y_index(p2) == 3 );
      return CGAL::compare(p1.point().x(),p2.point().x());
    }
  };
  
  class Compare_y_on_boundary_2: public Functor_base {    
    typedef Functor_base Base; 
  public:
    Compare_y_on_boundary_2(const CT_3* ct_3):Base(ct_3){}
    Comparison_result operator()
      (const Point_2 & p1, const Point_2 & p2) const
    {
      // should  not be called since this boundary is open 
      CGAL_assertion(false);
    }
  };
  
  class Make_x_monotone_2: public Functor_base {    
    typedef Functor_base Base; 
  public:  
    Make_x_monotone_2(const CT_3* ct_3):Base(ct_3){}
    
    template<typename OutputIterator>
    OutputIterator operator()( const Curve_2 & c, OutputIterator oi) const 
    {
      assert(false); // should not be called 
    }
  };
  
  class Split_2: public Functor_base {    
    typedef Functor_base Base; 
  public:
    Split_2(const CT_3* ct_3):Base(ct_3){}
    void operator()(
        const X_monotone_curve_2 & xc, const Point_2 & p,
        X_monotone_curve_2 & xc1, X_monotone_curve_2 & xc2) const
    {
      CGAL_assertion( this->ct_3()->y_index(p) == this->ct_3()->y_index(xc) );
      PX_monotone_curve_2 xc1_, xc2_;
      this->ct_3()->svcet_3()->
        split_2_object()(xc.xcurve(),p.point(),xc1_,xc2_);
      xc1 = X_monotone_curve_2(xc1_,xc.state());
      xc2 = X_monotone_curve_2(xc2_,xc.state());
    }
  };
  
  class Intersect_2  : public Functor_base {    
    typedef Functor_base Base; 
  public:
    Intersect_2(const CT_3* ct_3):Base(ct_3){}
    template<typename OutputIterator>
    OutputIterator operator()(
        const X_monotone_curve_2 & xc1,
        const X_monotone_curve_2 & xc2,
        OutputIterator oi) const{

      // first part deals with curves on the same plane 
      if(xc1.state()==xc2.state()){
        std::list<CGAL::Object> objs;
        this->ct_3()->svcet_3()->
          intersect_2_object()(
              xc1.xcurve(),xc2.xcurve(),std::back_inserter(objs));
        PX_monotone_curve_2 xc;
        PPoint_2 p;
        CGAL_assertion_code(int count = 0);
        for(std::list<CGAL::Object>::iterator it = objs.begin();
            it != objs.end(); it++){
          if(assign(xc,*it)){
            *oi++ = make_object(X_monotone_curve_2(xc,xc1.state()));
            CGAL_assertion_code(count ++);
          }
          if(assign(p,*it)){
            *oi++ = make_object(Point_2(p,xc1.state()));
            CGAL_assertion_code(count ++);
          }
        }
        CGAL_assertion(count == objs.size());
        return oi;
      }
      
      //curves lie on different plane, but curve ends may be equal
      Point_2 p1_min = this->ct_3()->construct_min_vertex_2_object()(xc1);
      Point_2 p1_max = this->ct_3()->construct_max_vertex_2_object()(xc1);
      Point_2 p2_min = this->ct_3()->construct_min_vertex_2_object()(xc2);
      Point_2 p2_max = this->ct_3()->construct_max_vertex_2_object()(xc2);
      
      if(this->ct_3()->equal_2_object()(p1_min,p2_min)) *oi++=make_object(p1_min);
      if(this->ct_3()->equal_2_object()(p1_min,p2_max)) *oi++=make_object(p1_min);
      if(this->ct_3()->equal_2_object()(p1_max,p2_min)) *oi++=make_object(p1_max);
      if(this->ct_3()->equal_2_object()(p1_max,p2_max)) *oi++=make_object(p1_max); 
      return oi++;
    }
  };
 
  class Are_mergeable_2 :public Functor_base {    
    typedef Functor_base Base; 
  public:
    Are_mergeable_2(const CT_3* ct_3):Base(ct_3){}
    bool operator()(
        const X_monotone_curve_2 & xc1,
        const X_monotone_curve_2 & xc2) const
    {
      if( this->ct_3()->y_index(xc1) != this->ct_3()->y_index(xc2)) return false; 
      return this->ct_3()->svcet_3()->
        are_mergeable_2_object()(xc1.xcurve(),xc2.xcurve());
    }
  };
  
  class Merge_2 :public Functor_base {    
   typedef Functor_base Base; 
  public:
    Merge_2(const CT_3* ct_3):Base(ct_3){}
    void operator()(
        const X_monotone_curve_2 & xc1,
        const X_monotone_curve_2 & xc2,
        X_monotone_curve_2 & xc) const
    {
      CGAL_precondition (this->ct_3()->are_mergeable_2_object()(xc1, xc2) == true);
      typename SVCET_3::Merge_2 merge_2(this->ct_3()->svcet_3()->merge_2_object());
      PX_monotone_curve_2 xc1_ = xc1.xcurve();
      PX_monotone_curve_2 xc2_ = xc2.xcurve();
      PX_monotone_curve_2 xc_;
      merge_2(xc1_,xc2_,xc_);
      xc=X_monotone_curve_2(xc_,xc1.state());
    }
  };

#define CT_3_snap_functor(_Name,_get_object)            \
  _Name _get_object() const  { return _Name(this); }    \

  CT_3_snap_functor(Compare_x_2,compare_x_2_object);
  CT_3_snap_functor(Compare_xy_2,compare_xy_2_object);
  CT_3_snap_functor(Construct_min_vertex_2,construct_min_vertex_2_object);
  CT_3_snap_functor(Construct_max_vertex_2,construct_max_vertex_2_object); 
  CT_3_snap_functor(Is_vertical_2,is_vertical_2_object);
  CT_3_snap_functor(Compare_y_at_x_2,compare_y_at_x_2_object);
  CT_3_snap_functor(Compare_y_at_x_left_2,compare_y_at_x_left_2_object);
  CT_3_snap_functor(Compare_y_at_x_right_2,compare_y_at_x_right_2_object);
  CT_3_snap_functor(Equal_2,equal_2_object);
  CT_3_snap_functor(Parameter_space_in_x_2,parameter_space_in_x_2_object);
  CT_3_snap_functor(Parameter_space_in_y_2,parameter_space_in_y_2_object);
  CT_3_snap_functor(Compare_x_near_boundary_2,compare_x_near_boundary_2_object);
  CT_3_snap_functor(Compare_y_near_boundary_2,compare_y_near_boundary_2_object);
  CT_3_snap_functor(Compare_x_on_boundary_2,compare_x_on_boundary_2_object);
  CT_3_snap_functor(Compare_y_on_boundary_2,compare_y_on_boundary_2_object);
  CT_3_snap_functor(Is_on_x_identification_2,is_on_x_identification_2_object);
  CT_3_snap_functor(Is_on_y_identification_2,is_on_y_identification_2_object);
  CT_3_snap_functor(Make_x_monotone_2,make_x_monotone_2_object);
  CT_3_snap_functor(Split_2,split_2_object);
  CT_3_snap_functor(Are_mergeable_2,are_mergeable_2_object);
  CT_3_snap_functor(Merge_2,merge_2_object);
  CT_3_snap_functor(Intersect_2,intersect_2_object);
#undef CT_3_snap_functor


  // Envelope_3 fucntors 
  // =======================================================================
  typedef typename SVCET_3::Surface_3                 Surface_3;
  typedef typename SVCET_3::Xy_monotone_surface_3     PXy_monotone_surface_3;
  struct Xy_monotone_surface_3{
    typedef std::pair<CGAL::Sign,PXy_monotone_surface_3> Component;
    typedef std::vector<Component> Component_container;
    typedef typename Component_container::const_iterator Component_const_iterator; 
    typedef typename Component_container::iterator Component_iterator; 
    std::vector<std::pair<CGAL::Sign,PXy_monotone_surface_3> > comps; 
    
    PXy_monotone_surface_3 xy_surface() const {return comps[0].second;}
  };


  
  class Make_xy_monotone_3:public Functor_base {    
   typedef Functor_base Base; 
  public:
    Make_xy_monotone_3(const CT_3* ct_3):Base(ct_3){}

    template <class OutputIterator>
    OutputIterator operator()
      (const Surface_3& l_s, bool is_lower, OutputIterator oi){
      std::vector<PXy_monotone_surface_3> pos_xy_surfaces;
      this->ct_3()->svcet_3(CGAL::POSITIVE)->make_xy_monotone_3_object()
        (l_s,is_lower,std::back_inserter(pos_xy_surfaces));
      std::vector<PXy_monotone_surface_3> neg_xy_surfaces;
      this->ct_3()->svcet_3(CGAL::NEGATIVE)->make_xy_monotone_3_object()
        (l_s,is_lower,std::back_inserter(neg_xy_surfaces));
      

      // the line is skew to the base line. 
      // the positive or negative plane is divided into two surfaces 
      if(pos_xy_surfaces.size() != neg_xy_surfaces.size()){

        // We can not use this code since the Envelope code currently assumes that 
        // boundary cuves are not equal. The problem is that the Surface wrapps around 
        // the cylinder and ends from both sides at the  infinite line.
        // this results in a assertion since the envelope code tries to insert 
        // nonintersecting curve, though the curves overlap. 

//         Xy_monotone_surface_3 xy_surface;
//         if(pos_xy_surfaces.size()==2){ 
//           xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,pos_xy_surfaces[0]));
//           xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,pos_xy_surfaces[1]));
//           xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,neg_xy_surfaces[0]));
//         }else{
//           CGAL_assertion(neg_xy_surfaces.size() == 2);
//           xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,pos_xy_surfaces[0]));
//           xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,neg_xy_surfaces[0]));
//           xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,neg_xy_surfaces[1]));
//         }
//         return *oi++ = xy_surface;
       
        // solution: split surface into two pieces
        
       
        if(pos_xy_surfaces.size()==2){ 
          {
            Xy_monotone_surface_3 xy_surface;
            xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,pos_xy_surfaces[0]));
            xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,pos_xy_surfaces[1]));
            *oi++ = xy_surface;
          }
          {
            Xy_monotone_surface_3 xy_surface;
            xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,pos_xy_surfaces[0]));
            xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,pos_xy_surfaces[1]));
            *oi++ = xy_surface;
          }
        }else{
          CGAL_assertion(neg_xy_surfaces.size() == 2);
          {
            Xy_monotone_surface_3 xy_surface;
            xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,neg_xy_surfaces[0]));
            xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,neg_xy_surfaces[1]));
            *oi++ = xy_surface;
          }
          {
            Xy_monotone_surface_3 xy_surface;
            xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,neg_xy_surfaces[0]));
            xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,neg_xy_surfaces[1]));
            *oi++ = xy_surface;
          }
        }
        return oi;
      }
     
      // the line is parallel to the base line each plane contains on part of the surface 
      if(pos_xy_surfaces.size()== 1 && neg_xy_surfaces.size() == 1){
        std::cout<< " parallel lines " << std::endl; 
        Xy_monotone_surface_3 xy_surface; 
        xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,pos_xy_surfaces[0]));
        xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,neg_xy_surfaces[0]));
        return *oi++ = xy_surface;
      }
      
      // the line intersects the base line, each plane is subdivided into 4 surfaces 
      // we have to glue them together. 
      // this part strongly depends on the order of in which svcet_3 writes the surfaces 
      CGAL_assertion(pos_xy_surfaces.size() == 4);
      CGAL_assertion(neg_xy_surfaces.size() == 4);
      
      {
        Xy_monotone_surface_3 xy_surface; 
        xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,pos_xy_surfaces[0]));
        xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,neg_xy_surfaces[2]));
        *oi++ = xy_surface;
      }{
        Xy_monotone_surface_3 xy_surface; 
        xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,pos_xy_surfaces[1]));
        xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,neg_xy_surfaces[3]));
        *oi++ = xy_surface;
      }{
        Xy_monotone_surface_3 xy_surface; 
        xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,pos_xy_surfaces[2]));
        xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,neg_xy_surfaces[0]));
        *oi++ = xy_surface;
      }{
        Xy_monotone_surface_3 xy_surface; 
        xy_surface.comps.push_back(make_pair(CGAL::POSITIVE,pos_xy_surfaces[3]));
        xy_surface.comps.push_back(make_pair(CGAL::NEGATIVE,neg_xy_surfaces[1]));
        *oi++ = xy_surface;
      }
      return oi;
    }
  };

   class Construct_projected_boundary_2 :public Functor_base {    
   typedef Functor_base Base; 
  public:
    Construct_projected_boundary_2(const CT_3* ct_3):Base(ct_3){}

    template <class OutputIterator>
    OutputIterator
    operator()(const Xy_monotone_surface_3& s, OutputIterator oi){
      typedef typename PXy_monotone_surface_3::Boundary_container Boundary;
      typedef typename Boundary::const_iterator BIT; 
      
      for(int i = 0; i < s.comps.size(); i++){
        CGAL::Sign state = s.comps[i].first;
        Boundary boundary = s.comps[i].second.boundary();
        for(BIT bit = boundary.begin(); bit != boundary.end(); bit++){
          *oi++= make_object(
              std::make_pair(X_monotone_curve_2(bit->first,state),bit->second));
        }
      }
      return oi;
    }
  };



  

  class Construct_projected_intersections_2 :public Functor_base {    
   typedef Functor_base Base; 
  public:
    Construct_projected_intersections_2(const CT_3* ct_3):Base(ct_3){}
 
    struct Object_creator: public std::unary_function<CGAL::Object,CGAL::Object>{
      CGAL::Sign m_state; 
      Object_creator(const CGAL::Sign& state):m_state(state){};
      CGAL::Object operator()(const Object& o) const
      {
        std::pair<PX_monotone_curve_2,Multiplicity> 
          p(object_cast<std::pair<PX_monotone_curve_2,Multiplicity> >(o));
        std::cout << "create intersection xc: " <<m_state<< " "<< p.first << std::endl;
        return CGAL::make_object(
            std::make_pair(X_monotone_curve_2(p.first,m_state),p.second));
      }      
    };
    
    template <class OutputIterator>
    OutputIterator
    operator()(
        const Xy_monotone_surface_3& s1,
        const Xy_monotone_surface_3& s2,
        OutputIterator o) const 
    {
      typedef typename Xy_monotone_surface_3::Component_const_iterator CIT;
      
      for(CIT cit1 = s1.comps.begin(); cit1 != s1.comps.end(); cit1++){
        for(CIT cit2 = s2.comps.begin(); cit2 != s2.comps.end(); cit2++){
          if(cit1->first == cit2->first){
            CGAL::Sign state = cit1->first;
            this->ct_3()->svcet_3(state)->construct_projected_intersections_2_object()
              (cit1->second,cit2->second, 
                  make_transformed_output_iterator(o,Object_creator(state)));
          }
        }
      }      
      return o;
    }
  };
  
  class Compare_z_at_xy_3:public Functor_base {    
    typedef Functor_base Base; 
  public:
    Compare_z_at_xy_3(const CT_3* ct_3):Base(ct_3){}
  
    Comparison_result operator()(
        const Point_2& p,
        const Xy_monotone_surface_3& s1,
        const Xy_monotone_surface_3& s2) const {
      CGAL::Sign state = p.state(); 
      const SVCET_3* svcet_3 = this->ct_3()->svcet_3(state);
      PPoint_2               pp  = p.point();    
      PXy_monotone_surface_3 ps1 = s1.xy_surface();      
      PXy_monotone_surface_3 ps2 = s2.xy_surface();
      return svcet_3->compare_z_at_xy_3_object()(pp,ps1,ps2);
    }
    Comparison_result operator()(
        const X_monotone_curve_2& xc,
        const Xy_monotone_surface_3& s1,
        const Xy_monotone_surface_3& s2) const {
      CGAL::Sign state = xc.state(); 
      const SVCET_3* svcet_3 = this->ct_3()->svcet_3(state);
      PX_monotone_curve_2    pxc = xc.xcurve(); 
      PXy_monotone_surface_3 ps1 = s1.xy_surface();      
      PXy_monotone_surface_3 ps2 = s2.xy_surface();
      return svcet_3->compare_z_at_xy_3_object()(pxc,ps1,ps2);
    }
    Comparison_result operator()(
        const Xy_monotone_surface_3& s1,
        const Xy_monotone_surface_3& s2) const {
      const SVCET_3* svcet_3 = this->ct_3()->svcet_3(POSITIVE);
      PXy_monotone_surface_3 ps1 = s1.xy_surface();      
      PXy_monotone_surface_3 ps2 = s2.xy_surface();
      Comparison_result result = svcet_3->compare_z_at_xy_3_object()(ps1,ps2);
      CGAL_assertion(this->ct_3()->svcet_3(NEGATIVE)->
          compare_z_at_xy_3_object()(ps1,ps2) == result);
      return result; 
    }
  };
  
  class Compare_z_at_xy_below_3 :public Functor_base {    
   typedef Functor_base Base; 
  public:
    Compare_z_at_xy_below_3(const CT_3* ct_3):Base(ct_3){}
 
    Comparison_result
    operator()(const X_monotone_curve_2& cv,
        const Xy_monotone_surface_3& s1,
        const Xy_monotone_surface_3& s2) const {
      assert(false);
      return Comparison_result(0);
    }
  };

  class Compare_z_at_xy_above_3 :public Functor_base {    
   typedef Functor_base Base; 
  public:
    Compare_z_at_xy_above_3(const CT_3* ct_3):Base(ct_3){}
  
    Comparison_result
    operator()(const X_monotone_curve_2& cv,
        const Xy_monotone_surface_3& s1,
        const Xy_monotone_surface_3& s2) const {
      assert(false);
      return Comparison_result(0);
    }
  };
   
#define CT_3_snap_functor(_Name,_get_object)            \
  _Name _get_object() const  { return _Name(this); }    \
    
  CT_3_snap_functor(Make_xy_monotone_3,make_xy_monotone_3_object);
  CT_3_snap_functor(Construct_projected_boundary_2,construct_projected_boundary_2_object);
  CT_3_snap_functor(Construct_projected_intersections_2,construct_projected_intersections_2_object);
  CT_3_snap_functor(Compare_z_at_xy_above_3,compare_z_at_xy_above_3_object);
  CT_3_snap_functor(Compare_z_at_xy_below_3,compare_z_at_xy_below_3_object);
  CT_3_snap_functor(Compare_z_at_xy_3,compare_z_at_xy_3_object);

#undef CT_3_snap_functor
};




} // VDOL_3
} //namespace CGAL

#endif // CGAL_VDL_CYLINDER_TRAITS_3_H
