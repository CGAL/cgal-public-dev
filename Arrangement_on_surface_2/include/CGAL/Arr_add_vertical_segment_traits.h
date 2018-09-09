// Copyright (c) 2011 Tel-Aviv University (Israel), INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/branches/features/Lines_through_segments_3-asafpor/Arrangement_on_surface_2/include/CGAL/Arr_extended_rational_arc_traits_d_1.h $
// $Id: Arr_extended_rational_arc_traits_d_1.h 64848 2011-07-21 07:29:28Z efif $
//
// Author(s)     : Oren Salzman <orenzalz@post.tau.ac.il >
//                 Michael Hemmer <Michael.Hemmer@sophia.inria.fr>

//TODO: somehow use the fact the the x-value is the same in all comparisons

#ifndef CGAL_ARR_ADD_VERTICAL_SEGMENT_TRAITS
#define CGAL_ARR_ADD_VERTICAL_SEGMENT_TRAITS

#include <CGAL/tags.h>
#include <CGAL/assertions.h>
#include <CGAL/Object.h>

#include <vector> 
#include "boost/variant.hpp"
#include <boost/function_output_iterator.hpp>

namespace CGAL {

template <class Traits_>
class Arr_add_vertical_segment_traits
{
public:
  typedef Traits_                                   Traits; 
  typedef Arr_add_vertical_segment_traits<Traits>   Self;
  typedef typename Traits::Point_2                  Point_2;
  typedef typename Traits::Multiplicity             Multiplicity;
public:
  typedef typename Traits::X_monotone_curve_2       XCURVE;
  typedef typename Traits::Curve_2                  CURVE;
public:  
  struct Vertical_segment{
    // store an x coordinate and two optional end points 
    Point_2 min_point; 
    Point_2 max_point;
    bool has_min_point; 
    bool has_max_point; 
    bool is_directed_right; 
    Vertical_segment(const Point_2& p)
      :min_point(p),max_point(p),has_min_point(false),has_max_point(false),is_directed_right(true){}
    Vertical_segment(const Point_2& p, bool upward)
      :min_point(p),max_point(p),has_min_point(upward),has_max_point(!upward),is_directed_right(true){}
    Vertical_segment(const Point_2& min, const Point_2& max)
      :min_point(min),max_point(max),has_min_point(true),has_max_point(true),is_directed_right(true){}
	
    friend std::ostream & operator<<(std::ostream &os, const Vertical_segment& vs) {
      os << "VS " << vs.has_min_point<< " " << vs.has_max_point << " " << vs.min_point << " " << vs.max_point;
      return os;
    }

  };
  
private:
  typedef Vertical_segment VSEGMENT; 
public:  
  typedef boost::variant<CURVE,VSEGMENT>  Curve_2;
  typedef boost::variant<XCURVE,VSEGMENT> X_monotone_curve_2; 
    
  //Category tags:
  typedef Tag_true              Has_left_category;
  typedef Tag_true              Has_merge_category;
  typedef Tag_true              Has_do_intersect_category;

  typedef Arr_open_side_tag     Left_side_category;
  typedef Arr_open_side_tag     Top_side_category;
  typedef Arr_open_side_tag     Bottom_side_category;
  typedef Arr_open_side_tag     Right_side_category;
 
private:
  typedef std::vector<CGAL::Object> Object_vector;


private:
  Traits _traits;
public:
  const Traits& traits() const {return _traits;}
  Traits& traits() {return _traits;}

public:
  //------------
  //Constructors
  //------------

  //---------------------
  // Default constructor.
  Arr_add_vertical_segment_traits(){}
  Arr_add_vertical_segment_traits(const Traits& traits):_traits(traits){}
    
  //------------------------
  //Functor definitions.
  //------------------------
  
  struct Construct_vertical_segment_2{
    typedef Vertical_segment result_type; 
    Vertical_segment operator()(const Point_2& p){
      return Vertical_segment(p);
    }
    Vertical_segment operator()(const Point_2& p, bool upward){
      return Vertical_segment(p,upward);
    }
    Vertical_segment operator()(const Point_2& min, const Point_2& max){
      return Vertical_segment(min,max);
    }
  };
  
  Construct_vertical_segment_2 construct_vertical_segment_2_object(){
    return Construct_vertical_segment_2();
  }
  
  //---------------------------------------------------------------
  //A functor that compares the x-coordinates of two points 
  typedef typename Traits::Compare_x_2 Compare_x_2;
  /*! Obtain a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object () const {
    return _traits.compare_x_2_object(); 
  }
  
  /*! A functor that compares two points lexigoraphically: by x, then by y. */
  typedef typename Traits::Compare_xy_2 Compare_xy_2;
  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const {
    return _traits.compare_xy_2_object(); 
  }
  

  
public:  
  /*! A functor that obtains the left endpoint of a curve. */
  class Construct_min_vertex_2
  {     
    Traits _traits; 
    
    struct Visitor : public boost::static_visitor <const Point_2&> {
      Traits _traits; 
      Visitor(const Traits& traits) :_traits(traits) {}
      const Point_2& operator() (const XCURVE & xc) const
      {
        return _traits.construct_min_vertex_2_object()(xc);
      }
      const Point_2& operator() (const VSEGMENT & vs) const
      {
        CGAL_precondition(vs.has_min_point);
        return (vs.min_point);
      }
    };  // struct Visitor
  public:
    Construct_min_vertex_2(const Traits& traits) : _traits(traits) {}   
    const Point_2& operator() (const X_monotone_curve_2 & cv) const{
      return boost::apply_visitor(Visitor(_traits),cv);
    }
  };  //Construct_min_vertex_2

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const {
    return Construct_min_vertex_2(_traits);
  }


  /*! A functor that obtains the left endpoint of a curve. */
  class Construct_max_vertex_2 
  {    
    const Traits& _traits; 
    
    struct Visitor : public boost::static_visitor <const Point_2&> {
      Traits _traits; 
      Visitor(const Traits& traits) :_traits(traits) {}
      const Point_2& operator() (const XCURVE & xc) const
      {
        return _traits.construct_max_vertex_2_object()(xc);
      }
      const Point_2& operator() (const VSEGMENT & vs) const
      {
        CGAL_precondition(vs.has_max_point);
        return (vs.max_point);
      }
    };  // struct Visitor
    
  public:
    Construct_max_vertex_2(const Traits& traits) : _traits(traits) {}   
    const Point_2& operator() (const X_monotone_curve_2 & cv) const{
      return boost::apply_visitor(Visitor(_traits),cv);
    }
  };  //Construct_max_vertex_2

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const {
    return Construct_max_vertex_2(_traits);
  }
  

  /*! A functor that checks whether a given curve is vertical. */
  struct Is_vertical_2 {
  private:
    struct Visitor : public boost::static_visitor <bool> {
      bool operator() (const XCURVE & xc) const { return false; }
      bool operator() (const VSEGMENT & vs) const { return true; }
    }; // Visitor
  public:
    bool operator() (const X_monotone_curve_2& cv) const {
      return boost::apply_visitor(Visitor(),cv);
    }
  };  //Is_vertical_2

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object () const {
    return Is_vertical_2();
  }

  /*! A functor that compares the y-coordinates of a point and a curve at
   * the point x-coordinate.
   */
  class Compare_y_at_x_2 {
  private:
    struct Visitor : public boost::static_visitor <Comparison_result>
    {
      const Traits& _traits;
      Point_2 _p;
      
      Visitor(const Traits& traits,const Point_2& p) : _traits(traits),_p(p){}
      
      Comparison_result operator() (const XCURVE & xc) const {
        return _traits.compare_y_at_x_2_object() (_p,xc);
      }
      
      Comparison_result operator() (const VSEGMENT & vs) const
      {
        CGAL_precondition(_traits.compare_x_2_object()(_p,vs.max_point) == CGAL::EQUAL);       

        typename Traits::Compare_xy_2 compare_xy_2 = _traits.compare_xy_2_object();

        if(vs.has_min_point){
          if(compare_xy_2(_p,vs.min_point) == CGAL::SMALLER) return CGAL::SMALLER;
        }

        if(vs.has_max_point){
          if(compare_xy_2(_p,vs.max_point) != CGAL::LARGER) return CGAL::LARGER;
        }        
        return CGAL::EQUAL; 
      }
    
    };  //Compare_y_at_x_2_visitor
    

  private:
    const Traits& _traits;    
  
  public:
    Compare_y_at_x_2(const Traits& traits) : _traits(traits) {}
    
    Comparison_result operator() (
        const Point_2& p, 
        const X_monotone_curve_2& cv) const 
    {
      return (boost::apply_visitor(Visitor(_traits,p),cv));
    }
  
  };  //Compare_y_at_x_2

  //*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return Compare_y_at_x_2(_traits);
  }

  class Compare_y_at_x_left_2
  {
  private:
    const Traits& _traits;
  public:
    Compare_y_at_x_left_2(const Traits& traits) : _traits(traits) {}
    Comparison_result operator() (
        const X_monotone_curve_2& cv1,
        const X_monotone_curve_2& cv2,
        const Point_2& p) const { 
      std::cout << "Compare_y_at_x_left_2" << std::endl; 
      return (boost::apply_visitor(Visitor(_traits,p),cv1,cv2));
    }
  private:
    struct Visitor
      : public boost::static_visitor <Comparison_result>
    {
      Visitor(const Traits& traits , const Point_2& p) 
        : _traits(traits), _p(p){}
      
      Comparison_result operator() ( const XCURVE & cv1, const XCURVE & cv2) const {
        return _traits.compare_y_at_x_left_2_object() (cv1,cv2,_p);
      }     
      Comparison_result operator() ( const VSEGMENT & cv1, const XCURVE & cv2) const {
        return CGAL::SMALLER;
      }
      Comparison_result operator() ( const XCURVE& cv1, const VSEGMENT  & cv2) const{
        return CGAL::LARGER;
      }
      Comparison_result operator() ( const VSEGMENT & cv1, const VSEGMENT & cv2) const {
        return CGAL::EQUAL;
      }
    private:
      const Traits& _traits;
      Point_2 _p;
    };  //Visitor
  }; //Compare_y_at_x_left_2

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const {return Compare_y_at_x_left_2(_traits);}
  /*! A functor that compares compares the y-coordinates of two curves
   * immediately to the left of their intersection point.
   */
  
  /*! A functor that checks whether two points and two curves are identical. */
  class Equal_2
  {
  private:
    const Traits& _traits;
  public:
    /*!
     * Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    Equal_2 (const Traits& traits) :_traits(traits) {}
    bool operator() ( const X_monotone_curve_2& cv1,
        const X_monotone_curve_2& cv2) const
    {
      return (boost::apply_visitor(Visitor(_traits),cv1,cv2));
    }
    
  private:
    class Visitor : public boost::static_visitor <bool>
    {
    private:
      typedef boost::static_visitor <bool> Base;      
      const Traits& _traits;
    public:
      Visitor(const Traits& traits) : _traits(traits), Base() {}
      bool operator() (const XCURVE& cv1, const XCURVE& cv2) const {
        return _traits.equal_2_object()  (cv1,cv2);
      }
      bool operator() (const XCURVE& cv1, const VSEGMENT & vs2) const {
        return false; 
      }
      bool operator() (const VSEGMENT & vs1, const XCURVE& cv2) const {
        return false; 
      }
      bool operator() (const VSEGMENT & vs1, const VSEGMENT & vs2) const      
      {
        if (&vs1 == &vs2) return (true);
        
        typename Traits::Compare_xy_2 compare_xy_2 = _traits.compare_xy_2_object();
        typename Traits::Compare_x_2  compare_x_2  = _traits.compare_x_2_object();
        
        if(vs1.has_min_point != vs2.has_min_point) return false; 
        if(vs1.has_max_point != vs2.has_max_point) return false; 
        
        if(vs1.has_min_point && compare_xy_2(vs1.min_point,vs2.min_point) != CGAL::EQUAL)
          return false; 
        
        if(vs1.has_max_point && compare_xy_2(vs1.max_point,vs2.max_point) != CGAL::EQUAL)
          return false; 
        
        if(vs1.has_min_point && vs1.has_max_point)
          return compare_x_2(vs1.min_point,vs2.min_point) == CGAL::EQUAL;
        
        return true;
      }
    };  //Equal_2_visitor
  };  //Equal_2

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object () const { return Equal_2(_traits); }

  
  /*! A functor that divides a curve into continues (x-monotone) curves. */
  class Make_x_monotone_2
  {
  private:    
    const Traits& _traits;
  public:
    Make_x_monotone_2(const Traits& traits) : _traits(traits) {}
    template<class OutputIterator>
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi) const {    
      return boost::apply_visitor(Visitor<OutputIterator>(_traits,oi),cv);
    }
  private:
    template< class OutputIterator> 
    class Visitor : public boost::static_visitor < OutputIterator >
    {     
    private:
      const Traits& _traits;
      mutable OutputIterator _oit; 
    public:
      Visitor(const Traits& traits, OutputIterator oit) : _traits(traits), _oit(oit) {}
      OutputIterator operator() (const CURVE& cv) const      
      {
        std::vector<CGAL::Object> objects; 
        _traits.make_x_monotone_2_object()(cv,std::back_inserter(objects));
        BOOST_FOREACH(const CGAL::Object& obj, objects){
          *_oit++ = Object_recast()(obj);
        }
        return _oit; 
      }
      OutputIterator operator() (const VSEGMENT & cv) const      
      {
        *_oit++ = make_object (X_monotone_curve_2(cv)); 
        return _oit;
      }
    };  //Visitor
  };  //Make_x_monotone_2

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () const {
    return Make_x_monotone_2(_traits);
  }

  /*! A functor that splits a curve at a point. */
  class Split_2
  {
  private:
    const Traits& _traits;
  public:
    Split_2 (const Traits& traits) : _traits(traits) {}
    void operator() ( 
        const X_monotone_curve_2& cv, 
        const Point_2 & p,
        X_monotone_curve_2& c1, 
        X_monotone_curve_2& c2) const {
      boost::apply_visitor(Visitor(_traits,p,c1,c2),cv);
    }
  private:
    class Visitor : public boost::static_visitor <void>
    { 
      const Traits& _traits;
      const Point_2& _p;
      X_monotone_curve_2& _c1;
      X_monotone_curve_2& _c2; 
      
    public:
      Visitor(
          const Traits       &traits,
          const Point_2      &p,
          X_monotone_curve_2 &xc1, 
          X_monotone_curve_2 &xc2) 
        : _traits(traits),_p(p), _c1(xc1), _c2(xc2) {}
      
      void operator() (const XCURVE& cv) const {
        XCURVE c1,c2; 
        _traits.split_2_object() (cv,_p,c1,c2);
        _c1=X_monotone_curve_2(c1);
        _c2=X_monotone_curve_2(c2); 
      }
      
      void operator() (const VSEGMENT & vs) const      
      {
        if(vs.has_min_point){
          _c1 = X_monotone_curve_2(Vertical_segment(vs.min_point,_p));
        }else{
          _c1 = X_monotone_curve_2(Vertical_segment(_p,false)); 
        } 
        
        if(vs.has_max_point){
          _c2 = X_monotone_curve_2(Vertical_segment(_p,vs.max_point));
        }else{
          _c2 = X_monotone_curve_2(Vertical_segment(_p,true)); 
        }
        
      }
    };  //Visitor
  };  //Split_2

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object () const { return Split_2(_traits); }

  /*! A functor that computes intersections between two curves. */
  class Intersect_2
  {
  private:
    const Traits& _traits;    
  public:
    Intersect_2 (const Traits& traits) : _traits(traits) {}
    template<class OutputIterator>
    OutputIterator operator() ( 
        const X_monotone_curve_2& cv1,
        const X_monotone_curve_2& cv2,
        OutputIterator oit) const 
    {
      return boost::apply_visitor(Visitor<OutputIterator>(_traits,oit),cv1,cv2);
    }
    
  private:
    template<class OutputIterator> 
    class Visitor : public boost::static_visitor <OutputIterator>
    {
    private:
      const Traits& _traits;
      mutable OutputIterator _oit; 
    public:
      Visitor(const Traits& traits,OutputIterator oit) 
        : _traits(traits), _oit(oit) {}
      
      //intersection of two XCURVE
      OutputIterator operator() (const XCURVE& cv1,const XCURVE& cv2) const {
        std::vector<CGAL::Object> objects; 
        _traits.intersect_2_object() (cv1,cv2,std::back_inserter(objects));
        BOOST_FOREACH(const CGAL::Object& obj, objects){
          *_oit++ = Object_recast()(obj);
        }
        return _oit; 
      }
      
      //intersection of a XCURVE and a VSEGMENT
      OutputIterator operator() (const VSEGMENT& vs, const XCURVE& cv) const {
        return this->operator()(cv,vs);
      }
      
      //intersection of a XCURVE and a VSEGMENT
      OutputIterator operator() (const XCURVE& cv, const VSEGMENT& vs) const      
      {
        typename Traits::Parameter_space_in_x_2 parameter_space_in_x_2 = _traits.parameter_space_in_x_2_object();
        typename Traits::Parameter_space_in_y_2 parameter_space_in_y_2 = _traits.parameter_space_in_y_2_object();
        typename Traits::Construct_max_vertex_2 max_vertex = _traits.construct_max_vertex_2_object(); 
        typename Traits::Construct_min_vertex_2 min_vertex = _traits.construct_min_vertex_2_object(); 
        typename Traits::Compare_x_2 compare_x_2 = _traits.compare_x_2_object(); 
        typename Traits::Compare_y_at_x_2 compare_y_at_x_2 = _traits.compare_y_at_x_2_object();
        typename Traits::Compare_x_at_limit_2 compare_x_at_limit_2 = _traits.compare_x_at_limit_2_object();
        
        // check x range of curve
        if(parameter_space_in_x_2(cv,ARR_MIN_END) == ARR_INTERIOR){
          if(parameter_space_in_y_2(cv,ARR_MIN_END) == ARR_INTERIOR){
            if(compare_x_2(vs.min_point,min_vertex(cv)) == CGAL::SMALLER) return _oit; 
          }else{
            if(compare_x_at_limit_2(vs.min_point,cv,ARR_MIN_END) == CGAL::SMALLER) return _oit; 
          }
        }        
        if(parameter_space_in_x_2(cv,ARR_MAX_END) == ARR_INTERIOR){
          if(parameter_space_in_y_2(cv,ARR_MAX_END) == ARR_INTERIOR){
            if(compare_x_2(vs.min_point,max_vertex(cv)) == CGAL::LARGER) return _oit; 
          }else{
            if(compare_x_at_limit_2(vs.min_point,cv,ARR_MAX_END) == CGAL::LARGER) return _oit; 
          }
        }
        // check if segment intersects curve 
        if(vs.has_min_point){
          if(compare_y_at_x_2(vs.min_point,cv) == CGAL::LARGER) return _oit; 
        }
        if(vs.has_max_point){
          if(compare_y_at_x_2(vs.max_point,cv) == CGAL::SMALLER) return _oit; 
        }
        // report the projected point         
        *_oit++ = make_object(
            std::make_pair(_traits.project_x_2_object()(vs.min_point,cv),(unsigned int)(1)));
        return _oit; 
        
      }
      //intersection of two VSEGMENT
      OutputIterator operator() (const VSEGMENT& vs1 , const VSEGMENT& vs2) const      
      {
        // segments must have same x coordinate in order to intersect 
        typename Traits::Compare_x_2 compare_x_2 = _traits.compare_x_2_object(); 
        typename Traits::Compare_xy_2 compare_xy_2 = _traits.compare_xy_2_object(); 
        
        if(compare_x_2( vs1.min_point, vs2.max_point) != EQUAL) 
          return _oit; 
        
        
        VSEGMENT result = vs1; 
        if(vs2.has_min_point){
          if (! result.has_min_point || compare_xy_2(result.min_point,vs2.min_point) == SMALLER)
            result.min_point = vs2.min_point;                  
          result.has_min_point = true; 
        }
        
        if(vs2.has_max_point){
          if (! result.has_max_point || compare_xy_2(result.max_point,vs2.max_point) == LARGER)
            result.max_point = vs2.max_point;                  
          result.has_max_point = true; 
        }
        
        if(result.has_max_point && result.has_min_point){
          Comparison_result comp = compare_xy_2(result.min_point,result.max_point);
          if(comp == LARGER) return _oit; // segments did not overlap 
          if(comp == EQUAL){ // segments intersect only in one point
            *_oit++ = make_object(result.min_point); 
            return _oit; 
          }
        }
        
        *_oit++ = make_object(X_monotone_curve_2(result));
        return _oit; 
      }
   
    }; //Visitor
  };  //Intersect_2

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const { return Intersect_2(_traits); }

  /*! A functor that tests whether two curves can be merged. */
  class Are_mergeable_2
  {
  private:
    const Traits& _traits;
  public:
    Are_mergeable_2 (const Traits& traits) : _traits(traits) {}
    bool operator() ( const X_monotone_curve_2& cv1, const X_monotone_curve_2& cv2) const
    {
      return (boost::apply_visitor(Are_mergeable_2_visitor(_traits),cv1,cv2));
    }
  private:
    class Visitor : public boost::static_visitor <bool>
    {     
      const Traits& _traits;
    public:
      Visitor(const Traits& traits) : _traits(traits){}
      
      bool operator() ( const XCURVE& cv1, const XCURVE& cv2) const {
        return  _traits.are_mergeable_2_object() (cv1,cv2);
      }
      bool operator() ( const VSEGMENT& cv1, const XCURVE& cv2) const {
        return false;
      }
      bool operator() ( const XCURVE& cv1, const VSEGMENT& cv2) const {
        return false;
      }
      bool operator() ( const VSEGMENT& vs1, const VSEGMENT& vs2) const
      {
        // overlapping curves are not mergable 
        CGAL_precondition(compare_x_2(vs1.min_point,vs2.min_point)==EQUAL);
        if( vs1.has_min_point && vs2.has_max_point){
          if(equal_2(vs1.min_point,vs2.max_point)) return true; 
        }
        if( vs1.has_max_point && vs2.has_min_point){
          if(equal_2(vs1.max_point,vs2.min_point)) return true;
        }
        return false;
      }
    };  //Visitor
  };  //Are_mergeable_2

  /*! Obtain an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const { return Are_mergeable_2(_traits);}

  /*! A functor that merges two curves into one. */
  class Merge_2
  {
  private:
    const Traits& _traits;
  public:
    Merge_2 (const Traits& traits) : _traits(traits) {}
    void operator() ( 
        const X_monotone_curve_2& cv1,
        const X_monotone_curve_2& cv2,
        X_monotone_curve_2& c) const
    {
      c = boost::apply_visitor(Visitor(_traits),cv1,cv2); 
      return;
    }
  private:
    class Visitor : public boost::static_visitor <X_monotone_curve_2>
    {
      const Traits& _traits;
    public:
      Visitor(const Traits& traits) : _traits(traits){}
      X_monotone_curve_2 operator() (const XCURVE& cv1,const XCURVE& cv2) const
      {
        CGAL_precondition(_traits.are_mergeable_2_object()(cv1,cv2));
        XCURVE _c;
        _traits.merge_2_object()(cv1,cv2,_c);
        return _c;
      }
      X_monotone_curve_2 operator() (const XCURVE& cv1,const VSEGMENT& cv2) const
      {
        CGAL_precondition(false);
        return X_monotone_curve_2();
      }
      X_monotone_curve_2 operator() (const VSEGMENT& cv1,const XCURVE& cv2) const
      {
        CGAL_precondition(false);
        return X_monotone_curve_2();
      }
      X_monotone_curve_2 operator() (const VSEGMENT& vs1,const VSEGMENT& vs2) const
      {
        typename Traits::Compare_xy_2 compare_xy_2 = _traits.compare_xy_2_object();
        VSEGMENT result;
        if( vs1.has_min_point && vs2.has_max_point){
          if(equal_2(vs1.min_point,vs2.max_point)){
            result.has_min_point = vs2.has_min_point;
            result.min_point = vs2.min_point;
            result.has_max_point = vs1.has_max_point;
            result.max_point = vs1.max_point;
            return result; 
          }            
        }
        if( vs1.has_max_point && vs1.has_min_point){
          if(equal_2(vs1.max_point,vs1.min_point)){
            result.has_min_point = vs1.has_min_point;
            result.min_point = vs1.min_point;
            result.has_max_point = vs2.has_max_point;
            result.max_point = vs2.max_point;
            return result; 
          }
        }
        assert(false);
        return result;
      }
    };  //Visitor
  };  //Merge_2

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object () const { return Merge_2(_traits); }

  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2
  {
  private:
    const Traits& _traits;
  public:
    Parameter_space_in_x_2(const Traits& traits) :_traits(traits) {};
    /*! Obtains the parameter space at the end of a line along the x-axis    */
    Arr_parameter_space operator()( 
        const X_monotone_curve_2 & xcv,
        Arr_curve_end ce) const
    {
      return (boost::apply_visitor(Visitor(_traits,ce),xcv));
    }
    
  private:
    class Visitor : public boost::static_visitor <Arr_parameter_space>
    {
      const Traits& _traits;
      Arr_curve_end _ce;
    public:
      Visitor(const Traits& traits,Arr_curve_end ce) 
        :_traits(traits), _ce(ce){}
      Arr_parameter_space operator()(const XCURVE& xcv ) const
      {
        return _traits.parameter_space_in_x_2_object()(xcv,_ce);
      }
      Arr_parameter_space operator()(const VSEGMENT& xcv) const
      {
        return ARR_INTERIOR;
      }

    };  // Visitor
  };  //Parameter_space_in_x_2

  /*! Obtain a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { 
    return Parameter_space_in_x_2(_traits); 
  }
  
  /*! A function object that obtains the parameter space of a geometric
   * entity along the y-axis
   */
  class Parameter_space_in_y_2 
  {
  private:
    const Traits& _traits;
  public:
    Parameter_space_in_y_2 (const Traits& traits) :_traits(traits) {};
    Arr_parameter_space operator()( 
        const X_monotone_curve_2 & xcv,
        Arr_curve_end ce) const
    {
      return (boost::apply_visitor(Visitor(_traits,ce),xcv));
    }
    
    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    //Arr_parameter_space operator()(const Point_2 ) const
    //{
    //  return ARR_INTERIOR;
    //}
      
  private:
    class Visitor : public boost::static_visitor <Arr_parameter_space>
    {
      const Traits& _traits;
      Arr_curve_end _ce;
    public:
      Visitor(const Traits& traits,Arr_curve_end ce) 
        :_traits(traits),_ce(ce){}
      Arr_parameter_space operator()(const XCURVE& xcv ) const
      {
        return _traits.parameter_space_in_y_2_object()(xcv,_ce);
      }
      Arr_parameter_space operator()(const VSEGMENT& vs) const
      {
        if (_ce ==  ARR_MIN_END)
          return (vs.has_min_point)?ARR_INTERIOR:ARR_BOTTOM_BOUNDARY; 
        else  //ce ==  ARR_MAX_END
          return (vs.has_max_point)?ARR_INTERIOR:ARR_TOP_BOUNDARY; 
      }
    }; //Visitor
  };  //Parameter_space_in_y_2
  
  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(_traits); }

  class Compare_y_at_x_right_2
  {
  private:
    const Traits& _traits;
  public:
    Compare_y_at_x_right_2 (const Traits& traits) : _traits(traits) {}
    Comparison_result operator() (
        const X_monotone_curve_2& cv1,
        const X_monotone_curve_2& cv2,
        const Point_2& p) const
    {
      std::cout << "Compare_y_at_x_right_2" << std::endl; 
      return (boost::apply_visitor(Visitor(_traits,p),cv1,cv2));
    }
  private:
    class Visitor
      : public boost::static_visitor <Comparison_result>
    {
      const Traits& _traits;
      Point_2 _p;
    public:
      Visitor(const Traits& traits,const Point_2& p) 
        : _traits(traits), _p(p){}
      Comparison_result operator() (const XCURVE& cv1, const XCURVE& cv2) const      
      {
        return _traits.compare_y_at_x_right_2_object() (cv1,cv2,_p);
      }
      Comparison_result operator() (const XCURVE& cv1, const VSEGMENT & cv2) const      
      {
        return CGAL::SMALLER;
      }
      Comparison_result operator() (const VSEGMENT & cv1, const XCURVE& cv2) const      
      {
        return CGAL::LARGER;
      }
      Comparison_result operator() (const VSEGMENT & cv1, const VSEGMENT & cv2) const      
      {
        return CGAL::EQUAL; //test bug fix        
      }
    }; //Visitor
  }; //Compare_y_at_x_right_2

  /*! Obtain a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return Compare_y_at_x_right_2(_traits);
  }
  class Compare_x_near_limit_2 
  {
  private:
    const Traits& _traits;
  public:
    Compare_x_near_limit_2 (const Traits& traits) : _traits(traits) {}
    Comparison_result operator()( 
        const X_monotone_curve_2& xcv1, 
        const X_monotone_curve_2& xcv2, 
        Arr_curve_end ce) const
    {
      std::cout << "Compare_x_near_limit_2" << std::endl; 
      return (boost::apply_visitor(Visitor(_traits,ce),xcv1,xcv2));
    }
  private:
    class Visitor
      : public boost::static_visitor <Comparison_result>
    {
      const Traits& _traits;
      Arr_curve_end _ce;
    public:
      Visitor(const Traits& traits,Arr_curve_end ce) 
        :_traits(traits), _ce(ce) {}
      
      Comparison_result operator()(const XCURVE& xcv1,const XCURVE& xcv2 ) const
      {
        return _traits.compare_x_near_limit_2_object()(xcv1,xcv2,_ce);
      }
      Comparison_result operator()(const XCURVE& xcv1,const VSEGMENT& xcv2 ) const
      {  
        return (_ce == ARR_MIN_END)? LARGER : SMALLER; 
      }
      Comparison_result operator()(const VSEGMENT& xcv1,const XCURVE& xcv2 ) const
      {
        return (_ce == ARR_MIN_END)? SMALLER : LARGER; 
      }
      Comparison_result operator()(const VSEGMENT& xcv1,const VSEGMENT& xcv2 ) const
      {
        return EQUAL; 
      }
      
    };  //Visitor
  }; //Compare_x_near_limit_2
  
  Compare_x_near_limit_2 compare_x_near_limit_2_object() const 
  {
    return Compare_x_near_limit_2(_traits);
  }

  class Compare_x_at_limit_2 
  {
  private:
    const Traits& _traits;
  public:
    Compare_x_at_limit_2 (const Traits& traits) :_traits(traits) {}
    Comparison_result operator()( 
        const Point_2 & p,
        const X_monotone_curve_2 & xcv,
        Arr_curve_end ce) const
    {
      std::cout << "Compare_x_at_limit_2 with point " << p << std::endl; 
      //std::cout << boost::get<XCURVE>(xcv) << " " << ce << std::endl;
      Comparison_result result(boost::apply_visitor(Visitor_1(_traits,ce,p),xcv)); 
      //std::cout << result << std::endl; 
      return result; 
    }
    Comparison_result operator()( 
        const X_monotone_curve_2 & xcv1,
        Arr_curve_end ce1,
        const X_monotone_curve_2 & xcv2,
        Arr_curve_end ce2) const
    { 
      std::cout << "Compare_x_at_limit_2 two curve ends" << std::endl; 
      //std::cout << boost::get<XCURVE>(xcv1) << " " << ce1 << std::endl;
      //std::cout << boost::get<XCURVE>(xcv2) << " " << ce2 << std::endl;
      Comparison_result result(boost::apply_visitor(Visitor_2(_traits,ce1,ce2),xcv1,xcv2));
      //std::cout << result << std::endl; 
      return result; 
    }
  private:
    class Visitor_1
      : public boost::static_visitor <Comparison_result>
    {
      const Traits& _traits;
      Point_2 _p;
      Arr_curve_end _ce;
    public:
      Visitor_1(const Traits& traits,Arr_curve_end ce,const Point_2 & p) 
        :_traits(traits), _ce(ce), _p(p){}
      Comparison_result operator()(const XCURVE& xcv ) const
      {
        return _traits.compare_x_at_limit_2_object()(_p,xcv,_ce);
      }
      Comparison_result operator()(const VSEGMENT& xcv) const
      {
        return _traits.compare_x_2_object()(_p,xcv.min_point);
      }
    }; //Visitor_1
        
    class Visitor_2
      : public boost::static_visitor <Comparison_result>
    { 
      const Traits& _traits;
      Arr_curve_end _ce1,_ce2;
    public:
      typedef boost::static_visitor <Comparison_result> Base;
      Visitor_2(const Traits& traits,Arr_curve_end ce1,Arr_curve_end ce2) 
        :_traits(traits), _ce1(ce1), _ce2(ce2), Base() {}
      Comparison_result operator()(const XCURVE& xcv1,const XCURVE& xcv2 ) const
      {
        return _traits.compare_x_at_limit_2_object()(xcv1,_ce1,xcv2,_ce2);
      }
      Comparison_result operator()(const XCURVE& xc1,const VSEGMENT& vs ) const
      {
        // INVERT RESULT DUE TO ORDER OF ARGUEMNTS 
        return - _traits.compare_x_at_limit_2_object()(vs.min_point,xc1,_ce1);
      }
      Comparison_result operator()(const VSEGMENT& vs,const XCURVE& xc2 ) const
      {
        return _traits.compare_x_at_limit_2_object()(vs.min_point,xc2,_ce2);
      }
      Comparison_result operator()(const VSEGMENT& vs1,const VSEGMENT& vs2 ) const
      {
        return _traits.compare_x_2_object() (vs1.min_point,vs2.min_point);
      }     
    };  //Visitor_2
  };  //Compare_x_at_limit_2

  /*! Obtain a Compare_x_at_limit_2 function object */
  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(_traits); }
  /*! A function object that compares the y-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 
  {
  private:
    const Traits& _traits;
  public:
    /*! Compare the y-coordinates of 2 lines at their ends near the boundary
     * of the parameter space at x = +/- oo.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the line end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the lines xcv1 and xcv2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     */
    Compare_y_near_boundary_2 (const Traits& traits) :_traits(traits) {}
    Comparison_result operator()( 
        const X_monotone_curve_2 & xcv1,
        const X_monotone_curve_2 & xcv2,
        Arr_curve_end ce) const
    {
      std::cout << "Compare_y_near_boundary_2" << std::endl; 
      return (boost::apply_visitor(Visitor(_traits,ce),xcv1,xcv2));
    }
  private:
    class Visitor
      : public boost::static_visitor <Comparison_result>
    {
      const Traits&       _traits;
      Arr_curve_end _ce;
    public:
      typedef boost::static_visitor <Comparison_result> Base;
      Visitor(const Traits& traits,Arr_curve_end ce) 
        :_traits(traits), _ce(ce), Base() {}
      Comparison_result operator()(const XCURVE& xcv1,const XCURVE& xcv2 ) const
      {
        return _traits.compare_y_near_boundary_2_object()(xcv1,xcv2,_ce);
      }
      Comparison_result operator()(const XCURVE& xcv1,const VSEGMENT& xcv2 ) const
      {
        //a vertical curve is not on the boundary at x = +-oo
        CGAL_precondition(false);
        return Comparison_result();
      }
      Comparison_result operator()(const VSEGMENT& xcv1,const XCURVE& xcv2 ) const
      {
        //a vertical curve is not on the boundary at x = +-oo
        CGAL_precondition (false);
        return Comparison_result();
      }
      Comparison_result operator()(const VSEGMENT& xcv1,const VSEGMENT& xcv2 ) const
      {
        //a vertical curve is not on the boundary at x = +-oo
        CGAL_precondition (false);
        return Comparison_result();
      }
    };  //Visitor
  };

  /*! Obtain a Compare_y_near_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(_traits); }
  //@}
  
  /// \name Functor definitions for the Boolean set-operation traits.
  //@{
  class Compare_endpoints_xy_2
  {
  public:

    /*!
     * Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param cv The curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv)
    {
      return (boost::apply_visitor(Visitor(),cv));
    }
  private:
    class Visitor
      : public boost::static_visitor <X_monotone_curve_2>
    {
    public:
      Comparison_result operator()(const XCURVE& xcv ) const
      {
        return _traits.construct_opposite_2_object()(xcv);
      }
      Comparison_result operator()(const VSEGMENT& xcv) const
      {
        return (xcv.is_directed_right)?SMALLER:LARGER;
      }
    };  //Compare_endpoints_xy_2_visitor
  }; //Compare_endpoints_xy_2
  
  /*! Obtain a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  {
    return Compare_endpoints_xy_2();
  }

  class Construct_opposite_2
  {
  public:
    X_monotone_curve_2 operator() (const X_monotone_curve_2& cv) const
    {
      return (boost::apply_visitor(Visitor(),cv));
    }
  private:
    class Visitor
      : public boost::static_visitor <X_monotone_curve_2>
    {
    public:
      X_monotone_curve_2 operator()(const XCURVE& xcv ) const
      {
        return _traits.construct_opposite_2_object()(xcv);
      }
      X_monotone_curve_2 operator()(const VSEGMENT& xcv) const
      {
        VSEGMENT result(xcv);
        result.is_directed_right = !result.is_directed_right;
        return result; 
      }
    };  //Visitor
  }; // Construct_opposite_2
  
  Construct_opposite_2 construct_opposite_2_object() const
  {
    return Construct_opposite_2();
  }

  struct Object_recast{
    typedef CGAL::Object result_type;
    CGAL::Object operator()(const CGAL::Object& obj){
      std::cout << "converted obj ";      
      {
        const XCURVE* pxc = object_cast<XCURVE> (&obj);
        if (pxc != NULL){
          std::cout << "is XCurve_2" << std::endl;
          CGAL_precondition(pxc->is_valid());
          return  make_object(X_monotone_curve_2(*pxc)); 
        }
      }{      
        const CURVE* pc = object_cast<CURVE> (&obj);
        if (pc != NULL){
          std::cout << "is Curve_2" << std::endl;
          CGAL_precondition(pc->is_valid());
          return  make_object(Curve_2(*pc)); 
         
        }
      }{      
        const Point_2* pp = object_cast<Point_2> (&obj);
        if (pp != NULL){
          std::cout << "is Point_2" << std::endl;
          return obj;
        }
      }{      
        typedef std::pair<Point_2,Multiplicity> PM;
        const PM* ppm = object_cast<PM> (&obj);
        if (ppm != NULL){
          std::cout << "is std::pair<Point_2,Multiplicity>" << std::endl;
          std::cout << ppm->first << " " << ppm->second << std::endl;
          return obj;
        }
      }
      std::cout << "something else" << std::endl;
      // in case it is a point or std::pair<point,mult>
      return obj; 
    }
  };
};

} //namespace CGAL {



#endif  //CGAL_ARR_RATIONAL_ARC_TRAITS_D_1_H
