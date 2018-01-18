//#define CGAL_SL_VERBOSE 1

#include <CGAL/basic.h>
#include <CGAL/CORE_BigInt.h>                      //NT
#include <CGAL/Algebraic_kernel_d_1.h>             //Algebraic Kernel
#include <CGAL/Arr_rational_function_traits_2.h>   //Traits
#include <CGAL/Arrangement_2.h>                    //Arrangement
#include <CGAL/Sweep_line_2_algorithms.h>
#include <boost/foreach.hpp>
#include <CGAL/Algebraic_kernel_2_1.h>

#include <CGAL/Arr_add_vertical_segment_traits.h> 

#if 1
typedef CORE::BigInt                               NT;
typedef CGAL::Algebraic_kernel_d_1<NT>             AK1; 
#else
typedef CORE::BigRat                               Rational;
typedef CGAL::Lazy_exact_nt<Rational>              FT; 
typedef CGAL::Algebraic_kernel_2_1<FT>             AK1;
#endif
typedef CGAL::Arr_rational_function_traits_2<AK1>  Rational_function_traits_2;
typedef CGAL::Arr_add_vertical_segment_traits<Rational_function_traits_2> Traits_2; 

typedef CGAL::Arrangement_2<Traits_2>              Arrangement_2;


int main(){
  CGAL::set_pretty_mode(std::cout);
  // testing that all types are present: 
  
  // Traits 
  typedef Rational_function_traits_2::Polynomial_1         Polynomial_1;
  typedef Rational_function_traits_2::Algebraic_real_1     Algebraic_real_1;
  typedef Rational_function_traits_2::Coefficient          Coefficient; 
  typedef Rational_function_traits_2::Bound                Bound; 
  typedef Rational_function_traits_2::Algebraic_kernel_d_1 Algebraic_kernel_d_1; 
  typedef Rational_function_traits_2::Construct_curve_2    Construct_curve_2;
  typedef Rational_function_traits_2::Construct_x_monotone_curve_2 
    Construct_x_monotone_curve_2; 
  typedef Rational_function_traits_2::Construct_point_2   Construct_point_2;

  typedef Rational_function_traits_2::X_monotone_curve_2 XCURVE;
  typedef Rational_function_traits_2::Curve_2            CURVE;

  // typedef induced by concept 

  typedef Traits_2::Curve_2 Curve_2;
  typedef Traits_2::X_monotone_curve_2 X_monotone_curve_2; 
  typedef Traits_2::Point_2 Point_2; 
  
  typedef Traits_2::Compare_x_2 Compare_x_2;
  typedef Traits_2::Compare_xy_2 Compare_xy_2;
  typedef Traits_2::Construct_min_vertex_2 Construct_min_vertex_2;
  typedef Traits_2::Construct_max_vertex_2 Construct_max_vertex_2; 
  typedef Traits_2::Is_vertical_2 Is_vertical_2;
  typedef Traits_2::Compare_y_at_x_2 Compare_y_at_x_2;
  typedef Traits_2::Compare_y_at_x_left_2 Compare_y_at_x_left_2;
  typedef Traits_2::Compare_y_at_x_right_2 Compare_y_at_x_right_2;
  typedef Traits_2::Equal_2 Equal_2; 
  typedef Traits_2::Parameter_space_in_x_2 Parameter_space_in_x_2; 
  typedef Traits_2::Parameter_space_in_y_2 Parameter_space_in_y_2;
  typedef Traits_2::Compare_x_at_limit_2 Compare_x_at_limit_2;
  typedef Traits_2::Compare_x_near_limit_2 Compare_x_near_limit_2; 
  typedef Traits_2::Compare_y_near_boundary_2 Compare_y_near_boundary_2; 
  typedef Traits_2::Intersect_2 Intersect_2; 
  typedef Traits_2::Split_2 Split_2; 
  typedef Traits_2::Are_mergeable_2 Are_mergeable_2; 
  typedef Traits_2::Merge_2 Merge_2;
  typedef Traits_2::Make_x_monotone_2 Make_x_monotone_2;
  typedef Traits_2::Construct_vertical_segment_2 Construct_vertical_segment_2;

  // construction traits 
  // default construction   
  {
    Traits_2 traits; 
  }
  
  Rational_function_traits_2 rat_func_traits;
  Traits_2 traits(rat_func_traits);
  Arrangement_2 arr(&traits); 
  
  Rational_function_traits_2::Construct_curve_2 
    construct_curve_2 = rat_func_traits.construct_curve_2_object(); 
  Rational_function_traits_2::Construct_point_2 
    construct_point_2 = rat_func_traits.construct_point_2_object(); 
  Construct_vertical_segment_2 construct_vertical_segment_2 = 
    traits.construct_vertical_segment_2_object(); 
  
  Polynomial_1 P(3);
  Polynomial_1 Q(1);
  Polynomial_1 x = CGAL::shift(Polynomial_1(1),1);
  Algebraic_real_1 left(-10);
  Algebraic_real_1 right(10);

  std::vector<Curve_2> curves; 
  curves.push_back(construct_curve_2(Polynomial_1(2),left,right));
  curves.push_back(construct_curve_2(Polynomial_1(0),left,right));
  curves.push_back(construct_curve_2(Polynomial_1(1),left,right));
  curves.push_back(construct_curve_2(Polynomial_1(2),left,right));
  curves.push_back(construct_curve_2(Polynomial_1(3),left,right));
  curves.push_back(construct_curve_2(Polynomial_1(3),left,right));
  curves.push_back(construct_curve_2(3*x-4,x-1,left,right));
  curves.push_back(construct_curve_2(4*x-5,x-2,left,right));
  curves.push_back(construct_curve_2(5*x-6,x-3,left,right));
  curves.push_back(construct_curve_2(6*x-7,x-4,left,right));
  
  curves.push_back(construct_vertical_segment_2(construct_point_2(2,5))); 
  curves.push_back(construct_vertical_segment_2(construct_point_2(1,5),true));
  curves.push_back(construct_vertical_segment_2(construct_point_2(4,5),false));
  curves.push_back(construct_vertical_segment_2(construct_point_2(3,-100),construct_point_2(3,100)));    
  
  CGAL::insert(arr,curves.begin(),curves.end());
  
#if 0
  {     
    const Traits_2 traits;  
    traits.cleanup_cache(); 
    assert( traits.cache().rat_func_map().size()==0);
    assert( traits.cache().rat_pair_map().size()==0);

    Construct_curve_2 construct_curve_2 = traits.construct_curve_2_object();
    
    Polynomial_1 x = CGAL::shift(Polynomial_1(1),1);    
    {      
      std::vector<Curve_2> curves; 
      std::vector<X_monotone_curve_2> xcurves;
      std::vector<Point_2> points;
      curves.push_back(construct_curve_2(Polynomial_1(1)));
      curves.push_back(construct_curve_2(x-2,x-5));
      curves.push_back(construct_curve_2(5*x-5,x-4));
      curves.push_back(construct_curve_2(x-1,6*x-2));
      
     
      BOOST_FOREACH(const Curve_2& curve, curves){
        assert(CGAL::degree(curve.numerator())>=0);
      }   
      CGAL::compute_subcurves(curves.begin(),curves.end(),std::back_inserter(xcurves),false,traits); 
      BOOST_FOREACH(const X_monotone_curve_2& xcurve, xcurves){
        assert(CGAL::degree(xcurve.numerator())>=0);
      }    
	
      CGAL::compute_intersection_points (curves.begin(),curves.end(),std::back_inserter(points),false,traits);
      BOOST_FOREACH(const Point_2& point, points){
        assert(CGAL::degree(point.numerator())>=0);
       }  
    }  
    traits.cleanup_cache(); 
    assert( traits.cache().rat_func_map().size()==0);
    assert( traits.cache().rat_pair_map().size()==0);
    {      
      std::vector<Curve_2> curves; 
      std::vector<X_monotone_curve_2> xcurves;
      std::vector<Point_2> points;
      
      curves.push_back(construct_curve_2(Polynomial_1(1)));
      curves.push_back(construct_curve_2(x-2,x-5));
      curves.push_back(construct_curve_2(5*x-5,x-4));
      curves.push_back(construct_curve_2(x-1,6*x-2));
      
      traits.cleanup_cache();
      BOOST_FOREACH(const Curve_2& curve, curves){
        assert(CGAL::degree(curve.numerator())>=0);
      }   
      CGAL::compute_subcurves(curves.begin(),curves.end(),std::back_inserter(xcurves),false,traits); 
      BOOST_FOREACH(const X_monotone_curve_2& xcurve, xcurves){
        assert(CGAL::degree(xcurve.numerator())>=0);
      }    
      CGAL::compute_intersection_points (curves.begin(),curves.end(),std::back_inserter(points),false,traits);
      BOOST_FOREACH(const Point_2& point, points){
        assert(CGAL::degree(point.numerator())>=0);
      }  
      std::sort(points.begin(),points.end(), traits.compare_xy_2_object());
      std::sort(points.begin(),points.end(), traits.compare_x_2_object());
    }
  }
#endif 

  return 0;
}
