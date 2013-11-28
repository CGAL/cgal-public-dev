
//! \file examples/Arrangement_on_surface_2/polylines.cpp
// Constructing an arrangement of polylines.

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>
#include "arr_print.h"

#include <CGAL/Arr_polycurve_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>



typedef CGAL::Quotient<CGAL::MP_Float>                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;




//////////////
//line segment traits
//////////////
typedef CGAL::Arr_polycurve_traits_2<Segment_traits_2>    Polycurve_segment_traits_2;
typedef CGAL::Arrangement_2<Polycurve_segment_traits_2>   Segment_arrangment_2;
//typedef Polycurve_segment_traits_2::Point_2               Curve_point_2;
//typedef Curve_traits_2::Curve_2                         Polycurve_2;



///////////////
//circle segment traits
//////////////
typedef CGAL::Arr_circle_segment_traits_2<Kernel>       Arc_traits_2;
typedef CGAL::Arr_polycurve_traits_2<Arc_traits_2>      Polycurve_arc_traits_2;
//typedef Arc_traits_2::Point_2                           Arc_point_2;
typedef CGAL::Arrangement_2<Polycurve_arc_traits_2>     Arc_arrangment_2;



////////////////////
//conic traits
////////////////////
typedef CGAL::CORE_algebraic_number_traits                                Nt_traits;
typedef Nt_traits::Rational                                               Rational;
typedef Nt_traits::Algebraic                                              Algebraic;
typedef CGAL::Cartesian<Rational>                                         Rat_kernel;
typedef Rat_kernel::Point_2                                               Rat_point_2;
typedef Rat_kernel::Segment_2                                             Rat_segment_2;
typedef Rat_kernel::Circle_2                                              Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                                        Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>       Conic_Traits_2;
typedef Conic_Traits_2::Point_2                                           Conic_Point_2;
typedef Conic_Traits_2::Curve_2                                           Conic_arc_2;
typedef CGAL::Arr_polycurve_traits_2<Conic_Traits_2>                      Polycurve_conic_traits_2;
typedef CGAL::Arrangement_2<Polycurve_conic_traits_2>                     Conic_arrangment_2; 



///////////////////
//Bezier Curves
///////////////////
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>    Bezier_traits_2;
typedef CGAL::Arr_polycurve_traits_2<Bezier_traits_2>                         Polycurve_bezier_traits_2;
typedef CGAL::Arrangement_2<Polycurve_bezier_traits_2>                        Bezier_arrangment_2;
// typedef CGAL::Arrangement_with_history_2<Traits_2>        Arr_with_hist_2;
// typedef Arr_with_hist_2::Curve_handle                     Curve_handle;

int main ()
{
  //Arrangement_2         arr;



  Segment_arrangment_2          curve_arr;
  //Arc_arrangment_2          arc_arr;
  Conic_arrangment_2          conic_arr;
  //Bezier_arrangment_2         bezier_arr;


  // Point_2               points1[5];
  // points1[0] = Point_2 (0, 0);
  // points1[1] = Point_2 (2, 4);
  // points1[2] = Point_2 (3, 0);
  // points1[3] = Point_2 (4, 4);
  // points1[4] = Point_2 (6, 0);
  // Polyline_2            pi1 (&points1[0], &points1[5]);

  // std::list<Point_2>    points2;
  // points2.push_back (Point_2 (1, 3));
  // points2.push_back (Point_2 (0, 2));
  // points2.push_back (Point_2 (1, 0));
  // points2.push_back (Point_2 (2, 1));
  // points2.push_back (Point_2 (3, 0));
  // points2.push_back (Point_2 (4, 1));
  // points2.push_back (Point_2 (5, 0));
  // points2.push_back (Point_2 (6, 2));
  // points2.push_back (Point_2 (5, 3));
  // points2.push_back (Point_2 (4, 2));
  // Polyline_2            pi2 (points2.begin(), points2.end());

  // std::vector<Point_2>  points3 (4);
  // points3[0] = Point_2 (0, 2);
  // points3[1] = Point_2 (1, 2);
  // points3[2] = Point_2 (3, 6);
  // points3[3] = Point_2 (5, 2);
  // Polyline_2            pi3 (points3.begin(), points3.end());
  
  // insert (arr, pi1);
  // insert (arr, pi2);
  // insert (arr, pi3);
  
  // print_arrangement (arr);

  // Arr_with_hist_2::Curve_iterator            cit;
  // std::cout << "The arrangement contains " << arr.number_of_curves() << " curves:" << std::endl;
  
  // for (cit = arr.curves_begin(); cit != arr.curves_end(); ++cit)
  //   std::cout << "Curve [" << *cit << "] induces " << arr.number_of_induced_edges(cit) << " edges." << std::endl; 
  
  return 0;
}