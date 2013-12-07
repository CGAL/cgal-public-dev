
// Constructing an arrangement of polycurves.

#include <CGAL/basic.h>
#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl; 
  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <vector>
#include <list>
#include "arr_print.h"

#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>


typedef CGAL::Quotient<CGAL::MP_Float>                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;

//////////////
//line segment traits
//////////////
typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>     Polycurve_segment_traits_2;
typedef CGAL::Arrangement_2<Polycurve_segment_traits_2>   Segment_arrangment_2;
typedef Polycurve_segment_traits_2::Point_2               Segment_point_2;
typedef Polycurve_segment_traits_2::Curve_2               Segment_section_2;

///////////////
//circle segment traits
//////////////
typedef CGAL::Arr_circle_segment_traits_2<Kernel>       Arc_traits_2;
typedef CGAL::Arr_polyline_traits_2<Arc_traits_2>       Polycurve_arc_traits_2;
typedef Arc_traits_2::Point_2                           Arc_point_2;
typedef Arc_traits_2::Curve_2                           Arc_segment_2;
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
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>       Conic_traits_2;
typedef Conic_traits_2::Point_2                                           Conic_point_2;
typedef Conic_traits_2::Curve_2                                           Conic_arc_2;
typedef CGAL::Arr_polyline_traits_2<Conic_traits_2>                       Polycurve_conic_traits_2;
typedef CGAL::Arrangement_2<Polycurve_conic_traits_2>                     Conic_arrangment_2; 

///////////////////
//Bezier Curves
///////////////////
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>    Bezier_traits_2;
typedef CGAL::Arr_polyline_traits_2<Bezier_traits_2>                          Polycurve_bezier_traits_2;
typedef CGAL::Arrangement_2<Polycurve_bezier_traits_2>                        Bezier_arrangment_2;

void Make_linear_polycurve (Segment_arrangment_2* Curve_arr)
{

  Segment_point_2               points1[5];
  points1[0] = Segment_point_2 (0, 0);
  points1[1] = Segment_point_2 (2, 4);
  points1[2] = Segment_point_2 (3, 0);
  points1[3] = Segment_point_2 (4, 4);
  points1[4] = Segment_point_2 (6, 0);
  Segment_section_2            pi1 (&points1[0], &points1[5]);

  std::list<Segment_point_2>    points2;
  points2.push_back (Segment_point_2 (1, 3));
  points2.push_back (Segment_point_2 (0, 2));
  points2.push_back (Segment_point_2 (1, 0));
  points2.push_back (Segment_point_2 (2, 1));
  points2.push_back (Segment_point_2 (3, 0));
  points2.push_back (Segment_point_2 (4, 1));
  points2.push_back (Segment_point_2 (5, 0));
  points2.push_back (Segment_point_2 (6, 2));
  points2.push_back (Segment_point_2 (5, 3));
  points2.push_back (Segment_point_2 (4, 2));
  Segment_section_2            pi2 (points2.begin(), points2.end());

  std::vector<Segment_point_2>  points3 (4);
  points3[0] = Segment_point_2 (0, 2);
  points3[1] = Segment_point_2 (1, 2);
  points3[2] = Segment_point_2 (3, 6);
  points3[3] = Segment_point_2 (5, 2);
  Segment_section_2            pi3 (points3.begin(), points3.end());
  
  insert (*Curve_arr, pi1);
  insert (*Curve_arr, pi2);
  insert (*Curve_arr, pi3);
}

void Make_conic_polycurve(Conic_arrangment_2* Conic_arr)
{

  // Insert a hyperbolic arc, supported by the hyperbola y = 1/x
  // (or: xy - 1 = 0) with the endpoints (1/5, 4) and (2, 1/2).
  // Note that the arc is counterclockwise oriented.
  Conic_point_2       ps1 (Rational(1,4), 4);
  Conic_point_2       pt1 (2, Rational(1,2));
  Conic_arc_2         c1 (0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE, ps1, pt1);

  // Insert a parabolic arc that is supported by a parabola y = -x^2
  // (or: x^2 + y = 0) and whose endpoints are (-sqrt(3), -3) ~ (-1.73, -3)
  // and (sqrt(2), -2) ~ (1.41, -2). Notice that since the x-coordinates 
  // of the endpoints cannot be acccurately represented, we specify them
  // as the intersections of the parabola with the lines y = -3 and y = -2.
  // Note that the arc is clockwise oriented.
  Conic_arc_2   c2 =
    Conic_arc_2 (1, 0, 0, 0, 1, 0,       // The parabola.
                 CGAL::CLOCKWISE,
                 Conic_point_2 (-1.73, -3),    // Approximation of the source.
                 0, 0, 0, 0, 1, 3,       // The line: y = -3.
                 Conic_point_2 (1.41, -2),     // Approximation of the target.
                 0, 0, 0, 0, 1, 2);      // The line: y = -2.
  CGAL_assertion (c2.is_valid());
  
  // Insert the segment (1, 1) -- (0, -3).
  Rat_point_2   ps3 (1, 1);
  Rat_point_2   pt3 (0, -3);
  Conic_arc_2   c3 (Rat_segment_2 (ps3, pt3));
  
   insert (*Conic_arr, c1);
   insert (*Conic_arr, c2);
   insert (*Conic_arr, c3);
}

int main ()
{
 
  Segment_arrangment_2        Linear_segment_arr;
  //Arc_arrangment_2          arc_arr;
  Conic_arrangment_2          Conic_arr;
  //Bezier_arrangment_2         bezier_arr;

  Make_linear_polycurve ( & Linear_segment_arr );
  Make_conic_polycurve  ( & Conic_arr);
  
  print_arrangement (Linear_segment_arr);
  //print_arrangement (Conic_arr);

  
  return 0;
}
#endif
