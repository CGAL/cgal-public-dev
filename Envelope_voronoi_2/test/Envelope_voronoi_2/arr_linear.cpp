#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_circle_linear_traits_2.h>
#include <CGAL/Arr_tracing_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <CGAL/Envelope_voronoi_2.h>
#include <CGAL/Envelope_voronoi_traits_2/Moebius_diagram_traits_2.h>
#include <CGAL/Envelope_voronoi_2/Voronoi_diagram_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::FT                                    Number_type;
typedef Kernel::Circle_2                              Circle_2;
typedef Kernel::Segment_2                             Segment_2;
typedef Kernel::Line_2                                Line_2;
typedef Kernel::Ray_2                                 Ray_2;
typedef CGAL::Arr_tracing_traits_2< CGAL::Arr_circle_linear_traits_2<Kernel> >
Traits_2;
typedef Traits_2::CoordNT                             CoordNT;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Curve_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;


typedef CGAL::Moebius_diagram_traits_2<Kernel>        Moebius_diagram_traits_2;
typedef CGAL::Envelope_voronoi_2::Voronoi_diagram_2< Moebius_diagram_traits_2 >
                                                      Voronoi_diagram_2; 
typedef Moebius_diagram_traits_2::Site_2              Site_2;

int main()
{
  std::list<Site_2> sites;
  Voronoi_diagram_2 vor_diag;
  CGAL::voronoi_2(sites.begin(), sites.end(), vor_diag);

  return 0;
}

// int main ()
// {

//   Kernel::Point_2  s3 = Kernel::Point_2 (-2, -2);
//   Kernel::Point_2  t3 = Kernel::Point_2 (2, 2);
//   Ray_2        ray3 = Ray_2 (t3, s3);

//   curves.push_back (Curve_2 (ray3));

//   // Create a line segment with the same supporting line (y = x), but
//   // having one endpoint with irrational coefficients.
//   CoordNT          sqrt_15 = CoordNT (0, 1, 15); // = sqrt(15)
//   Point_2          s4 = Point_2 (3, 3);
//   Point_2          t4 = Point_2 (sqrt_15, sqrt_15);
//   Line_2           line4 = Line_2 (s3, t3);

//   curves.push_back (Curve_2 (line4, s4, t4));

//   // Create a circular arc that correspond to the upper half of the
//   // circle centered at (1,1) with squared radius 3. We create the
//   // circle with clockwise orientation, so the arc is directed from
//   // (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
//   Kernel::Point_2  c5 = Kernel::Point_2 (1, 1);
//   Circle_2         circ5 = Circle_2 (c5, 3, CGAL::CLOCKWISE);
//   CoordNT          one_minus_sqrt_3 = CoordNT (1, -1, 3);
//   CoordNT          one_plus_sqrt_3 = CoordNT (1, 1, 3);
//   Point_2          s5 = Point_2 (one_minus_sqrt_3, CoordNT (1));
//   Point_2          t5 = Point_2 (one_plus_sqrt_3, CoordNT (1));

//   curves.push_back (Curve_2 (circ5, s5, t5));

//   // Create a circular arc of the unit circle, directed clockwise from
//   // (-1/2, sqrt(3)/2) to (1/2, sqrt(3)/2). Note that we orient the
//   // supporting circle accordingly.
//   Kernel::Point_2  c6 = Kernel::Point_2 (0, 0);
//   CoordNT          sqrt_3_div_2 = CoordNT (0, Number_type(1)/Number_type(2), 3);
//   Point_2          s6 = Point_2 (Number_type(-1)/Number_type(2), sqrt_3_div_2);
//   Point_2          t6 = Point_2 (Number_type(1)/Number_type(2), sqrt_3_div_2);

//   curves.push_back (Curve_2 (c6, 1, CGAL::CLOCKWISE, s6, t6));

//   // Create a circular arc defined by two endpoints and a midpoint,
//   // all having rational coordinates. This arc is the upper-right
//   // quarter of a circle centered at the origin with radius 5.
//   Kernel::Point_2  s7 = Kernel::Point_2 (0, 5);
//   Kernel::Point_2  mid7 = Kernel::Point_2 (3, 4);
//   Kernel::Point_2  t7 = Kernel::Point_2 (5, 0);

//   curves.push_back (Curve_2 (s7, mid7, t7));

//   // Construct the arrangement of the curves.
//   Arrangement_2    arr;

//   insert (arr, curves.begin(), curves.end());

//   // Print the size of the arrangement.
//   std::cout << "The arrangement size:" << std::endl
//             << "   V = " << arr.number_of_vertices()
//             << ",  E = " << arr.number_of_edges()
//             << ",  F = " << arr.number_of_faces() << std::endl;

//   return (0);
// }
