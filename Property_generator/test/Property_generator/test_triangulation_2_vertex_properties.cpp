#include <iostream>
#include <limits>
#include <cmath>
#include <cassert>
#include "test_triangulation_2_types.h"
#include <CGAL/properties/triangulation_2.h>

using namespace CGAL::Properties::Triangulation_2;

bool approx_eq(double a, double b)
{
  static double epsilon = 1e-3;
  std::cout << a << " " << b << std::endl;
  return std::fabs(a-b) <= epsilon;
}

using namespace std;

int main()
{
  double infinity = std::numeric_limits<double>::infinity(); 

  // Triangulations for testing. We consider three test cases,
  // an empty triangulation, a triangulation with one point, and a 
  // triangulation with five points.
  Delaunay dt_0;
  Delaunay dt_1;
  Delaunay dt_5;

  // Create the test triangulations, and give short-hand names for
  // the vertices of interest.

  Vertex_handle v_0 = dt_0.infinite_vertex();
  Vertex_handle v_1 = dt_1.insert(Point( 0.0, 0.0 ));

  // All vertices in the triangulation with five points.
  std::vector<Vertex_handle> v_5(6);
  v_5[0] = dt_5.insert(Point(  0.0,  0.0  ));
  v_5[1] = dt_5.insert(Point( -0.5,  0.5  ));
  v_5[2] = dt_5.insert(Point(  1.0,  1.0  ));
  v_5[3] = dt_5.insert(Point(  1.0, -2.0  ));
  v_5[4] = dt_5.insert(Point( -1.0, -2.0  ));
  v_5[5] = dt_5.infinite_vertex();


  // Check helper functions compile as expected.
  make_degree(dt_0);

  make_dual_area(dt_0);
  make_dual_area(dt_0, CGAL::Finite_test_tag());
  make_dual_area(dt_0, CGAL::No_finite_test_tag());

  make_link_length(dt_0);
  make_link_length(dt_0, CGAL::Finite_test_tag());
  make_link_length(dt_0, CGAL::No_finite_test_tag());

  make_max_star_angle(dt_0);
  make_max_star_angle(dt_0, CGAL::Finite_test_tag());
  make_max_star_angle(dt_0, CGAL::No_finite_test_tag());

  make_min_star_angle(dt_0);
  make_min_star_angle(dt_0, CGAL::Finite_test_tag());
  make_min_star_angle(dt_0, CGAL::No_finite_test_tag());


  //-- Degree ----------------------------------------------------------------//

  Degree<Delaunay> degree;

  assert(  degree(v_0)    ==  0  );
  assert(  degree(v_1)    ==  0  );
  assert(  degree(v_5[0]) ==  4  );
  assert(  degree(v_5[1]) ==  4  );
  assert(  degree(v_5[5]) ==  4  );

  //-- Dual_area -------------------------------------------------------------//

  Dual_area<Delaunay, CGAL::No_finite_test_tag> 
    dual_area_0_a(dt_0),
    dual_area_1_a(dt_1),
    dual_area_5_a(dt_5);
  Dual_area<Delaunay, CGAL::Finite_test_tag>    
    dual_area_0_b(dt_0),
    dual_area_1_b(dt_1),
    dual_area_5_b(dt_5);
  Dual_area<Delaunay>                           
    dual_area_0_c(dt_0),
    dual_area_1_c(dt_1),
    dual_area_5_c(dt_5);

  // The Voronoi cell is unbounded in each of these cases.
  assert(  dual_area_0_b(v_0)    == infinity  );
  assert(  dual_area_1_b(v_1)    == infinity  );
  assert(  dual_area_0_c(v_0)    == infinity  );
  assert(  dual_area_5_c(v_5[1]) == infinity  );
  
  // Compare dual area with manually computed value.
  assert(  approx_eq(dual_area_5_a(v_5[0]), 2.64583)  );
  assert(  approx_eq(dual_area_5_b(v_5[0]), 2.64583)  );


  //-- Star_area -------------------------------------------------------------//

  Star_area<Delaunay, CGAL::No_finite_test_tag> 
    star_area_0_a(dt_0),
    star_area_1_a(dt_1),
    star_area_5_a(dt_5);
  Star_area<Delaunay, CGAL::Finite_test_tag>    
    star_area_0_b(dt_0),
    star_area_1_b(dt_1),
    star_area_5_b(dt_5);
  Star_area<Delaunay>                           
    star_area_0_c(dt_0),
    star_area_1_c(dt_1),
    star_area_5_c(dt_5);

  // Unbounded tests.
  assert(  star_area_0_b(v_0)    == infinity  );
  assert(  star_area_1_c(v_1)    == infinity  );
  assert(  star_area_5_b(v_5[2]) == infinity  );

  // Bounded tests.
  assert(  approx_eq(star_area_5_b(v_5[0]), 4.75)  );
  assert(  approx_eq(star_area_5_a(v_5[0]), 4.75)  );


  //-- Link_length -----------------------------------------------------------//

  Link_length<Delaunay, CGAL::No_finite_test_tag> 
    link_length_0_a(dt_0),
    link_length_1_a(dt_1),
    link_length_5_a(dt_5);
  Link_length<Delaunay, CGAL::Finite_test_tag>    
    link_length_0_b(dt_0),
    link_length_1_b(dt_1),
    link_length_5_b(dt_5);
  Link_length<Delaunay>                           
    link_length_0_c(dt_0),
    link_length_1_c(dt_1),
    link_length_5_c(dt_5);

  // The link length around the infinite vertex is zero for an empty
  // triangulation.
  assert(  approx_eq(link_length_0_b(v_0), 0)  );
  assert(  approx_eq(link_length_1_b(v_1), 0)  );

  // The link length is infinite if the link contains the infinite vertex.
  assert(  link_length_5_b(v_5[1]) == infinity  );

  // The link length of the infinite vertex is the same as the perimeter
  // of the convex hull in this case.
  assert(  approx_eq(link_length_5_b(v_5[0]), 9.130648586)  );
  assert(  approx_eq(link_length_5_a(v_5[5]), 9.130648586)  );


  //-- Max_star_angle --------------------------------------------------------//

  Max_star_angle<Delaunay, CGAL::No_finite_test_tag> 
    max_star_angle_0_a(dt_0),
    max_star_angle_1_a(dt_1),
    max_star_angle_5_a(dt_5);
  Max_star_angle<Delaunay, CGAL::Finite_test_tag>    
    max_star_angle_0_b(dt_0),
    max_star_angle_1_b(dt_1),
    max_star_angle_5_b(dt_5);
  Max_star_angle<Delaunay>                           
    max_star_angle_0_c(dt_0),
    max_star_angle_1_c(dt_1),
    max_star_angle_5_c(dt_5);

  
  // We have chosen to treat an infinite triangle as having an angle of zero
  // adjacent to the infinite point.
  assert(  approx_eq(max_star_angle_0_b(v_0), 3.14159)  );
  assert(  approx_eq(max_star_angle_1_c(v_1), 3.14159)  );

  // Compare with manually computed values.
  assert(  approx_eq(max_star_angle_5_a(v_5[0]), 1.89254)  );
  assert(  approx_eq(max_star_angle_5_b(v_5[0]), 1.89254)  );

  //-- Min_star_angle --------------------------------------------------------//

  Min_star_angle<Delaunay, CGAL::No_finite_test_tag> 
    min_star_angle_0_a(dt_0),
    min_star_angle_1_a(dt_1),
    min_star_angle_5_a(dt_5);
  Min_star_angle<Delaunay, CGAL::Finite_test_tag>    
    min_star_angle_0_b(dt_0),
    min_star_angle_1_b(dt_1),
    min_star_angle_5_b(dt_5);
  Min_star_angle<Delaunay>                           
    min_star_angle_0_c(dt_0),
    min_star_angle_1_c(dt_1),
    min_star_angle_5_c(dt_5);

  // We have chosen to treat an infinite triangle as having an angle of zero
  // adjacent to the infinite point.
  assert(  approx_eq(min_star_angle_0_b(v_0), 3.14159)  );
  assert(  approx_eq(min_star_angle_1_c(v_1), 3.14159)  );

  // Compare with manually computed values.
  assert(  approx_eq(min_star_angle_5_a(v_5[0]), 0.92729)  );
  assert(  approx_eq(min_star_angle_5_b(v_5[0]), 0.92729)  );
}
