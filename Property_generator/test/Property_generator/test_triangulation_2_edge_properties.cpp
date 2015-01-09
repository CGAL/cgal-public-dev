#include <iostream>
#include <limits>
#include <algorithm>
#include <utility>
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

  // Triangulations for testing. We consider two test cases.
  // A triangulation with one point, and a triangulation with five points.d
  Delaunay dt_2;
  dt_2.insert(Point(  0.0,  0.0  ));
  dt_2.insert(Point(  0.0,  1.0  ));

  Delaunay dt_5;
  Vertex_handle v_0 = dt_5.insert(Point(  0.0,  0.0  ));
  dt_5.insert(Point( -0.5,  0.5  ));
  dt_5.insert(Point(  1.0,  1.0  ));
  dt_5.insert(Point(  1.0, -2.0  ));
  dt_5.insert(Point( -1.0, -2.0  ));

  // Choose an infinite edge.  
  Edge e_1;

  for (Delaunay::All_edges_iterator e = dt_2.all_edges_begin();
       dt_2.is_infinite(*e); ++e)
    e_1 = *e;

  // Choose some edges to consider for the five point triangulation.
  Face_handle f_1 = dt_5.locate(Point(0, 0.1));
  Face_handle f_2 = dt_5.locate(Point(0,-0.1));

  Edge e_5_1 = std::make_pair( f_1,  f_1->index(v_0)           );
  Edge e_5_2 = std::make_pair( f_1, (f_2->index(v_0) + 1) % 3  );
  Edge e_5_3 = std::make_pair( f_2,  f_2->index(v_0)           );

  // Check that functions are valid.
  make_length(dt_2);
  make_length(dt_2, CGAL::Finite_test_tag());
  make_length(dt_2, CGAL::No_finite_test_tag());

  make_neighbor_area(dt_2);
  make_neighbor_area(dt_2, CGAL::Finite_test_tag());
  make_neighbor_area(dt_2, CGAL::No_finite_test_tag());

  make_dual_length(dt_2);
  make_dual_length(dt_2, CGAL::Finite_test_tag());
  make_dual_length(dt_2, CGAL::No_finite_test_tag());

  //-- Length ----------------------------------------------------------------//

  Length<Delaunay, CGAL::No_finite_test_tag> 
    length_a;
  Length<Delaunay, CGAL::Finite_test_tag>    
    length_2_b(dt_2),
    length_5_b(dt_5);
  Length<Delaunay>
    length_2_c(dt_2),
    length_5_c(dt_5);

  assert(  length_2_b(e_1) == infinity  );
  assert(  length_2_c(e_1) == infinity  );

  assert(  approx_eq(length_5_c(e_5_1), 1.58113)  );
  assert(  approx_eq(length_5_c(e_5_2), 0.70710)  );
  assert(  approx_eq(length_5_c(e_5_3), 2      )  );
  assert(  approx_eq(length_a(e_5_3),   2      )  );
  assert(  approx_eq(length_5_b(e_5_1), 1.58113)  );
  assert(  approx_eq(length_5_b(e_5_2), 0.70710)  );

  //-- Neighbor_area ---------------------------------------------------------//

  Neighbor_area<Delaunay, CGAL::No_finite_test_tag> 
    neighbor_area_5_a(dt_5);
  Neighbor_area<Delaunay, CGAL::Finite_test_tag>    
    neighbor_area_2_b(dt_2),
    neighbor_area_5_b(dt_5);
  Neighbor_area<Delaunay>
    neighbor_area_2_c(dt_2),
    neighbor_area_5_c(dt_5);

  assert(  neighbor_area_5_b(e_5_1) == infinity  );
  assert(  neighbor_area_5_c(e_5_1) == infinity  );
  assert(  neighbor_area_2_c(e_1)   == 0  );
  assert(  approx_eq(neighbor_area_5_b(e_5_2), 1.25)  );
  assert(  approx_eq(neighbor_area_5_a(e_5_2), 1.25)  );

  //-- Dual_length -----------------------------------------------------------//

  Dual_length<Delaunay, CGAL::No_finite_test_tag> 
    dual_length_5_a(dt_5);
  Dual_length<Delaunay, CGAL::Finite_test_tag>    
    dual_length_2_b(dt_2),
    dual_length_5_b(dt_5);
  Dual_length<Delaunay>
    dual_length_2_c(dt_2),
    dual_length_5_c(dt_5);

  assert(  approx_eq(dual_length_5_a(e_5_2), 2.00347)  );
  assert(  approx_eq(dual_length_5_b(e_5_2), 2.00347)  );
  assert(  approx_eq(dual_length_5_c(e_5_2), 2.00347)  );

  assert(  approx_eq(dual_length_2_b(e_1), 0)  );
  assert(  approx_eq(dual_length_2_c(e_1), 0)  );

  assert(  dual_length_5_b(e_5_1) == infinity  );
  assert(  dual_length_5_c(e_5_1) == infinity  );
  assert(  dual_length_5_b(e_5_3) == infinity  );

}
