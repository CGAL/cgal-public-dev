#include<CGAL/Linear_cell_complex_for_generalized_map.h>
#include<cstdlib>
#include<iostream>

typedef CGAL::Linear_cell_complex_for_generalized_map<3>    LCC_3
typedef LCC_3::Dart_handle                                  Dart_handle;
typedef LCC_3::Dart_const_handle                                  Dart_const_handle;

double findangle(const LCC_3& lcc, Dart_const_handle dh){
  Point p0 = lcc.point(lcc.alpha(lcc.alpha(lcc.alpha(dh,0),1),1));
  Point p1 = lcc.point(dh);
  Point p2 = lcc.point(lcc.alpha(dh,1));
  Point p3 = lcc.point(lcc.alpha(dh,1),1);
  Vector_3 v1 = p1 - p0;
  Vector_3 v2 = p2 - p0;
  Vector_3 v3 = CGAL::cross_product(v1, v2);
  v2 = CGAL::cross_product(v3, v1);
  Vector_3 v4 = v3-v1;
  Vector_3 v5 = CGAL::cross_product(v1, v4);
  v4 = CGAL::cross_product(v5, v1);
  return acos(std::max(std::min(CGAL::scalar_product(v2, v4), 1.0),-1.0));
}

int calculate_edge_valence(const LCC_3& lcc, Dart_const_handle dh){
  double angle = 0;
  for (LCC_3::Dart_of_cell_range<2>::const_iterator it((lcc.one_dart_per_incident_cell<2>(dh)).begin()), itend((lcc.one_dart_per_incident_cell<2>(dh)).end()); it!=itend; ++it){
    if(!(is_edge_boundary(lcc, *it))){
      Dart_const_handle dh1= *it;
      //Dart_handle dh2 = lcc.alpha(dh1, 2);
      double newangle = findangle(lcc, dh1);
      if(isCellFlipped(ch))
        angle -= newangle;
      else angle += newangle;
    }
  }
}

bool is_edge_boundary(const LCC_3& lcc, Dart_const_handle dh){
 return lcc.is_free(dh, 3); 
}





