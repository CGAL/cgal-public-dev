#include"typedefs.h"
#include<fstream>
#include<CGAL/Linear_cell_complex_constructors.h>

#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#endif

int main(){
    LCC_3 lcc;
  // Create two hexahedra.
  //Dart_handle d1 = lcc.make_hexahedron(Point(-1, 0, 0), Point(0, 0, 0), Point(0, 1, 0), Point(-1, 1, 0), Point(-1, 1, 1), Point(-1, 0, 1), Point(0, 0, 1), Point(0, 1, 1));
 // lcc.reverse_orientation_connected_component(d1);
   Dart_handle d1 = lcc.make_hexahedron(Point(-1, 0, 0), Point(-1, 1, 0), Point(0, 1, 0), Point(0, 0, 0), Point(0, 0, 1), Point(-1, 0, 1), Point(-1, 1, 1), Point(0, 1, 1));
  Dart_handle d2 = lcc.make_hexahedron(Point(0, 0, 0), Point(0, 1, 0),  Point(1, 1, 0), Point(1, 0, 0), Point(1, 0, 1), Point(0, 0, 1), Point(0, 1, 1), Point(1, 1, 1));
  std::ofstream outfile;
  outfile.open("twohexes.off");
  CGAL::write_off(lcc, outfile);
  outfile.close();

  std::cout<<"Before 3-sew"<<std::endl;
  lcc.display_characteristics(std::cout);
#ifdef CGAL_LCC_USE_VIEWER
  display_lcc(lcc);
#endif
  
  Dart_handle dh1 = lcc.beta<2,1,1,2>(d1);
  Dart_handle dh2 = (d2);
  lcc.sew<3>(dh1, dh2);
  outfile.open("dart_sewn.off");
  CGAL::write_off(lcc, outfile);
  outfile.close();
/*  lcc.sew3_same_facets();
  outfile.open("same_facet_sew.off");
  CGAL::write_off(lcc, outfile);
  outfile.close();*/

  std::cout<<"After 3-sew"<<std::endl;
  lcc.display_characteristics(std::cout);
#ifdef CGAL_LCC_USE_VIEWER
  display_lcc(lcc);
#endif

  return 0;
}
