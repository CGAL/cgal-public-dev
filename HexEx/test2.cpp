#include"typedefs.h"
#include<fstream>
#include<CGAL/Linear_cell_complex_constructors.h>
int main(){
    LCC_3 lcc;
  // Create two hexahedra.
  Dart_handle d1 = lcc.make_hexahedron(Point(-1, 0, 0), Point(0, 0, 0), Point(0, 1, 0), Point(-1, 1, 0), Point(-1, 1, 1), Point(-1, 0, 1), Point(0, 0, 1), Point(0, 1, 1));
  Dart_handle d2 = lcc.make_hexahedron(Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0), Point(0, 1, 1), Point(0, 0, 1), Point(1, 0, 1), Point(1, 1, 1));
  std::ofstream outfile;
  outfile.open("twohexes.off");
  CGAL::write_off(lcc, outfile);
  outfile.close();
  Dart_handle dh1 = lcc.beta<1,1,2>(d1);
  Dart_handle dh2 = lcc.beta<2>(d2);
  //lcc.sew<3>(dh1, dh2);
  lcc.sew3_same_facets();
  outfile.open("2hexessewn.off");
  CGAL::write_off(lcc, outfile);
  outfile.close();
  return 0;
}
