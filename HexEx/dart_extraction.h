#include"typedefs.h"
#include"hexextr.h"
#include <CGAL/intersections.h>

bool B2(Point_3 p, CGAL::Tetrahedron_3 tet){
  if(tet.has_on_bounded_side(p)||(tet.is_on_boundary(p)&&tet.has_on_bounded_side(p)))
}

bool does_intersect(CGAL::Tetrahedron_3 tet, Point_3 p){
  if(!tet.has_on_bounded_side(p)) return false;
  return(tet.has_on_bounded_side(p)||tet.has_on_bounded_side(p+Vector_3(1,0,0))||tet.has_on_bounded_side(p+Vector_3(0,1,0))||tet.has_on_bounded_side(p+Vector_3(0,0,1))||tet.has_on_bounded_side(p+Vector_3(1,1,0))||tet.has_on_bounded_side(p+Vector_3(1,0,1))||tet.has_on_bounded_side(p+Vector_3(0,1,1))||tet.has_on_bounded_side(p+Vector_3(1,1,1));
}

void extract_darts(HexExtr& h){
  vector<Direction> directions; directions[0] = Direction(1,0,0); directions[1] = Direction(0,1,0); directions[2] = Direction(0,0,1);
  Direction dir;
  for(LCC_3::Vertex_attribute_range::iterator it = ((h.output_mesh).vertex_attributes).begin(), 
itend =((h.output_mesh).vertex_attributes).end(); it != itend; it++){ //all the vertices
    for(LCC_3::One_dart_per_cell_range<3>::iterator it1 = (h.input_tet_mesh).one_dart_per_cell<3>().begin(), it1end = (h.input_tet_mesh).one_dart_per_cell<3>().end();; it1 !=  it1end; it1++){ //all the tets in the input mesh
      std::vector<Point_3> points;
      for(LCC_3::One_dart_per_incident_cell_range<0,3>::iterator it2 = (h.input_tet_mesh).one_dart_per_incident_cell<0,3>(it1).begin(), it2end = (h.input_tet_mesh).one_dart_per_incident_cell<0,3>(it1).end(); it2 != it2end;it2++){
        Point p = lcc.point(it2);
        points.push_back(p);
      }
      CGAL::Tetrahedron_3 tet(points[0], points[1], points[2], points[3]);
      Point_3 p = (h.output_mesh).point_of_vertex_attribute(it);
      if(does_intersect(tet, p)){
        Dart_handle dh = lcc.make_hexahedron(p, p+Vector_3(1,0,0), p+Vector_3(1,1,0), p+Vector_3(0,1,0), p+Vector_3(1,0,1), p+Vector_3(1,1,1), p+Vector_3(0,1,1), p+Vector_3(0,0,1));
        
      }
      /*if(B2(*it, tet) && B3(*it, tet) && B4(*it, tet)){
        //create_dart(); 
      }
*/
    }
  }
}
