#include"typedefs.h"
#include <CGAL/intersections.h>

bool does_intersect(Tetrahedron_3 tet, Point_3 p){
  if(!tet.has_on_bounded_side(p)) return false;
  return(tet.has_on_bounded_side(p)||tet.has_on_bounded_side(p+Vector_3(1,0,0))||tet.has_on_bounded_side(p+Vector_3(0,1,0))||tet.has_on_bounded_side(p+Vector_3(0,0,1))||tet.has_on_bounded_side(p+Vector_3(1,1,0))||tet.has_on_bounded_side(p+Vector_3(1,0,1))||tet.has_on_bounded_side(p+Vector_3(0,1,1))||tet.has_on_bounded_side(p+Vector_3(1,1,1)));
}

void extract_darts(LCC_3& input_tet_mesh, LCC_3& output_mesh, std::vector<Aff_transformation>& parametrization_matrix){
  for(LCC_3::Vertex_attribute_range::iterator it = ((output_mesh).vertex_attributes()).begin(), 
itend =((output_mesh).vertex_attributes()).end(); it != itend; it++){ //all the vertices known in output mesh so far
    for(LCC_3::One_dart_per_cell_range<3>::iterator it1 = (input_tet_mesh).one_dart_per_cell<3>().begin(), it1end = (input_tet_mesh).one_dart_per_cell<3>().end(); it1 !=  it1end; it1++){ //all the tets in the input mesh
      std::vector<Point_3> points;
      for(LCC_3::One_dart_per_incident_cell_range<0,3>::iterator it2 = (input_tet_mesh).one_dart_per_incident_cell<0,3>(it1).begin(), it2end = (input_tet_mesh).one_dart_per_incident_cell<0,3>(it1).end(); it2 != it2end;it2++){
        Point p = (input_tet_mesh.info(it2)).parameters;
        points.push_back(p);
      }
      Tetrahedron_3 tet(points[0], points[1], points[2], points[3]);
      Aff_transformation at = parametrization_matrix[(input_tet_mesh.info(it1)).cell_no];
      Point_3 p = ((output_mesh).point_of_vertex_attribute(it)).transform(at);
     // if(does_intersect(tet, p)){
        at = at.inverse();
        output_mesh.make_hexahedron(p.transform(at), (p+Vector_3(1,0,0)).transform(at), (p+Vector_3(1,1,0)).transform(at), (p+Vector_3(0,1,0)).transform(at), (p+Vector_3(1,0,1)).transform(at), (p+Vector_3(1,1,1)).transform(at), (p+Vector_3(0,1,1)).transform(at), (p+Vector_3(0,0,1)).transform(at)); 
        
      //}
  
    }
  }
}
