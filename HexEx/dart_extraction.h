#include"typedefs.h"
#include <CGAL/intersections.h>
#include<fstream>
#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#endif
bool does_intersect(Tetrahedron_3 tet, Point_3 p){
 // return (tet.has_on_bounded_side(p));
  if(!(tet.has_on_bounded_side(p)||tet.has_on_boundary(p))) return false;
  return(tet.has_on_bounded_side(p)||(tet.has_on_boundary(p)&&(tet.has_on_bounded_side(p+Vector_3(1,0,0))||tet.has_on_bounded_side(p+Vector_3(0,1,0))||tet.has_on_bounded_side(p+Vector_3(0,0,1))||tet.has_on_bounded_side(p+Vector_3(1,1,0))||tet.has_on_bounded_side(p+Vector_3(1,0,1))||tet.has_on_bounded_side(p+Vector_3(0,1,1))||tet.has_on_bounded_side(p+Vector_3(1,1,1)))));
}

void extract_darts(LCC_3& input_tet_mesh, LCC_3& output_mesh, std::vector<Aff_transformation>& parametrization_matrix){
  for(LCC_3::Vertex_attribute_range::iterator it = ((output_mesh).vertex_attributes()).begin(), 
itend =((output_mesh).vertex_attributes()).end(); it != itend; it++){ //all the vertices known in output mesh so far
    for(LCC_3::One_dart_per_cell_range<3>::iterator it1 = (input_tet_mesh).one_dart_per_cell<3>().begin(), it1end = (input_tet_mesh).one_dart_per_cell<3>().end(); it1 !=  it1end; it1++){ //all the tets in the input mesh 
      std::vector<Point_3> points, params;
      /*for(LCC_3::One_dart_per_incident_cell_range<0,3>::iterator it2 = (input_tet_mesh).one_dart_per_incident_cell<0,3>(it1).begin(), it2end = (input_tet_mesh).one_dart_per_incident_cell<0,3>(it1).end(); it2 != it2end;it2++){
        Point p = (input_tet_mesh.info(it2)).parameters;
        points.push_back(p);
      }
      std::cout<<points[0]<<" "<<points[1]<<" "<<points[2]<< " " << points[3]<<std::endl; points.clear();*/
      Dart_handle temp = it1;
      Point p = (input_tet_mesh.info(temp)).parameters; params.push_back(p); p = input_tet_mesh.point(temp); points.push_back(p);
      temp = input_tet_mesh.alpha(it1,0);
      p = (input_tet_mesh.info(temp)).parameters; params.push_back(p); p = input_tet_mesh.point(temp); points.push_back(p);
      temp = input_tet_mesh.alpha(it1, 0, 1, 0);
      p = (input_tet_mesh.info(temp)).parameters; params.push_back(p); p = input_tet_mesh.point(temp); points.push_back(p);
      temp = input_tet_mesh.alpha(it1, 2, 0, 1, 0);
      p = (input_tet_mesh.info(temp)).parameters; params.push_back(p); p = input_tet_mesh.point(temp); points.push_back(p);
      Tetrahedron_3 tet(params[0], params[1], params[2], params[3]);
      Aff_transformation at_inv = (parametrization_matrix[(input_tet_mesh.info(it1)).cell_no]).inverse();
      //p = points[0];//((output_mesh).point_of_vertex_attribute(it)).transform(at);
     // tet.has_on_bounded_side(p);
      if(does_intersect(tet, p)){
        //at = at.inverse();
        //Point p0(), p1, p2, p3, p4, p5, p6, p7;
        Point p0 = points[0]; Point param = params[0];
        Point p1(param[0]+1, param[1], param[2]); p1 = p1.transform(at_inv);
        Point p2(param[0]+1, param[1]+1, param[2]); p2 = p2.transform(at_inv);
        Point p3(param[0], param[1]+1, param[2]); p3 = p3.transform(at_inv);
        Point p4(param[0], param[1]+1, param[2]+1); p4 = p4.transform(at_inv);
        Point p5(param[0], param[1], param[2]+1); p5 = p5.transform(at_inv);
        Point p6(param[0]+1, param[1], param[2]+1); p6 = p6.transform(at_inv);
        Point p7(param[0]+1, param[1]+1, param[2]+1); p7 = p7.transform(at_inv);
        output_mesh.make_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7);
        
      }
      
      /*std::ofstream of;
      of.open("output.off");
      CGAL::write_off(output_mesh, of); 
      of.close();

#ifdef CGAL_LCC_USE_VIEWER
  display_lcc(output_mesh);
#endif // CGAL_LCC_USE_VIEWER*/
    }
  }

std::ofstream of;
      of.open("output.off");
      CGAL::write_off(output_mesh, of); 
      of.close();
}
