#include"hexextr.h"

void HexExtr::extract_hexes(){/**
* This function combines the steps geometry extraction and topology extraction from the paper HexEx. 
* We iterate through each volume (in parametric space, stored in dart_info) in the input mesh, and find all the prospective vertices with integer coordinates in parametric space. If unit cubes formed by joining adjacent integer coordinates in the parametric space intersects the tetrahedron, we make a hexahedron (in our output mesh) corresponding to the inverse parametrization of the integer vertices of the cube. 
*/
  int num_hex = 0;
  output_mesh.clear();
//iterating through each tetrahedron
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = (input_tet_mesh.one_dart_per_cell<3>()).begin(), itend = (input_tet_mesh.one_dart_per_cell<3>()).end(); it != itend; it++){

    double minx = std::numeric_limits<double>::max(), miny = std::numeric_limits<double>::max(), minz = std::numeric_limits<double>::max(); 
    double maxx = std::numeric_limits<double>::min(), maxy = std::numeric_limits<double>::min(), maxz = std::numeric_limits<double>::min();
    std::vector<Point_3> params;
int i = 0;
for(LCC_3::One_dart_per_incident_cell_range<0,3>::iterator it1 = input_tet_mesh.one_dart_per_incident_cell<0,3>(it).begin(), it1end = input_tet_mesh.one_dart_per_incident_cell<0,3>(it).end(); it1 != it1end; it1++){
  Point param = (input_tet_mesh.info(it1)).parameters; params.push_back(param); i++;
}

    int orientation = calculate_cell_type(input_tet_mesh, it); 

//making a tetrahedron in the parametrized space
    Tetrahedron_3 tet(params[0], params[1], params[2], params[3]); 

//bounding box of the tetrahedron to find the minimum and maximum of x, y, z coordinates
    Bbox_3 box = tet.bbox(); 
    minx = round(box.xmin());
    maxx = round(box.xmax());
    miny = round(box.ymin());
    maxy = round(box.ymax());
    minz = round(box.zmin());
    maxz = round(box.zmax());

// inverse parametrization of the given tet
    Aff_transformation at_inv = (parametrization_matrices[(input_tet_mesh.info(it)).cell_no]).inverse(); 


/** iterating through the possible integer vertices. 
* If the integer vertex lies outside the tet, nothing happens.
* If the integer vertex lies on the boundary of the tet, and if the corresponding cube along positive x, y and z axis intersect with the tet, a hexahedron using inverse parametrization is created.
* If the integer vertex lies inside the tet, the corresponding cube is sure to intersect with the tet, so a hexahedron using inverse parametrization is created.
*/
    
    for(int i = minx; i<= maxx; i++){
      for(int j = miny; j<=maxy; j++){
        for(int k = minz; k<=maxz; k++){ 
          Point_3 p(i,j,k);
          if(does_intersect(tet, p)){ /*
* We create a map between integer grid points and inverse-parametrized points called output_points, to avoid repeated calculations and so that each grid point maps to a unique inverse parametrization (this need not happen due to numerical inefficiencies, leading to overlapping hexes) 
* The paper uses the transition function to find the parametrization matrix of adjacent tet- this could be done too, but we have the parametrization matrix itself readily available in parametrization_matrices- so we use that to know parametrization in other matrices (not necessarily adjacent)
*/          
            Point p0, param; 
            if(output_points.find(p) == output_points.end()){
              p0 = p.transform(at_inv); param = p;
              output_points.emplace(std::make_pair(p, p0));
            }
            else{
              p0 = output_points[p]; param = p;
            }

            Point p1(param[0]+1, param[1], param[2]); p = p1;
            if(output_points.find(p) == output_points.end()){
              if(tet.has_on_unbounded_side(p)){
                Aff_transformation *at_new = find_tet_parametrization(p);
                if(at_new == nullptr) p1 = p.transform(at_inv);
                else p1 = p.transform((*at_new).inverse());
              } 
              else p1 = p.transform(at_inv);
              output_points.emplace(std::make_pair(p, p1));
            }
            else{
              p1 = output_points[p]; 
            }

            Point p2(param[0]+1, param[1]+1, param[2]); p = p2;
            if(output_points.find(p) == output_points.end()){
              if(tet.has_on_unbounded_side(p)){
                Aff_transformation *at_new = find_tet_parametrization(p);
                if(at_new == nullptr) p2 = p.transform(at_inv);
                else p2 = p.transform((*at_new).inverse());
              } 
              else p2 = p.transform(at_inv);
              output_points.emplace(std::make_pair(p, p2));
            }
            else{
              p2 = output_points[p]; 
            }

            Point p3(param[0], param[1]+1, param[2]); p = p3;
            if(output_points.find(p) == output_points.end()){
              if(tet.has_on_unbounded_side(p)){
                Aff_transformation *at_new = find_tet_parametrization(p);
                if(at_new == nullptr) p3 = p.transform(at_inv);
               else p3 = p.transform((*at_new).inverse());
              } 
              else p3 = p.transform(at_inv);
              output_points.emplace(std::make_pair(p, p3));
            }
            else{
              p3 = output_points[p]; 
            }

            Point p4(param[0], param[1]+1, param[2]+1); p = p4;
            if(output_points.find(p) == output_points.end()){
              if(tet.has_on_unbounded_side(p)){
                Aff_transformation *at_new = find_tet_parametrization(p);
                if(at_new == nullptr) p4 = p.transform(at_inv);
                else p4 = p.transform((*at_new).inverse());
              } 
              else p4 = p.transform(at_inv);
              output_points.emplace(std::make_pair(p, p4));
            }
            else{
              p4 = output_points[p]; 
            }

            Point p5(param[0], param[1], param[2]+1); p = p5;
            if(output_points.find(p) == output_points.end()){
              if(tet.has_on_unbounded_side(p)){
                Aff_transformation *at_new = find_tet_parametrization(p);
                if(at_new == nullptr) p5 = p.transform(at_inv);
                else p5 = p.transform((*at_new).inverse());
              } 
              else p5 = p.transform(at_inv);
              output_points.emplace(std::make_pair(p, p5));
            }
            else{
              p5 = output_points[p]; 
            }

            Point p6(param[0]+1, param[1], param[2]+1); p = p6;
            if(output_points.find(p) == output_points.end()){
              if(tet.has_on_unbounded_side(p)){
                Aff_transformation *at_new = find_tet_parametrization(p);
                if(at_new == nullptr) p6 = p.transform(at_inv);
                else p6 = p.transform((*at_new).inverse());
              } 
              else p6 = p.transform(at_inv);
              output_points.emplace(std::make_pair(p, p6));
            }
            else{
              p6 = output_points[p]; 
            }

            Point p7(param[0]+1, param[1]+1, param[2]+1); p = p7;
            if(output_points.find(p) == output_points.end()){
              if(tet.has_on_unbounded_side(p)){
                Aff_transformation *at_new = find_tet_parametrization(p);
                if(at_new == nullptr) p7 = p.transform(at_inv);
                else p7 = p.transform((*at_new).inverse());
              } 
              else p7 = p.transform(at_inv);
              output_points.emplace(std::make_pair(p, p7));
            }
            else{
              p7 = output_points[p]; 
            }
 
            if(hex_handles.find(p0) == hex_handles.end()){ 
            Dart_handle dh = output_mesh.make_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7); num_hex++;
            for(LCC_3::Dart_of_cell_range<3>::iterator it5 = output_mesh.darts_of_cell<3>(dh).begin(),
it5end = output_mesh.darts_of_cell<3>(dh).end(); it5 != it5end; it5++){
              (output_mesh.info(it5)).flipped = (orientation == -1)? true : false;
            }
            hex_handles.emplace(std::make_pair(p0, dh));
}
          
          }
        }
      }
    }
  }


//display the characteristics of the output LCC.
  std::cout<<"Number of hexes created: "<<num_hex<<std::endl;
  std::cout<<"***Output mesh***"<<std::endl; output_mesh.display_characteristics(std::cout); std::cout<<std::endl;

/*
// output.off is used to visualize the output_mesh after this step. So LCC is written to output.off file.
  #ifdef CGAL_LCC_USE_VIEWER
    display_lcc(output_mesh);
  #endif // CGAL_LCC_USE_VIEWER
  std::ofstream of;
  of.open("output.off");
  CGAL::write_off(output_mesh, of); 
  of.close();
*/

} 



Aff_transformation* HexExtr::find_tet_parametrization(Point p){
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = input_tet_mesh.one_dart_per_cell<3>().begin(), itend = input_tet_mesh.one_dart_per_cell<3>().end(); it != itend; it++){
    std::vector<Point_3> point;
    for(LCC_3::One_dart_per_incident_cell_range<0,3>::iterator it1 = input_tet_mesh.one_dart_per_incident_cell<0,3>(it).begin(), it1end = input_tet_mesh.one_dart_per_incident_cell<0,3>(it).end(); it1 != it1end; it1++){
      point.push_back((input_tet_mesh.info(it1)).parameters);
    }

    Tetrahedron_3 tet(point[0], point[1], point[2], point[3]);
    if(tet.is_degenerate()) continue;
    if(!tet.has_on_unbounded_side(p)) return &(parametrization_matrices[(input_tet_mesh.info(it)).cell_no]);
  }
  return nullptr;
}

bool HexExtr::does_intersect(Tetrahedron_3 tet, Point_3 p){/**
* If p does not lie on the unbounded side of the tet, this function checks if a prospective unit volume cube in the positive x, y and z directions, with origin at p, would intersect the given tet. 
* If it lies on the unbounded side, we don't consider the volumes tobe intersecting: so the function returns false. 
*/
  if(tet.is_degenerate()) return false;
  else if(tet.has_on_bounded_side(p)) return true;
  else if(tet.has_on_boundary(p)&&(tet.has_on_bounded_side(Point_3(p[0]+1, p[1], p[2]))||tet.has_on_bounded_side(Point_3(p[0], p[1]+1, p[2]))||tet.has_on_bounded_side(Point_3(p[0], p[1], p[2]+1))||tet.has_on_bounded_side(Point_3(p[0]+1, p[1]+1, p[2]))||tet.has_on_bounded_side(Point_3(p[0]+1, p[1], p[2]+1))||tet.has_on_bounded_side(Point_3(p[0], p[1]+1, p[2]+1))||tet.has_on_bounded_side(Point_3(p[0]+1, p[1]+1, p[2]+1)))) return true;
  else return false;
}
