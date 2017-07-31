#include"typedefs.h"
#include<unordered_map>

void connect_darts(LCC_3& lcc, Dart_handle dh1, Dart_handle dh2){
  for(LCC_3::Dart_of_cell_range<2>::iterator it1 = lcc.darts_of_cell<2>(dh1).begin(), it2 = lcc.darts_of_cell<2>(dh2).begin(), it1end = lcc.darts_of_cell<2>(dh1).end(), it2end = lcc.darts_of_cell<2>(dh2).end(); it1!= it1end && it2!=it2end; it1++, it2++){
    if(DEBUG) std::cout<<lcc.point(it1)<<" "<<lcc.point(it2)<<std::endl; 
    //lcc.link_beta<3>(it1, it2);
    lcc.beta<3>(it1) = it2; lcc.beta<3>(it2) = it1;
  }
}

void extract_connections(LCC_3& output_mesh, std::unordered_map<Point_3, Dart_handle>& hex_handles, std::unordered_map<Point_3, Point_3> output_points){/*
* output_points is a map from points in parametrized domain to points in hexahedral mesh domain
* hex_handles is a map from points in hex mesh domain to one dart handle to the hexahedron to which the point is incident.
* This was done to make sure that the different parametrizations in different tetrahedra didn't give rise to different points in hex domain for the same integer point.
* We connect each hexahedron at [x,y,z] in parametrized domain (say) to the next hexahedron (if present) at [x+1, y,z], [x,y+1,z] and [x,y,z+1].
*/
  for(std::unordered_map<Point_3, Point_3>::iterator it = output_points.begin(), itend = output_points.end(); it != itend; it++){
    //Dart_handle dh =  (*it).second;
    Point_3 p = (*it).second; // The point mapped to hex-domain
    if(hex_handles.find(p) == hex_handles.end()) continue; //This point does not hold a pointer to a hexahedron
    Dart_handle dh =  hex_handles[p]; // dh is the dart incident to the hex
    
    //along x:
    Point temp((*it).first.x()+1, (*it).first.y(), (*it).first.z()); // Point temp in integer point (parametrized) domain, which is at [x+1, y, z]
    if(output_points.find(temp) == output_points.end() || hex_handles.find(output_points[temp]) == hex_handles.end())/* move to the next point, if temp is not a point in the hex mesh (as found through extract_hexes function), or if temp is a point in the hex mesh, but does not have a hexahedron incident to it (a boundary point).*/ continue;

    Dart_handle x_dh = output_mesh.beta<1,1,2>(dh);;//alpha<0, 1, 0, 1, 2, 0>(dh); // dart in hexahedron at [x,y,z] corresponding to the shared face/edge/vertex between the two hexahedra 
    Dart_handle x_to_be_sewn = hex_handles[output_points[temp]]; //dart to the next hexahedron at [x+1, y,z]
    x_to_be_sewn = output_mesh.beta<2>(x_to_be_sewn); //alpha<2>(x_to_b3_sewn);// dart in hexahedron at [x+1,y,z] corresponding to the shared face/edge/vertex between the two hexahedra 
    //connect_darts(output_mesh, x_dh, x_to_be_sewn);
    //output_mesh.sew<3>(x_dh, x_to_be_sewn);
    output_mesh.sew<3>(x_to_be_sewn, x_dh);
    //along y:

    Point temp2((*it).first.x(), (*it).first.y()+1, (*it).first.z()); // Point temp in integer point (parametrized) domain, which is at [x, y+1, z]
    if(output_points.find(temp2) == output_points.end() || hex_handles.find(output_points[temp2]) == hex_handles.end())// move to the next point, if temp is not a point in the hex mesh (as found through extract_hexes function), or if temp is a point in the hex mesh, but does not have a hexahedron incident to it (a boundary point).
continue;
    Dart_handle y_dh = output_mesh.beta<2,1,1,2>(dh);//alpha<2, 0, 1, 0, 1, 0, 2>(dh); // dart in hexahedron at [x,y,z] corresponding to the shared face/edge/vertex between the two hexahedra 
    Dart_handle y_to_be_sewn = hex_handles[output_points[temp2]]; //dart to the next hexahedron at [x, y+1,z] - already corresponds to the dart having shared face/edge/vertex with y_dh
    //connect_darts(output_mesh, y_dh, y_to_be_sewn);
    //output_mesh.sew<3>(y_dh, y_to_be_sewn);
    output_mesh.sew<3>(y_to_be_sewn, y_dh);
    //along z:
   
    Point temp3((*it).first.x(), (*it).first.y(), (*it).first.z()+1); // Point temp in integer point (parametrized) domain, which is at [x, y, z+1]
    if(output_points.find(temp3) == output_points.end() || hex_handles.find(output_points[temp3]) == hex_handles.end())// move to the next point, if temp is not a point in the hex mesh (as found through extract_hexes function), or if temp is a point in the hex mesh, but does not have a hexahedron incident to it (a boundary point). 
      continue;
    Dart_handle z_dh = output_mesh.beta<1,2>(dh);//alpha<0, 1, 2>(dh); // dart in hexahedron at [x,y,z] corresponding to the shared face/edge/vertex between the two hexahedra 
    Dart_handle z_to_be_sewn = hex_handles[output_points[temp3]]; //dart to the next hexahedron at [x, y,z+1]
    z_to_be_sewn = output_mesh.beta<1,1,1,2>(z_to_be_sewn);//beta<1, 2>(z_to_be_sewn); // dart in hexahedron at [x,y,z+1] corresponding to the shared face/edge/vertex between the two hexahedra 
    //connect_darts(output_mesh, z_dh, z_to_be_sewn);
   // output_mesh.sew<3>(z_dh, z_to_be_sewn);
 output_mesh.sew<3>(z_to_be_sewn, z_dh);


   /*TODO- Previous method (not needed?)
    Point_3 px = output_mesh.point(x_dh);
    if(hex_handles.find(px) != hex_handles.end()){
      Dart_handle x_to_be_sewn = hex_handles[px];
      x_to_be_sewn = output_mesh.beta(x_to_be_sewn, 2);
      output_mesh.sew<3>(x_dh, x_to_be_sewn);
    }*/

 /*   //along y:
    Dart_handle y_dh = output_mesh.beta(dh, 2, 0, 1, 0, 1, 0, 2);
    Point_3 py = output_mesh.point(y_dh);
    if(hex_handles.find(py) != hex_handles.end()){
      Dart_handle y_to_be_sewn = hex_handles[py];
      output_mesh.sew<3>(y_dh, y_to_be_sewn);
    }
    //along z:
    Dart_handle z_dh = output_mesh.beta(dh, 0, 1, 2);
    Point_3 pz = output_mesh.point(z_dh);
    if(hex_handles.find(pz) != hex_handles.end()){
      Dart_handle z_to_be_sewn = hex_handles[pz];
      z_to_be_sewn = output_mesh.beta(z_to_be_sewn, 1, 2);
      output_mesh.sew<3>(z_dh, z_to_be_sewn);
    }
   */ 
  }
}
