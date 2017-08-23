#include"hexextr.h"

void HexExtr::extract_connections(){/*
* output_points is a map from points in parametrized domain to points in hexahedral mesh domain
* hex_handles is a map from points in hex mesh domain to one dart handle to the hexahedron to which the point is incident.
* This was done to make sure that the different parametrizations in different tetrahedra didn't give rise to different points in hex domain for the same integer point.
* We connect each hexahedron at [x,y,z] in parametrized domain (say) to the next hexahedron (if present) at [x+1, y,z], [x,y+1,z] and [x,y,z+1].
*/

// output_mesh.sew3_same_facets(); return;
 
  for(std::unordered_map<Point_3, Point_3>::iterator it = output_points.begin(), itend = output_points.end(); it != itend; it++){

 // The point mapped to hex-domain
    Point_3 p = (*it).second;

    if(hex_handles.find(p) == hex_handles.end()) continue; //This point does not hold a pointer to a hexahedron

    Dart_handle dh =  hex_handles[p]; // dh is the dart incident to the hex
    

/****along x:****/
// Point temp in integer point (parametrized) domain, which is at [x+1, y, z]
    Point temp((*it).first.x()+1, (*it).first.y(), (*it).first.z()); 


//if temp is not a point in the hex mesh (as found through extract_hexes function), or if temp is a point in the hex mesh, but does not have a hexahedron incident to it (a boundary point), then  move to the next point.
    if(output_points.find(temp) == output_points.end() || hex_handles.find(output_points[temp]) == hex_handles.end()) continue;

// dart in hexahedron at [x,y,z] corresponding to the shared face/edge/vertex between the two hexahedra 
    Dart_handle x_dh = output_mesh.beta<1,1,2>(dh);

//dart to the next hexahedron at [x+1, y,z]
    Dart_handle x_to_be_sewn = hex_handles[output_points[temp]];

// dart in hexahedron at [x+1,y,z] corresponding to the shared face/edge/vertex between the two hexahedra
    x_to_be_sewn = output_mesh.beta<2>(x_to_be_sewn); 
   
//3-sew the two darts:
    output_mesh.sew<3>(x_to_be_sewn, x_dh);


/****along y:****/
// Point temp in integer point (parametrized) domain, which is at [x, y+1, z]
    Point temp2((*it).first.x(), (*it).first.y()+1, (*it).first.z());


//if temp is not a point in the hex mesh (as found through extract_hexes function), or if temp is a point in the hex mesh, but does not have a hexahedron incident to it (a boundary point), then  move to the next point.
    if(output_points.find(temp2) == output_points.end() || hex_handles.find(output_points[temp2]) == hex_handles.end()) continue;

// dart in hexahedron at [x,y,z] corresponding to the shared face/edge/vertex between the two hexahedra 
    Dart_handle y_dh = output_mesh.beta<2,1,1,2>(dh);

//dart to the next hexahedron at [x, y+1,z] - already corresponds to the dart having shared face/edge/vertex with y_dh
    Dart_handle y_to_be_sewn = hex_handles[output_points[temp2]];

//3-sew the two darts:
    output_mesh.sew<3>(y_to_be_sewn, y_dh);
    

/****along z:****/
// Point temp in integer point (parametrized) domain, which is at [x, y, z+1]
    Point temp3((*it).first.x(), (*it).first.y(), (*it).first.z()+1);

//if temp is not a point in the hex mesh (as found through extract_hexes function), or if temp is a point in the hex mesh, but does not have a hexahedron incident to it (a boundary point), then  move to the next point.
    if(output_points.find(temp3) == output_points.end() || hex_handles.find(output_points[temp3]) == hex_handles.end()) continue;

    Dart_handle z_dh = output_mesh.beta<1,2>(dh);

//dart to the next hexahedron at [x, y,z+1]
    Dart_handle z_to_be_sewn = hex_handles[output_points[temp3]];

    z_to_be_sewn = output_mesh.beta<1,1,1,2>(z_to_be_sewn);

//3-sew the two darts:
    output_mesh.sew<3>(z_to_be_sewn, z_dh);

  }
}
