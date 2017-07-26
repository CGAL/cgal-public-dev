#include"typedefs.h"
#include<unordered_map>
void extract_connections(LCC_3& output_mesh, std::unordered_map<Point_3, Dart_handle>& hex_handles){
  for(std::unordered_map<Point_3, Dart_handle>::iterator it = hex_handles.begin(), itend = hex_handles.end(); it != itend; it++){
    Dart_handle dh =  (*it).second;
    //along x:
    Dart_handle x_dh = output_mesh.alpha(dh, 0, 1, 0, 1, 2, 0);
    Point_3 px = output_mesh.point(x_dh);
    if(hex_handles.find(px) != hex_handles.end()){
      Dart_handle x_to_be_sewn = hex_handles[px];
      x_to_be_sewn = output_mesh.alpha(x_to_be_sewn, 2);
      output_mesh.sew<3>(x_dh, x_to_be_sewn);
    }
    //along y:
    Dart_handle y_dh = output_mesh.alpha(dh, 2, 0, 1, 0, 1, 0, 2);
    Point_3 py = output_mesh.point(y_dh);
    if(hex_handles.find(py) != hex_handles.end()){
      Dart_handle y_to_be_sewn = hex_handles[py];
      output_mesh.sew<3>(y_dh, y_to_be_sewn);
    }
    //along z:
    Dart_handle z_dh = output_mesh.alpha(dh, 0, 1, 2);
    Point_3 pz = output_mesh.point(z_dh);
    if(hex_handles.find(pz) != hex_handles.end()){
      Dart_handle z_to_be_sewn = hex_handles[pz];
      z_to_be_sewn = output_mesh.alpha(z_to_be_sewn, 1, 2);
      output_mesh.sew<3>(z_dh, z_to_be_sewn);
    }
    
  }
}
