#include<typedefs.h>

bool are_quads(LCC_3& lcc){ //are all faces quads?
  for(LCC_3::One_dart_per_cell_range<2>::iterator it = lcc.one_dart_per_cell<2>().begin(), itend = lcc.one_dart_per_cell<2>().end(); it!=itend; it++){
    /*int i = 0;
    for(LCC_3::Dart_of_cell_range<2>::iterator it1 = lcc.darts_of_cell<2>(it).begin(), it1end = lcc.darts_of_cell<2>(it).end(), it1 != it1end; it1++){
      i++;
    }*/
    int i = lcc.one_dart_per_incident_cell<0,2>(it).size();
//std::cout<<i<<std::endl;
    if(i == 4) continue;
    else return false;
    }
  return true;
}

bool are_quad_strips(LCC_3& lcc){
  for(LCC_3::One_dart_per_cell_range<2>::iterator it = lcc.one_dart_per_cell<2>().begin(), itend = lcc.one_dart_per_cell<2>().end(); it != itend; it++){
    int i = 0;
    for(LCC_3::One_dart_per_incident_cell_range<1,2>::iterator it1 = lcc.one_dart_per_incident_cell<1,2>(it).begin(), it1end = lcc.one_dart_per_incident_cell<1,2>(it).end(); it1 != it1end; it1++){
      if(!lcc.is_free<2>(it1)) i++;
    }
    if(i == 4) continue;
    else return false;
    }
  return true;
}


bool are_hexes(LCC_3& lcc){
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = lcc.one_dart_per_cell<3>().begin(), itend = lcc.one_dart_per_cell<3>().end(); it != itend; it++){
    int i = lcc.one_dart_per_incident_cell<2,3>(it).size();
    if(i == 6) continue;
    else return false;
    }
  return true;
}



int flipped_cells(LCC_3& lcc){
  int count = 0;
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = lcc.one_dart_per_cell<3>().begin(), itend =  lcc.one_dart_per_cell<3>().end(); it != itend; it++){
    if(calculate_cell_type(lcc, it) == -1) count++;
  }
  return count;
}


void annihilate_darts(LCC_3& output_mesh){
  for(LCC_3::Dart_range::iterator it = output_mesh.darts().begin(), itend = output_mesh.darts().end(); it != itend; it++){
    Dart_handle dh1 = it, dh2 = (output_mesh.beta<3>(it)); //all the darts of a flipped tet are marked as flipped, so we can find a dart-antidart pair only in different hexes.
    if((output_mesh.info(dh1)).flipped == true){
      if((output_mesh.info(dh2)).flipped == false){
        output_mesh.beta<2,2>(dh2) = output_mesh.beta<2>(dh1);
        output_mesh.beta<2,2>(dh1) = output_mesh.beta<2>(dh2);
        output_mesh.beta<1,1>(dh2) = output_mesh.beta<1>(output_mesh.opposite(dh1));
        output_mesh.beta<1,1>(dh1) = output_mesh.beta<1>(output_mesh.opposite(dh2));
        output_mesh.opposite(output_mesh.opposite(dh2)) =  output_mesh.opposite(dh1);
        output_mesh.opposite(output_mesh.opposite(dh1)) =  output_mesh.opposite(dh2);
        output_mesh.erase_dart(dh1); output_mesh.erase_dart(dh2);
      }
    }
  } 
}

void merge_vertices(LCC_3& output_mesh, std::unordered_map<Point_3, Point_3>& output_points){/*
we might not need this since we eliminate degenerate tets in the parametrization, and because we always check if a point already exists before adding it to the output_mesh. But assuming there might be duplicate vertices, we proceed with the steps 
*/
  for(LCC_3::Vertex_attribute_range::iterator v = (output_mesh.vertex_attributes()).begin(), vend = (output_mesh.vertex_attributes()).end(); v != vend; v++){
    for(LCC_3::Vertex_attribute_range::iterator v1 = v, v1end = (output_mesh.vertex_attributes()).end(); v1 != v1end; v1++){
      if(v!=v1){
        Point_3 p1 = output_mesh.point_of_vertex_attribute(v), p2 = output_mesh.point_of_vertex_attribute(v1); 
        if(p1 == p2){
          Point_3 p((p1.x()+p2.x())/2, (p1.y()+p2.y())/2, (p1.z()+p2.z())/2);
          output_mesh.erase_vertex_attribute(v); output_mesh.erase_vertex_attribute(v1);
          output_mesh.create_vertex_attribute(p);
        }
      }
    }
  }  
}

bool post_processing_req(LCC_3& lcc, LCC_3& input_mesh){
 // std::cout<<"Number of flipped cells: "<<flipped_cells(input_mesh)<<std::endl;
  return(!(are_quads(lcc)&& are_quad_strips(lcc) && are_hexes(lcc))||flipped_cells(input_mesh)>0);
}

void post_processing(LCC_3& lcc, std::unordered_map<Point_3, Point_3>& output_points){
  annihilate_darts(lcc);
  merge_vertices(lcc, output_points);
}
