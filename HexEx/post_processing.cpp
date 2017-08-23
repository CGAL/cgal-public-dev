#include"hexextr.h"

void HexExtr::refine(){ // extract() gives the correct hex mesh most of the time, post_processing() is mostly not even needed since Mesh_3 doesn't output mesh with flipped tets. So we check if post_processing is required, and execute if needed.
  if(post_processing_req()){
    post_processing();
  }
  return; //final mesh done
}

bool HexExtr::are_quads(){ //are all faces quads?
  for(LCC_3::One_dart_per_cell_range<2>::iterator it = output_mesh.one_dart_per_cell<2>().begin(), itend = output_mesh.one_dart_per_cell<2>().end(); it!=itend; it++){
    int count = output_mesh.one_dart_per_incident_cell<0,2>(it).size(); //counting the number of vertices incident on each face
    if(count == 4) continue;
    else return false; //some face has more than/less than 4 vertices
    }
  return true;
}

bool HexExtr::are_quad_strips(){ //Checks if only quad strips exist in the mesh - HexEx paper for definition
  for(LCC_3::One_dart_per_cell_range<2>::iterator it = output_mesh.one_dart_per_cell<2>().begin(), itend = output_mesh.one_dart_per_cell<2>().end(); it != itend; it++){
    int i = 0;
    for(LCC_3::One_dart_per_incident_cell_range<1,2>::iterator it1 = output_mesh.one_dart_per_incident_cell<1,2>(it).begin(), it1end = output_mesh.one_dart_per_incident_cell<1,2>(it).end(); it1 != it1end; it1++){
      if(!output_mesh.is_free<2>(it1)) i++; //counting if every edge incident on a face is shared bytwo adjacent faces- this would form a quad strip
    }
    if(i == 4) continue;
    else return false;
    }
  return true;
}


bool HexExtr::are_hexes(){ //Are all the 3-cells in the output_mesh hexahedral, i.e. have 6 incident faces?
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = output_mesh.one_dart_per_cell<3>().begin(), itend = output_mesh.one_dart_per_cell<3>().end(); it != itend; it++){
    int i = output_mesh.one_dart_per_incident_cell<2,3>(it).size();
    if(i == 6) continue;
    else return false;
    }
  return true;
}



int HexExtr::flipped_cells(){ //finds the number of flipped tet cells in input_tet_mesh to find if the parametrization was faulty- in which case post-processing is required.
  int count = 0;
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = input_tet_mesh.one_dart_per_cell<3>().begin(), itend =  input_tet_mesh.one_dart_per_cell<3>().end(); it != itend; it++){
    if(calculate_cell_type(input_tet_mesh, it) == -1) count++;
  }
  return count;
}


void HexExtr::annihilate_darts(){ // find dart-antdart pair in different hexes and erase them after adjusting the connections.

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

void HexExtr::merge_vertices(){/*
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

bool HexExtr::post_processing_req(){ //is post-processing required or is the output_mesh already correct?
  if(DEBUG) std::cout<<"Number of flipped cells: "<<flipped_cells()<<std::endl;
  return(!(are_quads()&& are_quad_strips() && are_hexes())||flipped_cells()>0);
}

void HexExtr::post_processing(){
  annihilate_darts();
  merge_vertices();
}

