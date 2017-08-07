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


/*
int flipped_cells(LCC_3& lcc){
  int count = 0;
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = lcc.one_dart_per_cell<3>().begin(), itend =  lcc.one_dart_per_cell<3>().end(); it != itend; it++){
    if(calculate_cell_type(lcc, it) == -1) count++;
  }
  return count;
}
*/

void annihilate_darts(LCC_3& output_mesh){
  
}

bool post_processing_req(LCC_3& lcc, LCC_3& input_mesh){
 // std::cout<<"Number of flipped cells: "<<flipped_cells(input_mesh)<<std::endl;
  return(!(are_quads(lcc)&& are_quad_strips(lcc) && are_hexes(lcc)));
}
