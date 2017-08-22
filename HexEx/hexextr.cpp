#ifndef FUNC_H
#define FUNC_H
#pragma once
#include"typedefs.h"
#include"hexextr.h"
#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#endif
//class HexExtr;

HexExtr::HexExtr(std::string infilename): identity(1,0,0,0,1,0,0,0,1,1){
 
//a parametrization function to test the extraction:
// dummy_parametrize(input_tet_mesh);

  load_mesh(infilename);
  
 // if(TIME) std::cout<<"Loading mesh: "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<" ms"<<std::endl;
//We would like to directly input parametrized meshes: so I am trying to first export LCC with parameters stored in dart_info using << in the lines 44 - 49, then try taking such a file as an input to our method in lines 51-57. If this works, parametrized_LCC file can be directly used to test.
  
}



void print_aff_transformation(Aff_transformation T){/**This function takes in an affine transformation as input and prints them: for debug purposes */
  for(int i=0; i<4; i++){
    for(int j = 0; j<4; j++)
      std::cout<<T.m(i,j)<<" ";
    std::cout<<std::endl; 
    }
  std::cout<<std::endl;
  return;
}

//namespace HexEx{



void dummy_parametrize(LCC_3& lcc){ /**dummy parametrization function: There is no guarantee that this could give good results, but for the sake of tests this might just work. 
* Every dart has an info as dart_info structure defined in typedefs.h. This contains information about the tet it belongs to (cell_no), and the parametric coordinates of the point. 
* This function updates the parametric coordinates. Any other relevant function can be used by changing the line marked with a comment below:*/
if(DEBUG) std::cout<<"Inside dummy_parameterize"<<std::endl;
  for(LCC_3::Dart_range::iterator it = (lcc.darts()).begin(), itend = (lcc.darts()).end(); it != itend; it++){
    //if(DEBUG) std::cout<<"Before creating info"<<std::endl;
    dart_info temp;
    //if(DEBUG) std::cout<<"After creating info"<<std::endl;
    temp.cell_no = 0;
  //  if(DEBUG) std::cout<<"Accessed info attribute"<<std::endl;
    Point_3 point(round(70*(lcc.point(it))[0]), round(70*(lcc.point(it))[1]), round(70*(lcc.point(it))[2]));// This line can be changed for trying different parametrizations.
//    if(DEBUG) std::cout<<"Created a point"<<std::endl;
    temp.parameters = point;
    //  if(DEBUG) std::cout<<"Accessed parameters"<<std::endl;
    temp.singular = false;
    //if(DEBUG) std::cout<<"Accessed singular"<<std::endl;
    temp.singular_edges = 0;
//    if(DEBUG) std::cout<<"Before accessing info"<<std::endl;
    lcc.info(it) = temp;
//std::cout<<lcc.point(it)<<"    "<<(lcc.info(it)).parameters<<std::endl;
  //  if(DEBUG) std::cout<<"After accessing info"<<std::endl;
  }
}

int HexExtr::set_dart_info(){
  for(LCC_3::Dart_range::iterator it = (input_tet_mesh.darts()).begin(), itend = (input_tet_mesh.darts()).end(); it != itend; it++){
    (input_tet_mesh.info(it)).cell_no = 0;
    (input_tet_mesh.info(it)).singular = false;
    (input_tet_mesh.info(it)).singular_edges = 0;
  }
//set_cell_nos();      

  int cell = 0;

//iterate through each tet in the input mesh:
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = input_tet_mesh.one_dart_per_cell<3>().begin(), itend = input_tet_mesh.one_dart_per_cell<3>().end(); it != itend; it++){ 

//iterate through all the darts of each input mesh, to assign a cell_no unique to darts of each tet. Two darts belonging to different tets will have different cell_no
    for(LCC_3::Dart_of_cell_range<3>::iterator it1 = input_tet_mesh.darts_of_cell<3>(it).begin(), it1end = input_tet_mesh.darts_of_cell<3>(it).end(); it1 != it1end; it1++){ 
      (input_tet_mesh.info(it1)).cell_no = cell;
    }
    cell++;
  }
  return cell;
}

void HexExtr::load_mesh(std::string infilename){
  std::ifstream in;

  input_tet_mesh.clear();
  in.open("tests/"+infilename);
  if (in.is_open())
  {  
    in>>input_tet_mesh;
    in.close();
  }
  else{
  std::cout<<"no"<<std::endl;
  }
}


void HexExtr::extract(){

//Preprocessing
  preprocess();  

//Sanitization
  sanitize(g);

//Topology Extraction: making hexahedrons in the output mesh incorporting vertex extraction and dart extraction in a single step:
  extract_hexes(parametrization_matrices);

//Connection extraction
  extract_connections();
    
}



int HexExtr::calculate_cell_type(LCC_3& lcc, Dart_handle dh){
  std::vector<Point> P;
  for(LCC_3::One_dart_per_incident_cell_range<0,3>::iterator it = lcc.one_dart_per_incident_cell<0,3>(dh).begin(), itend = lcc.one_dart_per_incident_cell<0,3>(dh).end(); it != itend; it++){
    P.push_back(lcc.point(it));
  }
  Vector_3 c1 = P[1] - P[0];
  Vector_3 c2 = P[2] - P[0];
  Vector_3 c3 = P[3] - P[0];
  double det = CGAL::determinant(c1, c2, c3);
  if(det == 0) return 0; //degenerate
  else if(det > 0) return 1; //proper
  else return -1; //flipped
}

//Comparison for finding points with  minimum and maximum of x, y and z
bool xfn(Point_3 i, Point_3 j){ return (i[0]<j[0]); }
bool yfn(Point_3 i, Point_3 j){ return (i[1]<j[1]); }
bool zfn(Point_3 i, Point_3 j){ return (i[2]<j[2]); }

//}
#endif

