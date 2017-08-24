#ifndef FUNC_H
#define FUNC_H
#pragma once
#include"typedefs.h"
#include"hexextr.h"
#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#endif

HexExtr::HexExtr(std::string infilename): identity(1,0,0,0,1,0,0,0,1,1){

//infilename is the name of xml file containing information about parametrized mesh.
  load_mesh(infilename);
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

void dummy_parametrize(LCC_3& lcc){ /**dummy parametrization function: There is no guarantee that this could give good results, but for the sake of tests this might just work. 
* Every dart has an info as dart_info structure defined in typedefs.h. This contains information about the tet it belongs to (cell_no), and the parametric coordinates of the point. 
* This function updates the parametric coordinates. Any other relevant function can be used by changing the line marked with a comment below:*/

  for(LCC_3::Dart_range::iterator it = (lcc.darts()).begin(), itend = (lcc.darts()).end(); it != itend; it++){

    dart_info temp;

    temp.cell_no = 0;

    Point_3 point(round(70*(lcc.point(it))[0]), round(70*(lcc.point(it))[1]), round(70*(lcc.point(it))[2]));// This line can be changed for trying different parametrizations.

    temp.parameters = point;

    temp.singular = false;

    temp.singular_edges = 0;

    lcc.info(it) = temp;

  }
}

void HexExtr::set_parametrization(LCC_3& lcc){
  for(LCC_3::Dart_range::iterator it = (lcc.darts()).begin(), itend = (lcc.darts()).end(); it != itend; it++){
    Point p(0,0,0);
    (lcc.info(it)).parameters = p;
  } 
}

int HexExtr::set_dart_info(LCC_3& lcc){/**
* initializing the dart_info assiciated with each dart to certain values, and setting cell_no according to the tet to which the dart handle belongs to.
*/
  for(LCC_3::Dart_range::iterator it = (lcc.darts()).begin(), itend = (lcc.darts()).end(); it != itend; it++){
    (lcc.info(it)).cell_no = 0;
    (lcc.info(it)).singular = false;
    (lcc.info(it)).singular_edges = 0;
  }

  int cell = 0;

//iterate through each tet in the mesh:
  for(LCC_3::One_dart_per_cell_range<3>::iterator it = lcc.one_dart_per_cell<3>().begin(), itend = lcc.one_dart_per_cell<3>().end(); it != itend; it++){ 

//iterate through all the darts of each tet, to assign a cell_no unique to all darts of each tet. Two darts belonging to different tets will have different cell_no
    for(LCC_3::Dart_of_cell_range<3>::iterator it1 = lcc.darts_of_cell<3>(it).begin(), it1end = lcc.darts_of_cell<3>(it).end(); it1 != it1end; it1++){ 
      (lcc.info(it1)).cell_no = cell;
    }
    cell++;
  }
  return cell;
}

void HexExtr::load_mesh(std::string infilename){ /**
* loads input_tet_mesh with lcc stored in an xml file in tests/ containing parametrized mesh.
*/

  std::ifstream in;

  input_tet_mesh.clear();
  in.open(infilename);
  if (in.is_open())
  {  
    in>>input_tet_mesh;
    in.close();
  }
  else{
  std::cout<<"Can't load file."<<std::endl;
  }
}


void HexExtr::extract(std::string str){ //executes the four stages of hex-mesh extraction, other than post processing

//Preprocessing
  preprocess();  

//Sanitization
  sanitize();

//Topology Extraction: making hexahedrons in the output mesh incorporating vertex extraction and dart extraction in a single step:
  extract_hexes();

//Connection extraction
  extract_connections();

  
    
}

void HexExtr::save_mesh(std::string str){
  std::ofstream os;
  os.open("results/"+str);
  os<<(output_mesh);
  os.close();
}

int HexExtr::calculate_cell_type(LCC_3& lcc, Dart_handle dh){/**
*This function finds if a tet in the LCC on which dh dart handle is incident is degenerate, or flipped, or proper.
*/
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

#endif

