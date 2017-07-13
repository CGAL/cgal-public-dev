#include<CGAL/Linear_cell_complex_for_generalized_map.h>
#include<cstdlib>
#include<iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Vector_3.h>
#include"hexextr.h"
#include"typedefs.h"
#include"func.h"
#include"handles.h"


void propagate_parameters(LCC_3 &lcc, Vector_3 p, Vertex_handle v){
  for(LCC_3::One_dart_per_incident_cell_range<3,0>::iterator c = lcc.one_dart_per_incident_cell<3,0>(v.incident_dart).begin(), cend = lcc.one_dart_per_incident_cell<3,0>(v.incident_dart).end(); c!=cend; c++){// for the vertex in each cell
      //update the parameters
      Vector_3 to_prop; // =get_transition_function(c).transform();//parameter transformed to adjacent cell
      
  }
}

void truncate_precision(LCC_3 &lcc, Dart_handle& dh){
  
  for(LCC_3::One_dart_per_incident_cell_range<0,3>::iterator v=lcc.one_dart_per_incident_cell<0,3>(dh).begin(), vend=lcc.one_dart_per_incident_cell<0,3>(dh).end(); v!=vend; ++v){
    double max = std::numeric_limits<double>::min();
    Point p = lcc.point(v); // No, you should look this vertex up in grid map, find the parameters, and use them instead. 
    for(LCC_3::One_dart_per_incident_cell_range<3,0>::iterator cell = lcc.one_dart_per_incident_cell<3,0>(v).begin(), cellend = lcc.one_dart_per_incident_cell<3,0>(v).end(); cell != cellend; cell++){
      max = std::max(max, std::abs(p[0])); max = std::max(max, std::abs(p[1])); max = std::max(max, std::abs(p[2]));
    }
  //  Point p =  lcc.point(v);
    double delta = std::pow(2, std::ceil(std::log(max)/std::log(2)));
    std::cout<<"Delta: "<<delta<<std::endl;
    
    for(unsigned int i = 0; i < 3; ++i){
       int sign = std::signbit(p[i]) ? -1 : 1;
        volatile double tmp = p[i]+sign*delta;
        //Now store tmp into parameters // this "sanitizes" the parameters.
    }
    // Now we can propagate these parameters
    //propagate_parameters(lcc, );  
  }
 
   


}


bool is_edge_on_boundary_face(const LCC_3& lcc, Dart_const_handle dh){
 return lcc.is_free(dh, 3); 
}


double find_angle(const LCC_3& lcc, Dart_const_handle dh){
  Point p0 = lcc.point(lcc.alpha(dh,0,1,1));
  Point p1 = lcc.point(dh);
  Point p2 = lcc.point(lcc.alpha(dh,1));
  Point p3 = lcc.point(lcc.alpha(dh,1,1));
  Vector_3 v1 = p1 - p0;
  Vector_3 v2 = p2 - p0;
  Vector_3 v3 = CGAL::cross_product(v1, v2);
  v2 = CGAL::cross_product(v3, v1);
  Vector_3 v4 = v3-v1;
  Vector_3 v5 = CGAL::cross_product(v1, v4);
  v4 = CGAL::cross_product(v5, v1);
  return acos(std::max(std::min(CGAL::scalar_product(v2, v4), 1.0),-1.0));
}

Direction find_normal(LCC_3 &lcc, Dart_handle dh){
  Point p0 = lcc.point(dh);
  dh = lcc.alpha(dh,0,1);
  Point p1 = lcc.point(dh);
  dh = lcc.alpha(dh,0,1);
  Point p2 = lcc.point(dh);
  Vector_3 v1 = p1 - p0;
  Vector_3 v2 = p2 - p0;
  Vector_3 v0 = CGAL::cross_product(v1, v2);
  return v0.direction();
}



int calculate_edge_valence(const LCC_3& lcc, Dart_const_handle dh){
  double angle = 0;
  for (LCC_3::One_dart_per_incident_cell_range<2,1>::const_iterator it = ((lcc.one_dart_per_incident_cell<2,1>(dh)).begin()), itend = ((lcc.one_dart_per_incident_cell<2,1>(dh)).end()); it!=itend; ++it){
    if(!(is_edge_on_boundary_face(lcc, it))){
      Dart_const_handle dh1= it;
      double new_angle = find_angle(lcc, dh1);
      if(calculate_cell_type(lcc, dh) == -1)
        angle -= new_angle;
      else angle += new_angle;
    }
  }
 return (int)std::round(angle);
}



bool is_face_degenerate(const LCC_3 &lcc, LCC_3::Dart_const_handle dh){
  Point p0 = lcc.point(dh);
  LCC_3::Dart_const_handle dh1 = lcc.alpha(dh,0,1);
  Point p1 = lcc.point(dh1);
  LCC_3::Dart_const_handle dh2 = lcc.alpha(dh1,0,1);
  Point p2 = lcc.point(dh2);
  Point ori(0,0,0);
  return CGAL::determinant(p0 - ori,p1 - ori,p2 - ori) == 0;

}

Face_handle find_face_handle(HexExtr &h, Dart_handle& dh){
if((h.dart_in_face).find(dh) != (h.dart_in_face).end()) return  (h.dart_in_face).at(dh); //h.faces[0];
  else{
    for(std::vector<Face_handle>::iterator f = (h.faces).begin(), fend = (h.faces).end(); f!=fend; f++){
      for(std::vector<Dart_handle>::iterator it = (f->incident_darts).begin(), itend = (f->incident_darts).end(); it != itend; it++){
        if(*it == dh) return (h.dart_in_face).at(f->a_dart);
      }  
    }
   std::cout<<"Dart wasn't found!"<<std::endl; assert(false);
  }

}


bool is_edge_singular(HexExtr& h, Dart_handle dh){
  LCC_3 &lcc = h.input_tet_mesh;
  if(is_edge_on_boundary_face(lcc, dh)){
    if(calculate_cell_type(lcc, dh) == 0){ //degenerate cell 
      return true;
    }
  }
  Aff_transformation tran(1, 0, 0, 0, 1, 0, 0, 0, 1, 1);
  Aff_transformation identity(1, 0, 0, 0, 1, 0, 0, 0, 1, 1);
  for(LCC_3::One_dart_per_incident_cell_range<1,2>::iterator it = ((lcc.one_dart_per_incident_cell<1,2>(dh)).begin()), itend = ((lcc.one_dart_per_incident_cell<1,2>(dh)).end()); it!=itend; ++it){
      //Dart_handle transition_dart_handle_for_face;
    if(is_face_degenerate(lcc, it)){
      return true;
    }
    Face_handle fh = find_face_handle(h, it);//(h.input_tet_mesh, dh,-1);//
    tran = (h.faces_with_transitions[fh])*tran;  //transition(lcc, dh, tran, h.all_faces_with_transitions);// TODO- DONE?
  }
  return (!(tran.m(0,0)==1 && tran.m(0,1)==0 && tran.m(0,2)==0 && tran.m(0,3)==0 && tran.m(1,0)==0 && tran.m(1,1)==1 && tran.m(1,2)==0 && tran.m(1,3)==0 && tran.m(2,0)==0 && tran.m(2,1)==0 && tran.m(2,2)==1 && tran.m(2,3)==0 && tran.m(3,0)==0 && tran.m(3,1)==0 && tran.m(3,2)==0 && tran.m(3,3)==1));
}

void sanitize(){


}

