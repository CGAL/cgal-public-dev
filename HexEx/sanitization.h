#include<CGAL/Linear_cell_complex_for_generalized_map.h>
#include<cstdlib>
#include<iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Vector_3.h>
#include"typedefs.h"
#include"func.h"
void truncate_precision(LCC_3 &lcc){
  double max = std::numeric_limits<double>::min();
  for(LCC_3::Vertex_attribute_const_range::iterator v=lcc.vertex_attributes().begin(), vend=lcc.vertex_attributes().end(); v!=vend; ++v){
    Point p = lcc.point_of_vertex_attribute(v);
    max = std::max(max, std::abs(p[0])); max = std::max(max, std::abs(p[1])); max = std::max(max, std::abs(p[2]));
    //std::cout << lcc.point_of_vertex_attribute(v) << std::endl;
    //std::cout << max << std::endl;
    Dart_handle dh = lcc.create_dart(v);
     std::cout<<lcc.point(dh)<<std::endl;
    //(lcc.point(dh))[0] = 1000;
    std::cout<<lcc.point(dh)<<std::endl<<lcc.point(lcc.alpha(dh,1))<<std::endl<<std::endl;
  }
  double delta = std::pow(2, std::ceil(std::log(max)/std::log(2)));
  std::cout<<"Delta: "<<delta<<std::endl;
   


}

bool is_edge_on_boundary_face(const LCC_3& lcc, Dart_const_handle dh){
 return lcc.is_free(dh, 3); 
}


double findangle(const LCC_3& lcc, Dart_const_handle dh){
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
  for (LCC_3::One_dart_per_incident_cell_range<2,2>::const_iterator it = ((lcc.one_dart_per_incident_cell<2,2>(dh)).begin()), itend = ((lcc.one_dart_per_incident_cell<2,2>(dh)).end()); it!=itend; ++it){
    if(!(is_edge_on_boundary_face(lcc, it))){
      Dart_const_handle dh1= it;
      //Dart_handle dh2 = lcc.alpha(dh1, 2);
      double newangle = findangle(lcc, dh1);
      if(calculate_cell_type(lcc, dh) == -1)
        angle -= newangle;
      else angle += newangle;
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
/*
?? rotate_clockwise(LCC_3, Dart_handle current_cell, Dart_handle half_edge){


}
*/
/*
Half_face_and_transition find_HFAT(const LCC_3 &lcc, Dart_handle dh, std::vector<Half_face_and_transition> all_faces_with_transitions){
  for(std::vector<Half_face_and_transition>::iterator it = all_faces_with_transitions.begin(), itend = all_faces_with_transitions.end(); it != itend; it++){
    if(dh == (*it).dart_handle) return (*it);
    else if(dh == lcc.alpha((*it).dart_handle, 0)) return (*it);
    else if(dh == lcc.alpha((*it).dart_handle, 1)) return (*it);
    else if(dh == lcc.alpha((*it).dart_handle, 1, 0)) return (*it);
    else if(dh == lcc.alpha((*it).dart_handle, 0, 1)) return (*it);
    else if(dh == lcc.alpha((*it).dart_handle, 0, 1, 0)) return (*it);
}
return Half_face_and_transition();
}


Aff_transformation transition(const LCC_3 &lcc, Dart_handle dh, Aff_transformation tr, std::vector<Half_face_and_transition> all_faces_with_transitions){
  Aff_transformation newtr = find_HFAT(lcc, dh, all_faces_with_transitions).min_transformation; //to be defined
  return newtr*tr;
}



bool is_edge_singular(HexExtr& h, Dart_handle dh){
  const LCC_3 &lcc = h.input_tet_mesh;
  if(is_edge_on_boundary_face(lcc, dh)){
    if(calculate_cell_type(lcc, dh) == 0){ //degenerate cell 
      return true;
    }///TODO: check!!!
  Aff_transformation tran(1, 0, 0, 0, 1, 0, 0, 0, 1, 1);
  Aff_transformation identity(1, 0, 0, 0, 1, 0, 0, 0, 1, 1);
  for(LCC_3::One_dart_per_incident_cell_range<1,2>::const_iterator it = ((lcc.one_dart_per_incident_cell<1,2>(dh)).begin()), itend = ((lcc.one_dart_per_incident_cell<1,2>(dh)).end()); it!=itend; ++it){
      //Dart_handle transition_dart_handle_for_face;
      if(is_face_degenerate(lcc, it)){
        return true;
      }
      tran = transition(lcc, dh, tran, h.all_faces_with_transitions);// TODO- DONE?
   }
  return (!(tran.m(0,0)==1 && tran.m(0,1)==0 && tran.m(0,2)==0 && tran.m(0,3)==0 && tran.m(1,0)==0 && tran.m(1,1)==1 && tran.m(1,2)==0 && tran.m(1,3)==0 && tran.m(2,0)==0 && tran.m(2,1)==0 && tran.m(2,2)==1 && tran.m(2,3)==0 && tran.m(3,0)==0 && tran.m(3,1)==0 && tran.m(3,2)==0 && tran.m(3,3)==1));
}  
   
** start comment
    Dart_handle boundary_dart_handle_1 = dh;
    while(!is_edge_on_boundary_face(lcc, boundary_dart_handle_1)){
      Dart_handle current_dart_handle_for_cell = lcc.alpha(boundary_dart_handle_1, 3);
      Dart_handle transition_dart_handle_for_face;// =  rotate_clockwise(lcc, current_dart_handle_for_cell, dh); //todo
      boundary_face_dart_handle_1 = transition_dart_handle_for_face;
      }
    boundary_face_dart_handle_2 = lcc.alpha(boundary_face_dart_handle_1, 2); //same cell adjacent face
    
  //inside the loop
    Aff_transformation transformation_func(1,0,0,0,1,0,0,0,1,1); 
    for(int i = 0; i < ((lcc.one_dart_per_incident_cell<2>(dh)).size())-2; i++){
      Dart_handle trans_face_dart_handle = boundary_face_dart_handle_2;
      if(is_face_degenerate(trans_face_dart_handle)){
        return true;
      }
      transformation_func = transition(trans_face_dart_handle, transformation_func); //define this. get the transition function of the face and multiply by the given transformation_function (accumulate) to see if multiplying them makes the normals (found below) equal.
      boundary_face_dart_handle_2 = lcc.alpha(boundary_face_dart_handle_2, 3, 2);//Is this correct? Need to check again.
    }
    Direction normal1 = find_normal(lcc, boundary_face_dart_handle_1); 
    Direction normal2 = find_normal(lcc, boundary_face_dart_handle_2);
    normal1 = normal1.transform(transformation_func);
    double scal_prod = CGAL::scalar_product(normal1, normal2);
    if(scal_prod < 0.5 && calculate_edge_valence(lcc, dh) == 2) return false; //about the threshold?
    return (scal_prod < 0.5);
  }
  else{
    Aff_transformation transformation_func(1,0,0,0,1,0,0,0,1,1); 
    Af_transformation identity(1,0,0,0,1,0,0,0,1,1);
    for(LCC_3::Dart_of_cell_range<2>::const_iterator it((lcc.one_dart_per_incident_cell<2>(dh)).begin()), itend((lcc.one_dart_per_incident_cell<2>(dh)).end()); it!=itend; ++it){
      //Dart_handle transition_dart_handle_for_face;
      if(is_face_degenerate(*it)){
        return true;
      }
      transition(trans_face_dart_handle, transformation_func);     //doTransition(transitionFace, currentCell, tranFun);
    }

    if(transformation_func == identity){
      if(calculate_edge_valence(dh) > 6) return false;
      else return true;      
    }


   }



** end comment here

}
*/




