#ifndef GEO_EXTR_H
#define GEO_EXTR_H
#include"typedefs.h"
#include<vector>
#include"func.h"
#include"hexextr.h"
#include"handles.h"

void vertex_extraction(HexExtr &h){
  //through all vertex handles
  for(vector<Vertex_handle>::iterator it = (h.vertices).begin(), itend = (h.vertices).end(); it != itend; it++){
    Vector_3 param;  //using dart_handle find the parameters.
    if((param[0]-std::round(param[0]) == 0) && (param[1]-std::round(param[1]) == 0) && (param[2]-std::round(param[2]) == 0)){ // this is a h-vertex.
     
   }
  }

  for(vector<Edge_handle>::iterator it = (h.edges).begin(), itend = (h.edges).end(); it != itend; it++){
    Vector_3 u, v, temp;  //using dart_handle find the parameters of the two end points of the edge.
    double alpha = 0;
    while(alpha<=1){
      //temp1 = u*alpha; temp2 = v*(1-alpha);
      temp = (u*alpha) + (v*(1-alpha));
      if((temp[0]-std::round(temp[0]) == 0) && (temp[1]-std::round(temp[1]) == 0) && (temp[2]-std::round(temp[2]) == 0)){ // should it be this exact or do we allow for epsilon?
      
      }
      alpha += 0.0001;// is this okay?
    }
     
  }
 
  for(vector<Face_handle>::iterator it = (h.faces).begin(), itend = (h.faces).end(); it != itend; it++){
    Vector_3 u,v,w, temp;
    double alpha = 0, beta = 0;
    while(alpha<=1){
      while(beta<=1){
        temp = (u*alpha) + (v*beta) + (w*(1-alpha-beta));
        if((temp[0]-std::round(temp[0]) == 0) && (temp[1]-std::round(temp[1]) == 0) && (temp[2]-std::round(temp[2]) == 0)){ // should it be this exact or do we allow for epsilon?
      
        }
        beta += 0.001;
      }
      alpha += 0.001; 
    }

  }

  for(vector<Cell_handle>::iterator it = (h.cells).begin(), itend = (h.cells).end(); it != itend; it++){
    Vector_3 u,v,w,x,temp;
    double alpha = 0, beta = 0, gamma = 0;
    while(alpha<=1){
      while(beta<=1){
        while(gamma<=1){
          temp = (u*alpha) + (v*beta) + (w*gamma) + (x*(1-alpha-beta-gamma));
          if((temp[0]-std::round(temp[0]) == 0) && (temp[1]-std::round(temp[1]) == 0) && (temp[2]-std::round(temp[2]) == 0)){ // should it be this exact or do we allow for epsilon?
      
          }
          gamma += 0.001;
        }
        beta += 0.001;
      }
      alpha += 0.001; 
    }

  }

}

#endif
