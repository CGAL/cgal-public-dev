#ifndef FUNC_H
#define FUNC_H
#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include "Half_face_and_transition.h"
#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <cmath>
#include<cstdlib>
#include"typedefs.h"


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
Aff_transformation extract_transition_function(Dart_handle dh, const LCC_3& lcc, 
const std::vector<Aff_transformation>& G){ //the function returns the tranformation g applied to one tet to get to the vertices of the adjacent tet.
    Aff_transformation id(1,0,0,0,1,0,0,0,1,1);
    if(lcc.is_free(dh, 3)){//if the tet is on the boundary, it doesn't have a neighbouring tet. Returns identity.
      return id;
    }
    else{// Has an adjacent tet sharing a common face    

        Dart_const_handle dh1 = dh; 
        Dart_const_handle dh2 = lcc.alpha(dh1,3); //dart handle corresponding to the same edge but in the adjacent tet
    	std::vector<Point> face1, face2; 
        int i=0;
     
//adding the vertices of face1 (face of the first tet)
        face1.push_back(lcc.point(dh1));
        face1.push_back(lcc.point(lcc.alpha(dh1,1,0)));
        face1.push_back(lcc.point(lcc.alpha(dh1,1,0,1,0)));

//adding the vertices of face2 (face of the second tet)
        face2.push_back(lcc.point(dh2));
        face2.push_back(lcc.point(lcc.alpha(dh2,1,0)));
        face2.push_back(lcc.point(lcc.alpha(dh2,1,0,1,0)));
        //std::cout<<face1[0]<<" "<<face1[1]<<" "<<face1[2]<<std::endl;
        //std::cout<<face2[0]<<" "<<face2[1]<<" "<<face2[2]<<std::endl;
        if(face1[0] == face2[0] && face1[1] == face2[1] && face1[2] == face2[2]){// transition function is identity.
          return id;
        }
        Vector_3 c1 = face1[1] - face1[0];
        Vector_3 d1 = face1[2] - face1[1];
        Vector_3 c2 = face2[1] - face2[0];
        Vector_3 d2 = face2[2] - face2[1];
        auto min_dist = std::numeric_limits<double>::max();
        int min_trans_index = 0;  //the min transtion given by G[i]
        for(auto i = 0; i < 24; ++i){ //to find the transition function from cubical chiral symmetry group G
                Vector_3 transf_c1 = G[i].transform(c1);
                auto dist1 = CGAL::scalar_product((c2 - transf_c1),(c2 - transf_c1)); //Frobenius norm
                Vector_3 transf_d1 = G[i].transform(d1);
                auto dist2 = CGAL::scalar_product((d2 - transf_d1),(d2 - transf_d1)); //Frobenius norm
                if(dist1 + dist2 < min_dist){ 
                    min_dist = dist1+dist2;
                    min_trans_index = i;
                }
         }
         Point new_point = G[min_trans_index].transform(face1[0]);
        // Vector_3 t(std::round((face2[0])[0] - new_point[0]), std::round((face2[0])[1] - new_point[1]), std::round((face2[0])[2] - new_point[2])); //rounding to integer translation
         Vector_3 t((face2[0])[0] - new_point[0], (face2[0])[1] - new_point[1], (face2[0])[2] - new_point[2]);
       //Adding translation to the transformation matrix.
         Aff_transformation final_transform_for_dh1(G[min_trans_index].m(0,0), G[min_trans_index].m(0,1), G[min_trans_index].m(0,2), t[0], G[min_trans_index].m(1,0), G[min_trans_index].m(1,1), G[min_trans_index].m(1,2), t[1], G[min_trans_index].m(2,0), G[min_trans_index].m(2,1), G[min_trans_index].m(2,2), t[2], 1);
     
       return final_transform_for_dh1;

    }
  }

Aff_transformation get_parametrization_matrix(Point p, Point q, Point r, Point s, Point u, Point v, Point w, Point x){/**
* This fucntion returns the parametrization matrix. When this matrix is applied to four 3D points forming a tet, we get the corresponding coordinates of the tet in parametric space.
* p, q, r, s are the coordinates of the original points of the tet
* u, v, w, x  are the coordinates of the points of the tet in parametric space
*/
  Aff_transformation at1(q[0]-p[0], r[0]-p[0],s[0]-p[0], p[0], q[1]-p[1], r[1]-p[1],s[1]-p[1], p[1], q[2]-p[2], r[2]-p[2],s[2]-p[2], p[2], 1);
  Aff_transformation at2(v[0]-u[0], w[0]-u[0], x[0]-u[0], u[0], v[1]-u[1], w[1]-u[1], x[1]-u[1], u[1], v[2]-u[2], w[2]-u[2], x[2]-u[2], u[2], 1);
  return (at2*(at1.inverse()));
}

 //This will probably not be required:
int calculate_cell_type(const LCC_3& lcc, Dart_const_handle dh){
  std::vector<Point> P;
  for(int i=0;i<3;i++){
    P.push_back(lcc.point(dh));
    dh = lcc.alpha(dh, 0, 1);
  }
  Vector_3 c1 = P[1] - P[0];
  Vector_3 c2 = P[2] - P[0];
  Vector_3 c3 = P[3] - P[0];
  double det = CGAL::determinant(c1, c2, c3);
  if(det == 0) return 0; //degenerate
  else if(det > 0) return 1; //this is alright
  else return -1; //flipped
}


void dummy_parametrize(LCC_3& lcc){ /**dummy parametrization function: There is no guarantee that this could give good results, but for the sake of tests this might just work. 
* Every dart has an info as dart_info structure defined in typedefs.h. This contains information about the tet it belongs to (cell_no), and the parametric coordinates of the point. 
* This function updates the parametric coordinates. Any other relevant function can be used by changing the line marked with a comment below:*/
  for(LCC_3::Dart_range::iterator it = (lcc.darts()).begin(), itend = (lcc.darts()).end(); it != itend; it++){
    dart_info temp;
    temp.cell_no = 0;
    Point_3 point(round(10*(lcc.point(it))[0]), round(10*(lcc.point(it))[1]), round(10*(lcc.point(it))[2]));// This line can be changed for trying different parametrizations.
    temp.parameters = point;
    temp.singular = false;
    temp.singular_edges = 0;
    lcc.info(it) = temp;
  }
}

//}
#endif

