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
typedef CGAL::Exact_predicates_inexact_constructions_kernel K; 
typedef CGAL::Aff_transformation_3<K>                       Transformation;
typedef CGAL::Linear_cell_complex_for_generalized_map<3>    LCC_3;
typedef CGAL::Direction_3<K>                                Direction;
typedef CGAL::Linear_cell_complex_for_generalized_map<3>    LCC_3;
typedef LCC_3::Dart_handle                                  Dart_handle;
typedef LCC_3::Point                                        Point;
typedef CGAL::Vector_3<K>                                   Vector_3;
typedef CGAL::Point_3<K>				                            Point_3;
typedef HexEx::Half_face_and_transition                     Half_face_and_transition;

namespace HexEx{
  Half_face_and_transition extract_transition_function(LCC_3::Dart &d, LCC_3 lcc, std::vector<Transformation> G){ //: lcc(hexex.lcc()), G(hexex.G()){return hfat1;}

  //void operator()(LCC_3::Dart& d){
    Dart_handle dh1 = lcc.dart_handle(d);
    Dart_handle dh2 = lcc.alpha(dh1,3);
    Transformation id(1,0,0,0,1,0,0,0,1,1);
    if(dh2 == NULL){//boundary
      Half_face_and_transition hfat(dh1, id);
      return hfat;
    }
    else{    
      //if(!(lcc.is_marked(dh1, m))){	
        std::vector<Point> face1, face2;
        for (LCC_3::Dart_of_cell_range<2>::iterator it((lcc.darts_of_cell<2>(dh1)).begin()), itend((lcc.darts_of_cell<2>(dh1)).end()); it!=itend; ++it){
          //lcc.mark(it, m);
          face1.push_back(lcc.point(it));
        }
        for (LCC_3::Dart_of_cell_range<2>::iterator it((lcc.darts_of_cell<2>(dh2)).begin()), itend((lcc.darts_of_cell<2>(dh2)).end()); it!=itend; ++it){
          //lcc.mark(it, m);
          face2.push_back(lcc.point(it));
        }
        if(face1[0] == face2[0] && face1[1] == face2[1] && face1[2] == face2[2]){// transition function is identity.
          return Half_face_and_transition(dh1, id);
        }
        Vector_3 c1 = face1[1] - face1[0];
        Vector_3 d1 = face1[2] - face1[1];
        Vector_3 c2 = face2[1] - face2[0];
        Vector_3 d2 = face2[2] - face2[1];
        auto min_dist = std::numeric_limits<double>::max();
        int min_trans_index = 0;  //the min transtion given by G[i]
        for(auto i = 0; i < 24; ++i){
                Vector_3 transf_c1 = G[i].transform(c1);
                auto dist1 = CGAL::scalar_product((c2 - transf_c1),(c2 - transf_c1));
                Vector_3 transf_d1 = G[i].transform(d1);
                auto dist2 = CGAL::scalar_product((d2 - transf_d1),(d2 - transf_d1));
                if(dist1 + dist2 < min_dist){
                    min_dist = dist1+dist2;
                    min_trans_index = i;
                }
         }
         Point new_point = G[min_trans_index].transform(face1[0]);
         Vector_3 t(std::round((face2[0])[0] - new_point[0]), std::round((face2[0])[1] - new_point[1]), std::round((face2[0])[2] - new_point[2])); //rounding to integer translation

       //Adding translation to the transformation matrix.
         Transformation final_transform_for_dh1(G[min_trans_index].m(0,0), G[min_trans_index].m(0,1), G[min_trans_index].m(0,2), t[0], G[min_trans_index].m(1,0), G[min_trans_index].m(1,1), G[min_trans_index].m(1,2), t[1], G[min_trans_index].m(2,0), G[min_trans_index].m(2,1), G[min_trans_index].m(2,2), t[2], 1);
      // Transformation final_transform_for_dh2 = final_transform_for_dh1.inverse();
//Need to store these
         Half_face_and_transition hfat1(dh1, final_transform_for_dh1);// hfat2 = Half_face_and_transition(dh2, final_transition_dh2);
//         all_faces_with_transition.push_back(hfat1); //all_faces_with_transition.push_back(hfat2);
       return hfat1;
        
      //} 
    }
  }
}
#endif

