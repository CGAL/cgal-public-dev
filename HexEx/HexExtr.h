#ifndef HEXEXTR_H
#define HEXEXTR_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Aff_transformation_3.h>
#include <iostream>
#include <algorithm>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K; 
typedef CGAL::Cartesian<double>::Vector_3                   Vector_3;
typedef CGAL::Direction_3                                   Direction;
typedef CGAL::Linear_cell_complex_for_generalized_map<3>    LCC_3;
typedef LCC_3::Dart_handle                                  Dart_handle;
typedef LCC_3::Point                                        Point;
typedef LCC_3::Vertex_attribute_handle                      Vertex_attribute_handle;
typedef CGAL::Aff_transformation_3                          Transformation;
typedef LCC_3:: size_type				    size_type		


namespace HexEx{

class Half_face_and_transition{
public:
  Half_face_and_transition(Dart_handle, Transformation);
  Dart_handle dart_handle;
  Transformation min_transformation;  
};

Half_face_and_transition::Half_face_and_transition(Dart_handle &dh, Transformation &tr){
 dart_handle = dh;
 min_transformation = tr;
}


template<class T> 
struct extract_transition_function: std::binary_function<T, vector<Half_face_and transition>>
{
  void operator()(typename T::Dart& d, &vector<Half_face_and_transition> all_faces_with_transitions){ 
    Dart_handle dh1 = lcc.dart_handle(d);
    Dart_handle dh2 = lcc.alpha(dh1,3);
    if(dh2 == NULL){//boundary

      return;
    }
    else{    
      if(!(lcc.is_marked(dh1, m))){	
        vector<Point> face1, face2;
        for (LCC_3::Dart_of_cell_range<2>::iterator it(lcc.darts_of_cell<2>(dh1).begin()), itend(lcc.darts_of_cell<2>(dh1).end()); it!=itend; ++it){
          lcc.mark(it, m);
          face1.push_back(lcc.point(it));
        }
        for (LCC_3::Dart_of_cell_range<2>::iterator it(lcc.darts_of_cell<2>(dh2).begin()), itend(lcc.darts_of_cell<2>(dh2).end()); it!=itend; ++it){
          lcc.mark(it, m);
          face2.push_back(lcc.point(it));
        }
        if(face1[0] == face2[0] && face1[1] == face2[1] && face1[2] == face2[2]){// transition function is identity.
          return;
        }
        Vector_3 c1 = face1[1].vector() - face1[0].vector();
        Vector_3 d1 = face1[2].vector() - face1[1].vector();
        Vector_3 c2 = face2[1].vector() - face2[0].vector();
        Vector_3 d2 = face2[2].vector() - face2[1].vector();
        auto min_dist = std::numeric_limits<double>::max();
        int min_trans_index = 0;  //the min transtion given by G[i]
        for(auto i = 0; i < 24; ++i){
                Vector_3 transf_c1 = G.transform(c1);
                auto dist1 = CGAL::scalar_product((c2 - transf_c1),(c2 - transf_c1));
                Vector_3 transf_d1 = G.transform(d1);
                auto dist2 = CGAL::scalar_product((d2 - transf_d1),(d2 - transf_d1));
                if(dist1 + dist2 < min_dist){
                    min_dist = dist1+dist2;
                    min_trans_index = i;
                }
         }
       Point new_point = G[min_trans_index].transform(face1[0]);
       Vector_3 t; t[0] = (int)((face2[0].vector())[0] - new_point[0]); t[1] = (int)((face2[0].vector())[1] - new_point[1]); t[2] = (int)((face2[0].vector())[2] - new_point[2]); //rounding to integer translation

       //Adding translation to the transformation matrix.
       Transformation final_transform_for_dh1(G[min_trans_index].m(0,0), G[min_trans_index].m(0,1), G[min_trans_index].m(0,2), t[0], G[min_trans_index].m(1,0), G[min_trans_index].m(1,1), G[min_trans_index].m(1,2), t[1], G[min_trans_index].m(2,0), G[min_trans_index].m(2,1), G[min_trans_index].m(2,2), t[2], 1);
       Transformation final_transform_for_dh2 = final_transform_for_dh1.inverse();
       Half_face_and_transition hfat1, hfat2;
       hfat1 = Half_face_and_transition(dh1, final_transition_dh1); hfat2 = Half_face_and_transition(dh2, final_transition_dh2);
       all_faces_with_transition.push_back(hfat1); all_faces_with_transition.push_back(hfat2);
//Need to store these
        
      }
    }
  }
};


template<class T> 
struct add_points_per_cell: public std::binary_function<T, vector<vector<Point>>> //is this required?
{
  void operator()(typename T::Dart& d, vector<vector<Point>> &points_in_each_cell){ 
    vector<Point> temp;
    for(typename T::template One_dart_per_incident_cell_range<0,3>:: const_iterator it=lcc.template one_dart_per_incident_cell<0,3>(lcc.dart_handle(d)).begin(), itend=lcc.template one_dart_per_incident_cell<0,3> (lcc.dart_handle(d)).end(); it!=itend; ++it){
      temp.push_back(lcc.point(it));
    }
    points_in_each_cell.push_back(temp);
  }
};

class HexExtr{
  public:
    HexExtr();
    HexExtr(std::string infilename){
     //load the file into LCC
      identity(1,0,0,0,1,0,0,0,1,1);
      direction[0] = Direction(1,0,0);
      direction[1] = Direction(0,1,0);
      direction[2] = Direction(0,0,1);
      direction[3] = Direction(-1,0,0);
      direction[4] = Direction(0,-1,0);
      direction[5] = Direction(0,0,-1);
      for(int i = 0;i < 6;i++) 
        for(int j = 0;j < 6;j++) 
          for(int k = 0;k < 6;k++) 
            if(CGAL::cross_product(directions[i].vector(),directions[j].vector()) == directions[k].vector()) 
              G.push_back(Transformation(directions[i].dx(), directions[i].dy(), directions[i].dz(), directions[j].dx(), directions[j].dy(), directions[j].dz(), directions[k].dx(), directions[k].dy(), directions[k].dz(), 1));

     // std::for_each(lcc.one_dart_per_cell<3>().begin(),lcc.one_dart_per_cell<3>().end(),add_points_per_cell(lcc, points_in_each_cell));
      size_type m = lcc.get_new_mark(); //to mark if transition function has been calculated.
    
      std::for_each(lcc.one_dart_per_cell<2>().begin(), lcc.one_dart_per_cell<2>().end(), extract_transition_function(lcc, all_faces_with_transitions));//each face is shared by max. two cells. 
      lcc.free_mark(m);
    }

    vector<Half_face_and_transition> all_faces_with_transitions;
    vector<vector<Point>> points_in_each_cell;
    vector<Direction> directions(6);
    Transformation identity;
    vector<Transformation> G; //chiral cubical symmetry group

}//namespace HexEx
#endif
