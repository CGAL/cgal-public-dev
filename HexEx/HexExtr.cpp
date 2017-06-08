//#ifndef HEXEXTR_H
//#define HEXEXTR_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Aff_transformation_3.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <CGAL/Vector_3.h>
#include<CGAL/Point_3.h>
#include <triangulation_to_LCC.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K; 
typedef CGAL::Vector_3<K>                                   Vector_3;
typedef CGAL::Direction_3<K>                                Direction;
typedef CGAL::Linear_cell_complex_for_generalized_map<3>    LCC_3;
typedef LCC_3::Dart_handle                                  Dart_handle;
typedef LCC_3::Point                                        Point;
typedef LCC_3::Vertex_attribute_handle                      Vertex_attribute_handle;
typedef CGAL::Aff_transformation_3<K>                       Transformation;
typedef LCC_3:: size_type				    size_type;		
typedef CGAL::Point_3<K>				    Point_3;

namespace HexEx{

class Half_face_and_transition{
public:
  Half_face_and_transition(Dart_handle &dh, Transformation &tr){
    dart_handle = dh;
    min_transformation = tr;
  }
  Dart_handle dart_handle;
  Transformation min_transformation;  
};




Half_face_and_transition extract_transition_function(LCC_3::Dart &d, LCC_3 lcc, std::vector<Transformation> G){ //: lcc(hexex.lcc()), G(hexex.G()){return hfat1;}

  //void operator()(LCC_3::Dart& d){
    Dart_handle dh1 = lcc.dart_handle(d);
    Dart_handle dh2 = lcc.alpha(dh1,3);
    Transformation id(1,0,0,0,1,0,0,0,1,1);
    if(dh2 == NULL){//boundary

      return;
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
         Vector_3 t((int)((face2[0])[0] - new_point[0]), (int)((face2[0])[1] - new_point[1]), (int)((face2[0])[2] - new_point[2])); //rounding to integer translation

       //Adding translation to the transformation matrix.
         Transformation final_transform_for_dh1(G[min_trans_index].m(0,0), G[min_trans_index].m(0,1), G[min_trans_index].m(0,2), t[0], G[min_trans_index].m(1,0), G[min_trans_index].m(1,1), G[min_trans_index].m(1,2), t[1], G[min_trans_index].m(2,0), G[min_trans_index].m(2,1), G[min_trans_index].m(2,2), t[2], 1);
      // Transformation final_transform_for_dh2 = final_transform_for_dh1.inverse();
//Need to store these
         Half_face_and_transition hfat1(dh1, final_transform_for_dh1); //hfat2 = Half_face_and_transition(dh2, final_transition_dh2);
//         all_faces_with_transition.push_back(hfat1); //all_faces_with_transition.push_back(hfat2);
       return hfat1;
        
      //} 
    }
  }
/*
private:
  const LCC_3& lcc;
  const std::vector<Transformation> G;
  Half_face_and_transition hfat1;
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
*/


class HexExtr{
  public:
    HexExtr();
    HexExtr(std::string infilename): identity(1,0,0,0,1,0,0,0,1,1){
      lcc = load_off_to_LCC(infilename); //tetmesh to lcc
      directions.push_back(Direction(1,0,0));
      directions.push_back(Direction(0,1,0));
      directions.push_back(Direction(0,0,1));
      directions.push_back(Direction(-1,0,0));
      directions.push_back(Direction(0,-1,0));
      directions.push_back(Direction(0,0,-1));
      for(int i = 0;i < 6;i++) 
        for(int j = 0;j < 6;j++) 
          for(int k = 0;k < 6;k++) 
            if(CGAL::cross_product(directions[i].vector(),directions[j].vector()) == directions[k].vector()) 
              G.push_back(Transformation(directions[i].dx(), directions[i].dy(), directions[i].dz(), directions[j].dx(), directions[j].dy(), directions[j].dz(), directions[k].dx(), directions[k].dy(), directions[k].dz(), 1));

     // std::for_each(lcc.one_dart_per_cell<3>().begin(),lcc.one_dart_per_cell<3>().end(),add_points_per_cell(lcc, points_in_each_cell));
      //size_type m = lcc.get_new_mark(); //to mark if transition function has been calculated.
    
      //_faces_with_transitions.push_back(std::for_each(lcc.one_dart_per_cell<2>().begin(), lcc.one_dart_per_cell<2>().end(), extract_transition_function<HexExtr>(*this)));//each face is shared by max. two cells. 
      //lcc.free_mark(m);
    for(LCC_3::One_dart_per_cell_range<3>::iterator it = lcc.one_dart_per_cell<3>().begin(), itend = lcc.one_dart_per_cell<3>().end(); it != itend; it++){
		all_faces_with_transitions.push_back(extract_transition_function(it, lcc, G));
		
	} 
    
    }

    std::vector<Half_face_and_transition> all_faces_with_transitions;
    std::vector<std::vector<Point>> points_in_each_cell;
    std::vector<Direction> directions;
    Transformation identity;//(1,0,0,0,1,0,0,0,1,1);
    LCC_3 lcc;
    std::vector<Transformation> G; //chiral cubical symmetry group

}//namespace HexEx
//#endif

int main(){
  std::cout<<"Enter filename"<<std:endl;
  std::string str;
  std::cin>>str;
  HexExtr h(str);
  std::out<<h.direction[0]<<std::endl;
}
