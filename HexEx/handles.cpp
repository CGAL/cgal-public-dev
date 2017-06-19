#include"handles.h"
Face_handle::Face_handle(LCC_3& lcc, Dart_handle& dh, int i){
      //LCC_3 lcc = h->input_tet_mesh;
      a_dart = dh;
      //(h->dart_in_face).emplace(dh, this);
      for(LCC_3::Dart_of_cell_range<3>:: iterator it=lcc.darts_of_cell<3>(dh).begin(), 
itend=lcc.darts_of_cell<3>(dh).end(); it!=itend; ++it){
        incident_darts.push_back(it);
       // (h->dart_in_face).emplace(it, this);
      }
      enumeration = i;
}
