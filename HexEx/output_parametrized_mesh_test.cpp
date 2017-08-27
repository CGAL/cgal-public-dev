#include"hexextr.h" 
#include<iostream>
#include"typedefs.h"
#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#endif

void HexExtr::save_parametrized_mesh(){
for(LCC_3::Dart_range::iterator it = input_tet_mesh.darts().begin(), itend = input_tet_mesh.darts().end(); it != itend; it++){
input_tet_mesh.point(it) = (input_tet_mesh.info(it)).parameters;
}

}

int main(int argc, char** argv){
  std::string str;
  if (argc==1)
  {
    std::cout<<"Enter filename of input file from tests/: (eg: 2tets-1 (or) elephant-1 (or) fandisk-1)"<<std::endl;
    std::cin>>str;
  }
  else
  {
    str=argv[1];
  }

//To know how much time it takes to get the final output mesh, we start the timer:
  auto start_time = std::chrono::high_resolution_clock::now();

//creating HexExtr object by passing the name of the file- the parametrized mesh is loaded.
  HexExtr h(str);

//extract the mesh 
  std::cout<<str.substr(6)<<std::endl;
  
  h.extract(str.substr(6));


//post-processing is done, if required.
  h.refine();

//save the mesh
  h.set_parametrization(h.output_mesh);
  h.save_mesh(str.substr(6));

//is the output mesh valid?
  if((h.output_mesh).is_valid()) std::cout<<"Valid!"<<std::endl;
  else std::cout<<"Invalid!"<<std::endl;

//The timer after meshing done
  auto current_time = std::chrono::high_resolution_clock::now();

  std::cout << "Program has been running for " << std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count() << " ms" << std::endl;

// writing to a file
  h.save_parametrized_mesh();
  std::ofstream of;
  of.open("final_output.off");
  CGAL::write_off(h.output_mesh, of); 
  of.close();
  of.open("parametrized_input.off");
  CGAL::write_off(h.input_tet_mesh, of); 
  of.close();



  std::cout<<"*****FINAL OUTPUT MESH*****"<<std::endl;
  (h.output_mesh).display_characteristics(std::cout); std::cout<<std::endl;
  

  
  #ifdef CGAL_LCC_USE_VIEWER
    display_lcc(h.output_mesh);
  #endif

}

