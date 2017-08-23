#include"hexextr.h" 
#include<iostream>
#include"typedefs.h"
#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#endif

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
  auto current_time = std::chrono::high_resolution_clock::now();
   std::cout << "Loading mesh: " << std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count() << " ms" << std::endl;

//preprocessing
  auto pre_time = std::chrono::high_resolution_clock::now();
  h.preprocess();
  current_time = std::chrono::high_resolution_clock::now();
  std::cout << "Preprocessing: " << std::chrono::duration_cast<std::chrono::milliseconds>(current_time - pre_time).count() << " ms" << std::endl;

//sanitize the parametrization
  pre_time = std::chrono::high_resolution_clock::now(); 
  h.sanitize();
  current_time = std::chrono::high_resolution_clock::now();
  std::cout << "Sanitizing parametrization: " << std::chrono::duration_cast<std::chrono::milliseconds>(current_time - pre_time).count() << " ms" << std::endl;

//Extract hexes 
  pre_time = std::chrono::high_resolution_clock::now();
  h.extract_hexes();
  current_time = std::chrono::high_resolution_clock::now();
  std::cout << "Extracting hexes: " << std::chrono::duration_cast<std::chrono::milliseconds>(current_time - pre_time).count() << " ms" << std::endl;

//extract the connections
  pre_time = std::chrono::high_resolution_clock::now();
  h.extract_connections();
  current_time = std::chrono::high_resolution_clock::now();
  std::cout << "Connecting hexes: " << std::chrono::duration_cast<std::chrono::milliseconds>(current_time - pre_time).count() << " ms" << std::endl;

//post-processing, if required.
  pre_time = std::chrono::high_resolution_clock::now();
  h.refine();
  current_time = std::chrono::high_resolution_clock::now();
  std::cout << "Postprocessing: " << std::chrono::duration_cast<std::chrono::milliseconds>(current_time - pre_time).count() << " ms" << std::endl;

//is the output mesh valid?
  if((h.output_mesh).is_valid()) std::cout<<"Valid!"<<std::endl;
  else std::cout<<"Invalid!"<<std::endl;

//The timer after meshing done
  current_time = std::chrono::high_resolution_clock::now();

  std::cout << "Program has been running for " << std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count() << " ms" << std::endl;

// writing to a file
  std::ofstream of;
  of.open("final_output.off");
  CGAL::write_off(h.output_mesh, of); 
  of.close();

  std::cout<<"*****FINAL OUTPUT MESH*****"<<std::endl;
  (h.output_mesh).display_characteristics(std::cout); std::cout<<std::endl;

  
  #ifdef CGAL_LCC_USE_VIEWER
    display_lcc(h.output_mesh);
  #endif

}

