#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include"hexextr.h" 
#include<iostream>
#include<string>
#include<chrono>
//#include<CGAL/Aff_transformation_3.h>
#include"sanitization.h"
#include"typedefs.h"
//namespace HexEx
//#endif


int main(int argc, char** argv){
  std::string str;
  if (argc==1)
  {
    std::cout<<"Enter filename: (eg: data/elephant.off)"<<std::endl;
    std::cin>>str;
  }
  else
  {
    str=argv[1];
  }
  auto start_time = std::chrono::high_resolution_clock::now();
  HexExtr h(str);
  auto current_time = std::chrono::high_resolution_clock::now();
  std::cout << "Program has been running for " << std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count() << " ms" << std::endl;
  //for(int i = 0; i<24;i++) print_aff_tranformation(h.G[i]);

}

