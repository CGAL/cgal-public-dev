//#ifndef HEXEXTR_H
//#define HEXEXTR_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include"hexextr.h" //checking?
#include<iostream>
#include<string>
//#include<CGAL/Aff_transformation_3.h>
#include"sanitization.h"
#include"typedefs.h"
//namespace HexEx
//#endif


int main(int argc, char** argv){
  std::string str;
  if (argc==1)
  {
    std::cout<<"Enter filename"<<std::endl;
    std::cin>>str;
  }
  else
  {
    str=argv[1];
  }
  HexExtr h(str);
  //for(int i = 0; i<24;i++) print_aff_tranformation(h.G[i]);

}

