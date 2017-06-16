//#ifndef HEXEXTR_H
//#define HEXEXTR_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include"HexExtr.h"
#include<iostream>
#include<string>
#include<CGAL/Aff_transformation_3.h>
#include"sanitization.h"
//namespace HexEx
//#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel K; 
typedef CGAL::Aff_transformation_3<K>                       Transformation;

void print(Transformation T, int p){
  std::cout<<std::endl<<"******** "<<p<<std::endl;
  for(int i=0; i<4; i++){
    for(int j = 0; j<4; j++)
      std::cout<<T.m(i,j)<<" ";
    std::cout<<std::endl; 
    }
  return;
}

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
  HexEx::HexExtr h(str);
  truncate_precision(h.input_tet_mesh);
  //for(int i = 0; i<24;i++) print(h.G[i], i);

}

