//#ifndef HEXEXTR_H
//#define HEXEXTR_H
#include"HexExtr.h"
#include<iostream>
#include<string>
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
  HexEx::HexExtr h(str);
  std::cout<<h.directions[0]<<std::endl;
}
