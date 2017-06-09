//#ifndef HEXEXTR_H
//#define HEXEXTR_H
#include"HexExtr.h"
#include<iostream>
#include<string>
//namespace HexEx
//#endif

int main(){
  std::cout<<"Enter filename"<<std::endl;
  std::string str;
  std::cin>>str;
  HexEx::HexExtr h(str);
  std::cout<<h.directions[0]<<std::endl;
}
