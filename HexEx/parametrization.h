#include<iostream>
#include<string>
#include<fstream>
#include<cstring>
void parametrize_to_ply(std::string filename){
  std::ifstream infile;
  std::ofstream outfile;
  infile.open(filename);
  outfile.open("parameters.ply");
  std::sstring temp;
  int vertices;
  std::getline(infile, temp);
  while(temp.compare("end_header") != 0){
    if(temp.compare(0,14,"element vertex") == 0){
      vertices = atoi((temp.substr(15)).c_str());
    }
    outfile<<temp<<endl;
  }
  outfile<<end_header<<std::endl;
  int i = 0;
  while(i<vertices){
    i++;
    double a,b,c;
    infile>>a>>b>>c;
    outfile<<(int)round(a*1)<<" "<<(int)round(b*1)<<" "<<(int)round(c*1)<<" "<<endl;
  }

  while(std::getline(infile, temp)){
    outfile<<temp<<endl;
  }
}

