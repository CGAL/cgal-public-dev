#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <cmath>
#include "random.h"
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;

int main(int argc, char ** argv){
  if(argc != 3)
  {
    std::cout << "Usage: ./test <input filename> <percentage of points>" << std::endl;
	  return 0;
  }
    //creating tests with normals as input
  std::ifstream ifile(argv[1]);
  double n = std::stod(argv[2]);
  std::string str;
  std::vector<std::string> points;
  while(ifile){
    ifile >> str;
    points.push_back(str);
  }
  ifile.close();

  int m = std::floor(n/100.0 * (points.size()/6)); // number of points in output test

  std::ofstream ofile("new_test.xyz");
  for(int i = 0; i < m; i++)
  {
    int index = random_int(0, (points.size()/6) - 1);
    ofile << points[6 * index] << " " << points[6 * index + 1] << " " << points[6 * index + 2] << " " << points[6 * index + 3] << " " << points[6 * index + 4] << " " << points[6 * index + 5] << " " << std::endl;

  }
  ofile.close();
  return 0;
}
