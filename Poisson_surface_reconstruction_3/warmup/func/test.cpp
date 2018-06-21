#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include "random.h"
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;

int main(int argc, char ** argv){
  if(argc != 3)
  {
    std::cout << "Usage: ./test <output filename> <number of points>" << std::endl;

	return 0;
  }
    //creating tests with normals as input
  std::ofstream ofile(argv[1]);
  double n = std::stoi(argv[2]);
  ofile << n <<std::endl;
  for(int i = 0; i < n; i++)
  {
	const Vector vec = ::random_unit_vec<Vector>();
  const double rad = random_double(0, 0.8);
  //double value = 5.0 * rad - 5.0;
  const Point p1 = CGAL::ORIGIN + vec * rad;
  const Point p2 = CGAL::ORIGIN + vec * 1.5;
//  ofile << p1 << " " << value << std::endl;
//  ofile << p2 << " " << 1.0 << std::endl;
  double value1 = p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] - 1.0;
  double value2 = p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2] - 1.0;
  ofile << p1 << " " << value1 << std::endl;
  ofile << p2 << " " << value2 << std::endl;

/*  const double v1 = -1.0;
	const double v2 = 1.0;
	const Point p1 = CGAL::ORIGIN + vec * 0.8;
	const Point p2 = CGAL::ORIGIN + vec * 1.2;
	ofile << p1 << " " << v1 << std::endl;
	ofile << p2 << " " << v2 << std::endl;
*/
  }
  return 0;
}
