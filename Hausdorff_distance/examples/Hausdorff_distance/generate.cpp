
#include <CGAL/point_generators_2.h>
#include <CGAL/Simple_cartesian.h>
#include <boost/lexical_cast.hpp>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;

int main(int argc, char* argv[])
{
  CGAL::Random_points_in_square_2<Point_2> rp2;
  int n = (argc>1)? boost::lexical_cast<int>(argv[1]):1000;
  
  std::cout.precision(17);
  for(int i=0; i < n; ++i){
    std::cout << *rp2++ << std::endl;
  }
  return 0;
}
