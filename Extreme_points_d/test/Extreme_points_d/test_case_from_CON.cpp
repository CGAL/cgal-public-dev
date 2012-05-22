#include <vector>
#include <fstream>
#include <string>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <cassert>

#include <CGAL/config.h>
#include <CGAL/Timer.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Extreme_points_d.h>

#include <CGAL/Random.h>

typedef CGAL::Cartesian_d<double> Kernel_d;
typedef Kernel_d::Point_d Point_d;

template <class InputStream, class OutputIterator>
void read_input(InputStream &in, OutputIterator out, int &n, int &d) {
	// read *.con file
	std::string s;
	do {
		getline(in, s);
	} while (s[0] == '%'); // skipping comments
	std::stringstream ss(s);
	ss >> d >> n;
	
	std::vector<Point_d> points(n);
	for (int i = 0; i < n; ++i) {
		std::vector<double> p(d);
		for (int j = 0 ; j < d; ++j)
			in >> p[j];
		*out++ = Point_d(d, p.begin(), p.end());
	}
}


int main(int argc, char **argv) {
  
  CGAL::Timer timer;
  
	std::vector<Point_d> points;
	int n, d;
  
  read_input(std::cin, std::back_inserter(points), n, d);
  std::cout << "input read" << std::endl;
  
  std::vector<Point_d> extreme_points;
  timer.start();
  extreme_points_d_dula_helgason(points.begin(), points.end(), std::back_inserter(extreme_points));
  timer.stop();
  std::cout << "extreme_points_dula_helgason found " << extreme_points.size() << " extreme points" << std::endl;
  std::cout << "Time spent: " << timer.time() << std::endl;
  
	return 0;
}
