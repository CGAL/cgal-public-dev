#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <algorithm>
#include <iterator>

typedef CGAL::Simple_cartesian<CGAL::Gmpq> Kernel;
typedef CGAL::Delaunay_triangulation_3<Kernel> DT3;
typedef CGAL::Creator_uniform_3<CGAL::Gmpq,Kernel::Point_3>  Creator;


int main(int argc,char** argv)
{
  std::vector<Kernel::Point_3> points;
  CGAL::Random rand(0);
  CGAL::Random_points_in_sphere_3<Kernel::Point_3,Creator> g(1, rand);
  
  int nb_pts = 1000;
  if (argc!=1)
    nb_pts = atoi(argv[1]);
  
  std::cout << "Using " << nb_pts << " input points\n";
  
  
  CGAL::cpp11::copy_n( g, nb_pts, std::back_inserter(points));  
  
  CGAL::Timer timer;
  
  timer.start();
  DT3 dt3(points.begin(),points.end());
  timer.stop();
  std::cout << "Insert "<<  dt3.number_of_vertices() << " pts in "
            << timer.time() << "s.\n";
  
}
