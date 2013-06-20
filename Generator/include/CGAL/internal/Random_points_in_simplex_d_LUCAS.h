#include <iostream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <iterator>
#include <CGAL/enum.h>


typedef CGAL::Simple_cartesian<double>      K;
typedef K::Point_2                     Point_2;
typedef std::vector<Point_2>         Container;
typedef CGAL::Polygon_2<K, Container>   Polygon_2;
typedef K::Triangle_2               Triangle_2;

std::vector<double> generator(int d, CGAL::Random &rand) {
   std::vector<double> a;
   a.push_back(0.0);
   for(int i = 0; i < d; ++i) {
      a.push_back(rand.get_double(0,1));
   }
   a.push_back(1.0);
   std::sort(a.begin(),a.end());
   std::vector<double> b;
   for(int i = 0; i <= d; ++i) {
      b.push_back(a[i+1]-a[i]);
   }
   return b;
}

Point_2 operator*(const Point_2 &p, double c) {
   return CGAL::ORIGIN+((p-CGAL::ORIGIN)*c);
}

Point_2 operator+(const Point_2 &p, const Point_2 &q) {
   return Point_2(p.x()+q.x(),p.y()+q.y());
}

int main() {
   CGAL::Random rand;
   
   Point_2 pts[3] = {Point_2(50.5,60.7),Point_2(90.0,45.4),Point_2(42.0,34.5)};
   
   Point_2 int_points[3];
   
   for(int i = 0; i < 3; ++i) {
      int_points[i] = pts[i]+(pts[(i+1)%3]-pts[i])/2;
   }
   
   Triangle_2 out_triangle(pts[0],pts[1],pts[2]);
   Triangle_2 int_triangle(int_points[0],int_points[1],int_points[2]);

   
   int inside_small_triangle = 0;
   int total = 1000000;
   int summary[3]={0,0,0};
   for(int i = 0; i < total; ++i) {
      std::vector<double> random = generator(2,rand);
      Point_2 p = out_triangle[0]*random[0]+out_triangle[1]*random[1]+out_triangle[2]*random[2];
      switch(int_triangle.bounded_side(p)) {
         case CGAL::ON_BOUNDED_SIDE:
            ++inside_small_triangle;
            ++summary[0];
         break;
         
         case CGAL::ON_BOUNDARY:
            ++inside_small_triangle;
            ++summary[1];
         break;
         
         case CGAL::ON_UNBOUNDED_SIDE:
            ++summary[2];
         break;
      }   
      assert(out_triangle.bounded_side(p)==CGAL::ON_BOUNDED_SIDE ||
               out_triangle.bounded_side(p)==CGAL::ON_BOUNDARY);
   }
   std::cout<<"Points inside "<<summary[0]<<std::endl;
   std::cout<<"Points on boundary "<<summary[1]<<std::endl;
   std::cout<<"Points outside "<<summary[2]<<std::endl;
   std::cout<<"Fraction of points inside smaller triangle: "<<double(inside_small_triangle)/double(total)<<std::endl;
   std::cout<<"Expected value: "<<fabs(int_triangle.area())/fabs(out_triangle.area())<<std::endl;
   return 0;
}