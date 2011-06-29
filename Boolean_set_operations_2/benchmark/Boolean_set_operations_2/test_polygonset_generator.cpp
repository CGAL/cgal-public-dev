//Program to generate random polygons and write to a file 

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Polygon_2.h>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2                                 Point_2;
typedef std::list<Point_2>                         Container;
typedef CGAL::Polygon_2<K, Container>              Polygon_2;
typedef CGAL::Random_points_in_square_2< Point_2 > Point_generator;

int main() {
   
  //output file stream
  std::ofstream myfile;
  
  myfile.open ("test_8");

  //number of polygons to be generated
  int num_of_polygons = 30000;
  int count = 0;
  
  //first line of the file contains number of polygons
  myfile << num_of_polygons << std::endl;
  
  while(count < num_of_polygons)
  {
	  Polygon_2 polygon;
	  //first argument is the number of vertices in the polygon
	  CGAL::random_polygon_2(8, std::back_inserter(polygon),
							 Point_generator(100));
							 
	  //write polygon to file. First number is the number of vertices in the polygon, followed by pairs of coordinates
	  myfile << polygon;
	  myfile << std::endl;
	  count++;

  }
  
  myfile.close();
  return 0;
}
