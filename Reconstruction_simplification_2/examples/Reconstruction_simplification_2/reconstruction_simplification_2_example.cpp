#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>

#include <fstream>
// Example for Reconstruction_simplification_2

#include <iostream>
#include <string>
#include <iterator>
#include <utility>      // std::pair
#include <list>

#include <CGAL/property_map.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point;
typedef K::Segment_2                                        Segment;

typedef K::FT                                               FT;

typedef std::pair<Point, FT>                                PointMassPair;
typedef std::list<PointMassPair>                            PointMassList;

typedef CGAL::First_of_pair_property_map <PointMassPair>    Point_property_map;
typedef CGAL::Second_of_pair_property_map <PointMassPair>   Mass_property_map;

typedef CGAL::Reconstruction_simplification_2<
  K, Point_property_map, Mass_property_map>                 Rs_2;


void load_xym_file(const std::string& filename, PointMassList& points)
{
   std::ifstream ifs(filename);
	
   Point point;
	 FT mass;

   while (ifs >> point && ifs >> mass)
	   points.push_back(std::make_pair(point, mass));
   
   ifs.close();
}

int main ()
{
  PointMassList points;
  
  load_xym_file("data/stair.xym", points);
  
  Point_property_map point_pmap;
  Mass_property_map  mass_pmap;
  
  Rs_2 rs2(points, point_pmap, mass_pmap);
  
  rs2.run(100); // 100 steps
  
  std::vector<Point> isolated_vertices;
  std::vector<Segment> edges;
  
  rs2.extract_list_output(std::back_inserter(isolated_vertices), std::back_inserter(edges));
  
  std::cerr << "Isolated vertices" << std::endl;
  std::vector<Point>::iterator vit;
  for (vit = isolated_vertices.begin(); vit != isolated_vertices.end(); vit++) 
	std::cout  <<  *vit << std::endl;
  
  std::cerr << "Edges" << std::endl;
  std::vector<Segment>::iterator eit;
  for (eit = edges.begin(); eit != edges.end(); eit++) 
	std::cout << *eit << std::endl;
  
  return 0;
}
