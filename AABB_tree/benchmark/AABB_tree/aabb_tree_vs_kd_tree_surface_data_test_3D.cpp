// Author : 

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>  
#include <exception>
#include <fstream>

#include <CGAL/IO/read_off_points.h>
#include <utility> 

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_point_primitive.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/random_selection.h>

//search tree
#include <CGAL/Cartesian_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Timer.h>

typedef CGAL::Simple_cartesian<double> K;

//typedefs search tree
typedef K::Point_3 Point_d;
typedef CGAL::Search_traits_3<K> Traits_Search;
typedef CGAL::Kd_tree<Traits_Search> Tree1;
typedef CGAL::Fuzzy_sphere<Traits_Search> Fuzzy_sphere;

//tydefs AABB Tree
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Sphere_3 Circle;

typedef std::vector<Point>::iterator Iterator;
typedef CGAL::AABB_point_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;


template<typename vector>
void test_spatial_search(vector &lst,vector &centers)
{
    CGAL::Timer timer;
	timer.start();
	Tree1 search_tree(lst.begin(),lst.end());
	search_tree.build();
	timer.stop();
	std::cout << "building the KD tree took " << timer.time() << " sec." << std::endl;
	timer.reset();
	
	lst.clear();

	std::size_t number_of_points = 0;


	timer.start();

	vector::iterator itr = centers.begin();
	while(itr!=centers.end())
	{
	   vector primitives;
	   Fuzzy_sphere fs(*itr, 5.0, 0.0);

	   
	   search_tree.search(std::back_inserter(primitives), fs);
	  
	   number_of_points+=primitives.size();
	   primitives.clear();
	   ++itr;
	}
	timer.stop();
	std::cout << "KD Total Points Searched: " << number_of_points << std::endl;
	std::cout << "Searching took " << timer.time() << " sec." << std::endl;

	timer.reset();
	timer.start();
	search_tree.clear();
	timer.stop();
	std::cout << "Clearing took " << timer.time() << " sec.\n" << std::endl;

}






template<typename vector>
void test_aabb_tree(vector &lst,vector &centers)
{
	CGAL::Timer timer;    
	timer.start();
	Tree aabb_tree(lst.begin(), lst.end());
	aabb_tree.build();
	timer.stop();
	std::cout << "building the AABB tree took " << timer.time() << " sec." << std::endl;
	timer.reset();

 
	std::size_t number_of_points = 0;
	timer.start();

	vector::iterator itr = centers.begin();

	while(itr!=centers.end())
	{
		std::vector<typename Primitive::Id> primitives;
	   Circle circular_query(*itr, 5.0*5.0);
	
	   aabb_tree.all_contained_primitives(circular_query,std::back_inserter(primitives));
	   
	   number_of_points+=primitives.size();
	   primitives.clear();
		++itr;
	}

	timer.stop();
  
  std::cout << "AABB Total Points Searched: " << number_of_points << std::endl;
  std::cout << "Searching took " << timer.time() << " sec." << std::endl;
  timer.reset();
  timer.start();
  aabb_tree.clear();
  timer.stop();
  std::cout << "Clearing took " << timer.time() << " sec.\n" << std::endl;
}




int main()
{
   	
	
 
	//change the number of points in 10^n
	std::vector<Point> points;
	std::vector<Point> search_centers;

	std::ifstream in("./data/gargo1.7M-4.off");
	CGAL::read_off_points(in,std::back_inserter(points));

	std::size_t number_of_points = points.size();


	for(int i=0;i<10000;i++)
	{
		int iPos = i*number_of_points/10000;
		search_centers.push_back(points.at(iPos));
	}


	test_aabb_tree<std::vector<Point>>(points,search_centers);
	test_spatial_search<std::vector<Point>>(points,search_centers);


	points.clear();
	
   	
	return EXIT_SUCCESS;
}
