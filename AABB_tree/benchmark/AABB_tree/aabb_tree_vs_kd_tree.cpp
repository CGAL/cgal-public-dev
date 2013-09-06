
#include <iostream>
#include <vector>
#include <algorithm>
#include <exception>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_point_primitive.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

//search tree
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Timer.h>

typedef CGAL::Simple_cartesian<double> K;

//typedefs search tree
typedef K::Point_2 Point_d;
typedef CGAL::Search_traits_2<K> Traits_Search;
typedef CGAL::Kd_tree<Traits_Search> Tree1;
typedef CGAL::Fuzzy_sphere<Traits_Search> Fuzzy_sphere;

//tydefs AABB Tree
typedef K::FT FT;
typedef K::Point_2 Point;
typedef K::Circle_2 Circle;

typedef std::vector<Point>::iterator Iterator;
typedef CGAL::AABB_point_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

//tepedefs point generator
typedef CGAL::Creator_uniform_2<double,Point>  Creator;

int number_of_points, number_of_queries;

template<typename LIST>
void test_spatial_search(LIST &lst)
{
  CGAL::Timer timer;    
  timer.start();
  Tree1 search_tree(lst.begin(),lst.end());
  search_tree.build();
  timer.stop();
  std::cout << "building the KD tree took " << timer.time() << " sec." << std::endl;
  timer.reset();

  timer.start();
  std::size_t number_of_points = 0;

  for(int i=0;i<number_of_queries;i++)
    {
      LIST itr;
      Fuzzy_sphere fs(Point_d(0.05*i,0.05*i), 10.0, 0.0);


      search_tree.search(std::back_inserter(itr), fs);

      number_of_points+=itr.size();
      itr.clear();
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






template<typename LIST>
void test_aabb_tree(LIST &lst)
{
  CGAL::Timer timer;    
  timer.start();
  Tree aabb_tree(lst.begin(), lst.end());
  aabb_tree.build();
  timer.stop();
  std::cout << "building the AABB tree took " << timer.time() << " sec." << std::endl;
  timer.reset();

  timer.start();

  unsigned int number_of_range_queries = 0;
  std::size_t number_of_points = 0;

  for(int i=0;i<number_of_queries;i++)
    {
      std::vector<typename Primitive::Id> primitives;
      Circle circular_query(Point_d(0.05*i,0.05*i), 10.0*10.0);

      aabb_tree.all_contained_primitives(circular_query,std::back_inserter(primitives));

      number_of_points+=primitives.size();
      primitives.clear();
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

  std::cout << "Enter the number of points and the number of queries" << std::endl;
  std::cin >> number_of_points >> number_of_queries;
	
  std::vector<Point> points;

  CGAL::Random_points_in_disc_2<Point,Creator> g( 150.0);
  try
    {
      CGAL::cpp11::copy_n( g, number_of_points, std::back_inserter(points));
      std::cout << "Created " << number_of_points << " points" << std::endl;
    }
  catch (std::bad_alloc& )
    {
      std::cerr << "Can not create " << number_of_points << " points" << std::endl;
    }
		


  test_spatial_search(points);

  test_aabb_tree(points);


   	
  std::cerr << "done" << std::endl;
  return EXIT_SUCCESS;
}
