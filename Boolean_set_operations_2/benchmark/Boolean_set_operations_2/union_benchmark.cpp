/*! \file union_benchmark.cpp
 * Computing the union of a set of polygons read from a file.
 */

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <cstdlib>
#include <time.h>

//to enable user to enable or disable OpenMP only for Boolean_set_operations_2 package
#ifndef CGAL_BOOLEAN_SET_OPERATIONS_2_DONT_USE_OPENMP
#define CGAL_BOOLEAN_SET_OPERATIONS_2_USE_OPENMP 
#endif

//add OpenMP header
#if defined _OPENMP && defined CGAL_BOOLEAN_SET_OPERATIONS_2_USE_OPENMP
#include<omp.h>
#endif

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2                                 Point_2;
typedef std::vector<Point_2>                       Container;
typedef CGAL::Polygon_2<K, Container>              Polygon_2;
typedef CGAL::Arr_segment_traits_2<K>			   ArrSegmentTraits;
typedef CGAL::Gps_segment_traits_2<K,Container,ArrSegmentTraits> Traits_2;
typedef CGAL::General_polygon_set_2<Traits_2>      General_polygon_set_2;


static const int DEFAULT_GROUP_SIZE = 5;

int main() {
	// Open the input file.
	const char* filename = "test_4.txt";
	std::ifstream input_file(filename);
	if (! input_file.is_open())
	{
		std::cerr << "Failed to open the " << filename <<std::endl;
		return -1;
	}

	int     group_size = DEFAULT_GROUP_SIZE;
	std::vector<Polygon_2>	pgns;
	General_polygon_set_2   gps;
	
	// Read the polygons from the file.
	
	CGAL::Timer                    t_read;
      
	std::cout << "Reading <" << filename << "> ... " << std::flush;
	t_read.start();
       	
	int number_of_polygons;
	
	input_file >> number_of_polygons;
	
	int count = 0;
	
	while(count < number_of_polygons)
	{
		int number_of_vertices;
		
		input_file >> number_of_vertices;
		
		int vertex_count = 0;

		//create a polygon from the points read from the file
		Polygon_2 p;
		
		while(vertex_count < number_of_vertices)
		{
			  float p1,p2;
			  
			  input_file >> p1 >> p2;			  
			  p.push_back(Point_2(p1,p2));
			  	  
			  vertex_count++;
		}

		//Insert to the vector
		pgns.push_back(p);
		count++;
	}
	
	t_read.stop();
           
	std::cout << "Done! (" << t_read.time() << " seconds)." << std::endl;
	std::cout << std::distance (pgns.begin(), pgns.end()) << " polygons, " << std::endl;
    
   	input_file.close();
	
	// Compute their union.
	CGAL::Timer                    t_union;
 	time_t                         s,e;

	std::cout << "Computing the union ... " << std::flush;

	t_union.start();
     	s = time(NULL);

        #if defined _OPENMP && defined CGAL_BOOLEAN_SET_OPERATIONS_2_USE_OPENMP
        double start,end;
        start = omp_get_wtime(); 
        #endif     
	
	gps.join (pgns.begin(), pgns.end());
 
        #if defined _OPENMP && defined CGAL_BOOLEAN_SET_OPERATIONS_2_USE_OPENMP
        end = omp_get_wtime();
        #endif

	t_union.stop();
        e = time(NULL);

        #if defined _OPENMP && defined CGAL_BOOLEAN_SET_OPERATIONS_2_USE_OPENMP
        std::cout << "Done! time using omp_get_wtime() (" << end-start << " seconds)." << std::endl;
        #endif 
  
	std::cout << "Done! time using CGAL::Timer (using measure of clock cycles) (" << t_union.time() << " seconds).\n";

        std::cout << "Done! time using time() (" << difftime(e,s) << " seconds)." << std::endl;

	std::cout << "The result:"
			  << "  |V| = " << gps.arrangement().number_of_vertices()
			  << "  |E| = " << gps.arrangement().number_of_edges()
			  << "  |F| = " << gps.arrangement().number_of_faces() << std::endl;
	 
  return 0;
}
