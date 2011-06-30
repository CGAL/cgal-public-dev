//! \file compute_intersection_points.cpp
// Computing intersection points among set of polygons using the sweep line.

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <list>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef Kernel::Point_2                                 Point_2;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
typedef Traits::Curve_2									Segment_2;

int MAX_VERTICES = 4;

int main()
{
	// Open the input file.
	const char* filename = "test_4.txt";
	std::ifstream input_file(filename);
	if (! input_file.is_open())
	{
		std::cerr << "Failed to open the " << filename <<std::endl;
		return -1;
	}

	//resultant list of intersecting points
	std::list<Point_2>     	 intersecting_pts;	
	CGAL::Timer              t_read;
	  
	std::cout << "Reading <" << filename << "> ... " << std::flush;
	t_read.start();
	
	int number_of_polygons;
	
	input_file >> number_of_polygons;

	//array to hold segments in the test file
	Segment_2 * segments;

	//multiply number of polygons with max number of vertices in any polygon in the test data
	segments = new  Segment_2[number_of_polygons * MAX_VERTICES];

	//count number of segments
	int seg_count = 0;
	
	int count = 0;
	
	while(count < number_of_polygons)
	{
		int number_of_vertices;
		
		input_file >> number_of_vertices;
		
		int vertex_count = 0;

		//Array of points of the current polygon
		Point_2 * pts;
		pts = new Point_2[number_of_vertices];
	
		while(vertex_count < number_of_vertices)
		{
			  float p1,p2;
			  
			  input_file >> p1 >> p2;			  
			  Point_2 pt(p1,p2);
			  //store the point in the array
			  pts[vertex_count] = pt;
			  vertex_count++;
		}
		
		//Add line segments created by the points
		for(int i = 0 ; i < number_of_vertices ; i++)
		{
			Segment_2 seg(pts[i],pts[(i+1) % number_of_vertices]);
			segments[seg_count++] = seg;
		}

		count++;
	}
	
	t_read.stop();
	std::cout << "Done! (" << t_read.time() << " seconds)." << std::endl;
		
	input_file.close();
	
    CGAL::compute_intersection_points (segments, segments + seg_count,
                                     std::back_inserter (intersecting_pts));
  
	// Print the result.
	std::cout << "Found " << intersecting_pts.size() << " intersection points: " << std::endl; 
	
    return 0;
}
