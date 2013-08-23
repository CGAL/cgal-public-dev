// Author : 

#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/Bbox_3.h>

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Segment_3 Segment;
typedef K::Point_3 Point;
typedef CGAL::Bbox_3 Box;

typedef std::list<Segment>::iterator Iterator;
typedef CGAL::AABB_segment_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;


int main()
{
	Point a(0.0, 0.0, 0.0);
	Point b(2.0, 1.0, 1.0);
	Point c(3.0, 4.0, 2.0);
	Point d(1.0, 6.0, 3.0);
	Point e(0.0, 3.0, 1.0);

	std::list<Segment> segments;
	segments.push_back(Segment(a, b));
	segments.push_back(Segment(b, c));
	segments.push_back(Segment(c, d));
	segments.push_back(Segment(d, e));
	segments.push_back(Segment(e, a));

	//   // constructs the AABB tree 
	Tree tree(segments.begin(), segments.end());
	
	//   // construct the rectangle query object
	Box rectangle_query(0.0,0.0,0.0,3.0,4.0,2.0);

	//   // Check if the query contains primitives
	bool contain_check = tree.do_contain(rectangle_query);

	if(contain_check)
		std::cout << "The rectangular range fully contains at least one primitive" << std::endl;


	//   // Count the number of primitives fully contained
	std::cout << tree.number_of_contained_primitives(rectangle_query)
        << " primitives are fully contained inside the rectangular range" << std::endl;

	return EXIT_SUCCESS;
}
