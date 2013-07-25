
// Author(s)     :

#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Segment_2 Segment;
typedef K::Point_2 Point;

typedef std::vector<Segment>::iterator Iterator;
typedef CGAL::AABB_segment_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

//types to return primitive results, template argument is the query type
typedef boost::optional< Tree::Intersection_and_primitive_id<Segment>::Type > Segment_intersection;

//types to return primitive results for point queries
typedef Tree::Point_and_primitive_id Point_and_primitive_id;


int main()
{
	Point a(1,0);
	Point b(0,1);
	Point c(-1,0);
	Point d(0,-1);

	std::vector<Segment> segments;
	segments.push_back(Segment(a,b));
	segments.push_back(Segment(b,c));
	segments.push_back(Segment(c,d));
	segments.push_back(Segment(d,a));


	Tree tree(segments.begin(),segments.end());

	// segment intersection query
	Point p((FT) -2,  (FT) -2);
	Point q((FT) 2,  (FT) 2);
	Segment pq(p,q);

	//retrieve an intersection from the AABB tree
	//The return type contain the intersection by segment query and primitive that has been intersected.
	
	Segment_intersection intersection = tree.any_intersection(pq);

	//The first element of the pair Intersection_and_primitive_id contains the intersection
	//result. In this example it is a point.
	Point *point;
	if( point = boost::get<Point>(&(intersection->first)))
        std::cout << "Intersection point " << *point << std::endl;
   
	
	//Primitive::Id is the iterator type provided while constructing the AABB tree.
	//Second element of the pair Intersection_and_primitive_id contain the id of the primitive,
	//which is basically the index of the provided primitive vector.

	Iterator index = intersection->second;

	//Get the integer index
	unsigned int pindex = std::distance(segments.begin(), index);

	std::cout<<"Intersected primitive index: "<<pindex<<std::endl;
	std::cout<<"Intersected primitive: "<<*index<<std::endl;

	// Same as above the closest point and corresponding primitive can be retrieved.

	Point r((FT)-1.0, (FT)-1.0);

	//retrieve the closest point and primitive
	boost::optional<Point_and_primitive_id> point_primitive;
	point_primitive = tree.closest_point_and_primitive(r);

	Point_and_primitive_id pp = *point_primitive;
	Point result = pp.first;

	std::cout<<"Closest point to "<<r<< " is "<<result<<std::endl;


	//index of the primitive containing closest point
	index = pp.second;
	pindex = std::distance(segments.begin(), index);

	std::cout<<"Closest primitive index: "<< pindex<<std::endl;
	std::cout<<"Closest primitive: "<<*index<<std::endl;

	return EXIT_SUCCESS;
}

