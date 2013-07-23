// Copyright (c)  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
//
// Author(s)     :
//
//******************************************************************************
// File Description :
//
//******************************************************************************

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

//types to return primitive results
typedef Tree::Object_and_primitive_id Object_and_primitive_id;
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
	//The return type contain the exact intersected object and primitive that has been intersected.
	//AABB Object_and_primitive_id is used to get the result. 
	boost::optional<Object_and_primitive_id> any;
	any = tree.any_intersection(pq);

	//First element of the std::pair<Object,typename Primitive::Id> contain the intersection result which is a
	//CGAL::Object, in this example intersection is a point. 
	Object_and_primitive_id op = *any;
	CGAL::Object object = op.first;
	Point point;

	//Assign the object to a point. Note that this should be checked if the intersection object
	//type is not known.
	CGAL::assign(point,object);

	std::cout << "Intersection point: " << point << std::endl;

	//Primitive::Id is the iterator type provided while constructing the AABB tree.
	//Second element of the std::pair<Object,typename Primitive::Id> contain the id of the primitive,
	//which is basically the index of the provided primitive vector.

	Iterator index = op.second;

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

