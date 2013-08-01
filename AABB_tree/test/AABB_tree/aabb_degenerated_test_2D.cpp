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

#include <CGAL/Timer.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>

typedef CGAL::Simple_cartesian<double>  K;

template<typename K, typename LIST>
bool load_segment_set(LIST &list)
{
	typedef typename K::Point_2 Point;
	typedef typename K::Segment_2 Segment;

	Point a(1,0);
	Point b(0,1);
	Point c(-1,0);
	
	list.push_back(Segment(a,b));
	list.push_back(Segment(b,c));
	list.push_back(Segment(c,b));
	list.push_back(Segment(b,a));

	return true;
}

template <typename K>
int test()
{
	// types
	typedef typename K::FT FT;
	typedef typename K::Line_2 Line;
	typedef typename K::Point_2 Point;
	typedef typename K::Segment_2 Segment;
	typedef typename K::Triangle_2 Triangle;

	// load simple segment set
	typedef std::list<Segment> SegmentList;
	SegmentList segments;
	load_segment_set<K,SegmentList>(segments);

	// construct tree from segments
	typedef typename SegmentList::iterator Iterator;
	typedef typename CGAL::AABB_segment_primitive<K, Iterator> Primitive;
	typedef typename CGAL::AABB_traits<K,Primitive> Traits;
	typedef typename CGAL::AABB_tree<Traits> Tree;
	typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;
	typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;
	Tree tree(segments.begin(),segments.end());

	// degenerated segment intersection query
	Point p((FT) 0.0,  (FT) 1.0);
	Point q((FT) 0.0,  (FT) 1.0);
	Segment pq(p,q);

	if(!tree.do_intersect(pq))
	{
		std::cerr << "no intersection found" << std::endl;
		return EXIT_FAILURE;
	}

	if(tree.number_of_intersected_primitives(pq) != 4)
	{
		std::cerr << "number of intersections different than two" << std::endl;
		return EXIT_FAILURE;
	}

	boost::optional<Object_and_primitive_id> any;
	any = tree.any_intersection(pq);
	if(!any)
	{
		std::cerr << "did not find any intersection" << std::endl;
		return EXIT_FAILURE;
	}

	Object_and_primitive_id op = *any;
	CGAL::Object object = op.first;
	Point point;
	if(CGAL::assign(point,object))
	{
		std::cout << "Intersection point: " << point << std::endl;
	}
	else
	{
		std::cerr << "Intersection does not assign to a point" << std::endl;
		return EXIT_FAILURE;
	}

	// closest point query
	Point r((FT)0.0, (FT)1.0);
	Point closest((FT)0.0, (FT)1.0);
	Point result = tree.closest_point(r);
	if(result != closest)
	{
		std::cerr << "wrong closest point" << std::endl;
		return EXIT_FAILURE;
	}

	Point hint((FT)0.399, (FT)0.601);
	result = tree.closest_point(r,hint);
	if(result != closest)
	{
		std::cerr << "wrong closest point" << std::endl;
		return EXIT_FAILURE;
	}

	FT distance_known = 0.0;
	FT distance_calculated = tree.squared_distance(r);
	if(distance_calculated != distance_known)
	{
		std::cerr << "wrong distance" << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

int main()
{
	if(test<CGAL::Simple_cartesian<float> >() == EXIT_FAILURE)
		return EXIT_FAILURE;

	if(test<CGAL::Simple_cartesian<double> >() == EXIT_FAILURE)
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}
