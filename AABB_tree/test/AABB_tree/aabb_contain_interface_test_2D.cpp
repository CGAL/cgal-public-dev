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
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>

typedef CGAL::Simple_cartesian<double>  K;

template<typename K, typename List>
bool load_segment_set(List &list)
{
	typedef typename K::Point_2 Point;
	typedef typename K::Segment_2 Segment;

	Point a(1,0);
	Point b(0,1);
	Point c(-1,0);
	Point d(0,-1);

	list.push_back(Segment(a,b));
	list.push_back(Segment(b,c));
	list.push_back(Segment(c,d));
	list.push_back(Segment(d,a));

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
	typedef typename K::Iso_rectangle_2 Rectangle;

	// load simple segment set
	typedef std::list<Segment> Segments;
	Segments segments;
	load_segment_set<K,Segments>(segments);

	// construct tree from segments
	typedef typename Segments::iterator Iterator;
	typedef typename CGAL::AABB_segment_primitive<K, Iterator> Primitive;
	typedef typename CGAL::AABB_traits<K,Primitive> Traits;
	typedef typename CGAL::AABB_tree<Traits> Tree;
	typedef boost::optional< Tree::Intersection_and_primitive_id<Segment>::Type > Segment_intersection;
	typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;
	Tree tree(segments.begin(),segments.end());

	// rectangular range query
	Rectangle rectangular_query(-1.0,0.0,1.0,1.0);

	bool contain = tree.do_contain(rectangular_query);
	
	if(!contain)
	{
		std::cerr << "contain check error" << std::endl;
		return EXIT_FAILURE;
	}
	
	// number of contained primitives
	std::size_t count = tree.number_of_contained_primitives(rectangular_query);
	if(count != 2)
	{
		std::cerr << "wrong number of contained primitives" << std::endl;
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

	if(test<CGAL::Exact_predicates_inexact_constructions_kernel>() == EXIT_FAILURE)
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}

