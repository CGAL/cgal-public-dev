// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
#include <CGAL/AABB_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double>  K;

template<typename K, typename LIST>
bool load_triangle_set(LIST &list)
{
	typedef typename K::Point_2 Point;
	typedef typename K::Triangle_2 Triangle;

	Point a(0,0);
	Point b(1,-2);
	Point c(2,-1);
	Point d(2,1);
	Point e(1,2);
	Point f(-1,-2);
	Point g(-2,1);
	Point h(-2,-1);
	Point i(-1,-2);

	list.push_back(Triangle(a,b,c));
	list.push_back(Triangle(a,d,e));
	list.push_back(Triangle(a,f,g));
	list.push_back(Triangle(a,h,i));

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

	// load simple triangle set
	typedef std::list<Triangle> List_Triangle;
	List_Triangle triangles;
	load_triangle_set<K,List_Triangle>(triangles);

	// construct tree from triangles
	typedef typename List_Triangle::iterator Iterator;
	typedef typename CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
	typedef typename CGAL::AABB_traits<K,Primitive> Traits;
	typedef typename CGAL::AABB_tree<Traits> Tree;
	typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;
	Tree tree(triangles.begin(),triangles.end());

	// segment intersection query
	Point p((FT) 0.25,  (FT) 0.0);
	Point q((FT) 0.25,  (FT) 1.255);
	Segment pq(p,q);

	if(!tree.do_intersect(pq))
	{
		std::cerr << "no intersection found" << std::endl;
		return EXIT_FAILURE;
	}

	if(tree.number_of_intersected_primitives(pq) != 1)
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
	Segment segment;
	if(CGAL::assign(segment,object))
	{
		std::cout << "Intersection segment: " << segment << std::endl;
	}
	else
	{
		std::cerr << "Intersection does not assign to a segment" << std::endl;
		return EXIT_FAILURE;
	}

	//check the iterator 
	Iterator index = op.second;

	if(std::distance(triangles.begin(), index)!=1)
	{
	  std::cerr << "Iterator gives a wrong index" << std::endl;
	  return EXIT_FAILURE;
	}

	// closest point query
	Point r((FT)5.0, (FT)1.0);
	Point closest((FT)2.0, (FT)1.0);
	Point result = tree.closest_point(r);
	if(result != closest)
	{
		std::cerr << "wrong closest point" << std::endl;
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

