/*
 * Do_contain_traits.h
 *
 *  Created on: 19 août 2013
 *      Author: jayalatn
 */

#ifndef CGAL_DO_CONTAIN_TEST_TRAITS_H_
#define CGAL_DO_CONTAIN_TEST_TRAITS_H_

#include "Get_primitve_vertice_count.h"

using namespace internal;

namespace CGAL

{

template<typename AABBTraits, typename GeomTraits, typename ObjectType>
class Do_contain_test_traits
{
public:
	Do_contain_test_traits(){
		number_of_vertices = Get_primitive_vertice_count<GeomTraits,ObjectType>()();
	}

	template<typename Query>
	bool operator() (Query &query, ObjectType &object) const
	{
		unsigned int count = number_of_vertices;
		while(count>0)
		{
			AABBTraits::Point point = AABBTraits::Construct_vertex_d()(object,count);
			if(query.has_on_unbounded_side(point))
				return false;
			count--;
		}

		return true;
	}


private:
	unsigned int number_of_vertices;
};


template<typename AABBTraits, typename GeomTraits>
class Do_contain_test_traits<AABBTraits,GeomTraits,typename GeomTraits::Point_2>
{
public:
	Do_contain_test_traits(){
		number_of_vertices = Get_primitive_vertice_count<GeomTraits,GeomTraits::Point_2>()();
	}

	template<typename Query>
	bool operator() (Query &query, typename GeomTraits::Point_2 &object) const
	{
		unsigned int count = number_of_vertices;
		while(count>0)
		{
			AABBTraits::Point point = object;
			if(query.has_on_unbounded_side(point))
				return false;
			count--;
		}

		return true;
	}


private:
	unsigned int number_of_vertices;
};


//For sphere and circle


template<typename AABBTraits,typename GeomTraits>
class Do_contain_test_traits<AABBTraits,GeomTraits,typename GeomTraits::Sphere_3>
{
public:
	Do_contain_test_traits():number_of_vertices(Get_primitive_vertice_count<GeomTraits,Primitive_type>()()){}


private:
	unsigned int number_of_vertices;
};

template<typename AABBTraits,typename GeomTraits>
class Do_contain_test_traits<AABBTraits,GeomTraits,typename GeomTraits::Circle_2>
{
public:
	Do_contain_test_traits():number_of_vertices(Get_primitive_vertice_count<GeomTraits,Primitive_type>()()){}


private:
	unsigned int number_of_vertices;
};

}
#endif /* Do_contain_test_traits_H_ */
