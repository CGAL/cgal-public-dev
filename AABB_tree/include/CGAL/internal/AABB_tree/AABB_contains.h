/*
 * Do_contain_traits.h
 *
 *  Created on: 19 août 2013
 *      Author: jayalatn
 */

#ifndef Do_contain_traits_H_
#define Do_contain_traits_H_

#include "Get_primitve_vertice_size.h"

template<typename AABBTraits, typename GeomTraits,typename PrimitiveType>
class Do_contain_traits
{
public:
	Do_contain_traits():number_of_vertices(Get_primitive_vertice_size<GeomTraits,PrimitiveType>()()){};

	template<typename Query>
	bool operator() (Query &query, PrimitiveType &primitive) const
	{
		int itr = number_of_vertices;

		while(itr>0)
		{
			AABBTraits::Point point = AABBTraits::Construct_vertex_d()(primitive,itr);

			//if(!Query.contains(point))
			//	return false;
			itr++;
		}

		return true;
	}


private:
	unsigned int number_of_vertices;
};


//consider later


//template<typename GeomTraits>
//class Do_contain_traits<GeomTraits,GeomTraits::Sphere_3>
//{
//public:
//	Do_contain_traits():number_of_vertices(Get_primitive_vertice_size<GeomTraits,Primitive_type>()()){};
//
//
//
//private:
//	unsigned int number_of_vertices;
//};
//
//template<typename GeomTraits>
//class Do_contain_traits<GeomTraits,GeomTraits::Circle_2>
//{
//public:
//	Do_contain_traits():number_of_vertices(Get_primitive_vertice_size<GeomTraits,Primitive_type>()()){};
//
//
//
//private:
//	unsigned int number_of_vertices;
//};

#endif /* Do_contain_traits_H_ */
