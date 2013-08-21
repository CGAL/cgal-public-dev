/*
 * Get_primitve_vertice_size.h
 *
 *  Created on: 19 août 2013
 *      Author: jayalatn
 */

#ifndef CGAL_GET_PRIMITVE_VERTICE_COUNT_H_
#define CGAL_GET_PRIMITVE_VERTICE_COUNT_H_

class Get_primitive_vertice_count_base
{
public:
	Get_primitive_vertice_count_base(unsigned int count):number_of_vertices(count){};
   
	inline unsigned int operator()() const
   {
	   return number_of_vertices;
   }
private:
   unsigned int number_of_vertices;
};

template<typename K,typename PrimitiveType>
class Get_primitive_vertice_count
{

};

template<typename K>
class Get_primitive_vertice_count<K, typename K::Segment_2>:public Get_primitive_vertice_count
{
public:
   Get_primitive_vertice_count():Get_primitive_vertice_count(2){};
 };

template<typename K>
class Get_primitive_vertice_count<K, typename K::Triangle_2>:public Get_primitive_vertice_count
{
public:
   Get_primitive_vertice_count():Get_primitive_vertice_count(3){};
};

template<typename K>
class Get_primitive_vertice_count<K, typename K::Segment_3>:public Get_primitive_vertice_count
{
public:
   Get_primitive_vertice_count():Get_primitive_vertice_count(2){};
};

template<typename K>
class Get_primitive_vertice_count<K, typename K::Triangle_3>:public Get_primitive_vertice_count
{
public:
   Get_primitive_vertice_count():Get_primitive_vertice_count(3){};
};


#endif /* GET_PRIMITVE_VERTICE_SIZE_H_ */
