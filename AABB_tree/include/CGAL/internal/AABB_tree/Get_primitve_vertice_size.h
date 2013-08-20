/*
 * Get_primitve_vertice_size.h
 *
 *  Created on: 19 août 2013
 *      Author: jayalatn
 */

#ifndef GET_PRIMITVE_VERTICE_SIZE_H_
#define GET_PRIMITVE_VERTICE_SIZE_H_

template<typename K,typename PrimitiveType>
class Get_primitive_vertice_size
{

};

template<typename K>
class Get_primitive_vertice_size<K, typename K::Segment_2>
{
public:
   Get_primitive_vertice_size():number_of_vertices(2){};
   unsigned int operator()() const
   {

	   return number_of_vertices;
   }
private:
   unsigned int number_of_vertices;

};

template<typename K>
class Get_primitive_vertice_size<K, typename K::Triangle_2>
{
public:
   Get_primitive_vertice_size():number_of_vertices(3){};
   unsigned int operator()() const
   {

	   return number_of_vertices;
   }
private:
   unsigned int number_of_vertices;
};

template<typename K>
class Get_primitive_vertice_size<K, typename K::Segment_3>
{
public:
   Get_primitive_vertice_size():number_of_vertices(2){};
   unsigned int operator()() const
   {

	   return number_of_vertices;
   }
private:
   unsigned int number_of_vertices;

};

template<typename K>
class Get_primitive_vertice_size<K, typename K::Triangle_3>
{
public:
   Get_primitive_vertice_size():number_of_vertices(3){};
   unsigned int operator()() const
   {

	   return number_of_vertices;
   }
private:
   unsigned int number_of_vertices;

};


#endif /* GET_PRIMITVE_VERTICE_SIZE_H_ */
