/*
 * AABB_test_util_2D.h
 *
 *  Created on: 21 juil. 2013
 *      Author: Mehdious
 */

#ifndef AABB_TEST_UTIL_2D_H_
#define AABB_TEST_UTIL_2D_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <boost/mem_fn.hpp>
#include <CGAL/AABB_segment_primitive.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>


#define X_EXTENT 100
#define Y_EXTENT 100

enum Primitive_type {
  SEGMENT, TRIANGLE
};

enum Query_type {
	SEGMENT_QUERY=0, RAY_QUERY, LINE_QUERY, TRIANGLE_QUERY
};

const char Query_name[3][10]={
	"segment", "ray", "line"
};


double random_in(const double a,
                 const double b)
{
    double r = rand() / (double)RAND_MAX;
    return a + (b - a) * r;
}

template <typename K>
typename K::Point_2 random_point_in(const CGAL::Bbox_2& bbox)
{
    typedef typename K::FT FT;
    FT x = (FT)random_in(bbox.xmin(),bbox.xmax());
    FT y = (FT)random_in(bbox.ymin(),bbox.ymax());
    return typename K::Point_2(x,y);
}

template <typename K>
typename K::Vector_2 random_vector()
{
    typedef typename K::FT FT;
    FT x = (FT)random_in(0.0,1.0);
    FT y = (FT)random_in(0.0,1.0);
    return typename K::Vector_2(x,y);
}

template<typename K, typename Primitives>
void get_segment_set(const char *fileName, Primitives& primitives)
{
	std::ifstream primitiveFile(fileName);

	unsigned int size;

	primitiveFile >> size;
	std::cout<<size<<std::endl;

	while(size>0)
	{
		K::Segment_2 segment;
		primitiveFile >> segment;
		primitives.push_back(segment);
		size--;
	}

	std::cout<<"number"<<primitives.size()<<std::endl;
}

template <typename Tree, typename K>
void test_distance_speed(Tree& tree,
                         const double duration)
{
    typedef typename K::Point_2 Point;

    std::cout<<duration<<std::endl;
    CGAL::Timer timer;
    timer.start();
    unsigned int nb = 0;
    while(timer.time() < duration)
    {
            // picks a random point in the tree bbox
            Point query = random_point_in<K>(tree.bbox());
            std::cout<<tree.bbox()<<std::endl;
            Point closest = tree.closest_point(query);
            (void) closest;
            nb++;
    }
    double speed = (double)nb / timer.time();
    std::cout << speed << " distance queries/s" << std::endl;
    timer.stop();
}


template <typename Tree, typename K>
void test_do_intersection_speed(Tree& tree,
        const double duration, Query_type query_type)
{
	typedef typename K::Ray_2 Ray;
	typedef typename K::Line_2 Line;
	typedef typename K::Point_2 Point;
	typedef typename K::Vector_2 Vector;
	typedef typename K::Segment_2 Segment;

	CGAL::Timer timer;
	unsigned int nb = 0;
	timer.start();
	while(timer.time() < duration)
	{
		switch(query_type)
		{
			case RAY_QUERY:
				{
					Point source = random_point_in<K>(tree.bbox());
					Vector vec = random_vector<K>();
					tree.do_intersect(Ray(source, vec));
					break;
				}
			case SEGMENT_QUERY:
				{
					Point a = random_point_in<K>(tree.bbox());
					Point b = random_point_in<K>(tree.bbox());
					tree.do_intersect(Segment(a,b));
					break;
				}
				break;
			case LINE_QUERY:
				{
					Point a = random_point_in<K>(tree.bbox());
					Point b = random_point_in<K>(tree.bbox());
					tree.do_intersect(Line(a,b));
					break;
				}
		}
		nb++;
	}
	unsigned int speed = (unsigned int)(nb / timer.time());
	std::cout.precision(10);
	std::cout.width(15);
	std::cout << speed << " intersections/s with " << Query_name[query_type] << std::endl;
	timer.stop();
}



template<typename K, typename Tree, typename Primitive>
void test_impl(Tree& tree, const double duration);

template<typename K, Primitive_type Primitive>
struct Primitive_triats{};

template<typename K>
struct Primitive_triats<K, SEGMENT>
{
	typedef typename K::Segment_2 DataType;
};

template<typename K>
struct Primitive_triats<K, TRIANGLE>
{
	typedef typename K::Triangle_2 DataType;
};


template<typename K, typename PrimitiveTraits>
void test(const char* filename, const double duration)
{
	typedef typename PrimitiveTraits::DataType DataType;
	typedef typename std::vector<DataType> PrimitiveVector;
	typedef typename PrimitiveVector::iterator Iterator;
	typedef typename CGAL::AABB_segment_primitive<K, Iterator> Primitive;
	typedef CGAL::AABB_traits<K, Primitive> Traits;
	typedef CGAL::AABB_tree<Traits> Tree;

	//get primitive set

	PrimitiveVector primitives;

	get_segment_set<K,PrimitiveVector>(filename,primitives);
	// constructs AABB tree and internal search KD-tree with
	Tree tree(primitives.begin(),primitives.end());

	tree.accelerate_distance_queries();
	// call all tests
	test_impl<K,Tree,DataType>(tree,duration);
}

template<Primitive_type Primitive>
void test_kernels(const char *filename,
                  const double duration)
{
    std::cout << std::endl;
    std::cout << filename << std::endl;
    std::cout << "============================" << std::endl;

    std::cout << std::endl;
    std::cout << "Simple cartesian float kernel" << std::endl;
    test<CGAL::Simple_cartesian<float>,Primitive_triats<CGAL::Simple_cartesian<float>,Primitive>>(filename,duration);

    std::cout << std::endl;
    std::cout << "Cartesian float kernel" << std::endl;
    test<CGAL::Cartesian<float>,Primitive_triats<CGAL::Cartesian<float>,Primitive>>(filename,duration);

    std::cout << std::endl;
    std::cout << "Simple cartesian double kernel" << std::endl;
    test<CGAL::Simple_cartesian<double>,Primitive_triats<CGAL::Simple_cartesian<double>,Primitive>>(filename,duration);

    std::cout << std::endl;
    std::cout << "Cartesian double kernel" << std::endl;
    test<CGAL::Cartesian<double>,Primitive_triats<CGAL::Cartesian<double>,Primitive>>(filename,duration);

    std::cout << std::endl;
    std::cout << "Epic kernel" << std::endl;
    test<CGAL::Exact_predicates_inexact_constructions_kernel,Primitive_triats<CGAL::Exact_predicates_inexact_constructions_kernel,Primitive>>(filename,duration);
}


#endif /* AABB_TEST_UTIL_2D_H_ */
