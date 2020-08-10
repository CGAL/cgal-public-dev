#include <fstream>
#include <iostream>

#include "test_util.h"

enum Query_type {RAY_QUERY,
                 SEGMENT_QUERY};

template <class Tree, class K>
void test_speed_for_query(const Tree& tree,
                          const Query_type query_type,
                          const char *query_name,
                          const double duration)
{
    typedef typename K::Ray_3 Ray;
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef typename K::Segment_3 Segment;

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
                    Ray ray(source, vec);
                    tree.do_intersect(ray);
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
        }
        nb++;
    }
    unsigned int speed = (unsigned int)(nb / timer.time());
    std::cout.precision(10);
    std::cout.width(15);
    std::cout << speed << " intersections/s with " << query_name << std::endl;
    timer.stop();
}

template <class Tree, class K>
void test_speed(Tree& tree,
                const double duration)
{
    std::cout << "Test for speed" << std::endl;
    test_speed_for_query<Tree,K>(tree,RAY_QUERY,"ray",duration);
    test_speed_for_query<Tree,K>(tree,SEGMENT_QUERY,"segment",duration);
}

template<class K, class Tree, class Polyhedron>
void test_impl(Tree& tree, Polyhedron&, const double duration)
{
  test_all_intersection_query_types<Tree,K>(tree);
  test_speed<Tree,K>(tree,duration);
}

int main()
{
    const double duration = 0.1; // duration of each test
    test_kernels("./data/cube.off",duration);
    test_kernels("./data/coverrear.off",duration);
    test_kernels("./data/finger.off",duration);
    test_kernels("./data/pinion.off",duration);
    return EXIT_SUCCESS;
}
