#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>


#include <CGAL/Embree/AABB_tree.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Timer.h>

#include <CGAL/disable_warnings.h>

// #include <CGAL/internal/AABB_tree/Primitive_helper.h>
#include <CGAL/use.h>

#include <boost/mem_fn.hpp>

double random_in(const double a, const double b)
{
    double r = rand() / (double)RAND_MAX;
    return a + (b - a) * r;
}

template <class K>
typename K::Point_3 random_point_in(const CGAL::Bbox_3& bbox)
{
    typedef typename K::FT FT;
    FT x = (FT)random_in(bbox.xmin(),bbox.xmax());
    FT y = (FT)random_in(bbox.ymin(),bbox.ymax());
    FT z = (FT)random_in(bbox.zmin(),bbox.zmax());
    return typename K::Point_3(x,y,z);
}

template <class K>
typename K::Vector_3 random_vector()
{
    typedef typename K::FT FT;
    FT x = (FT)random_in(0.0,1.0);
    FT y = (FT)random_in(0.0,1.0);
    FT z = (FT)random_in(0.0,1.0);
    return typename K::Vector_3(x,y,z);
}


template <class Tree, class K>
void test_all_intersection_query_types(Tree& tree)
{
    std::cout << "Test all query types" << std::endl;

    typedef typename K::FT FT;
    typedef typename K::Ray_3 Ray;
    typedef typename K::Line_3 Line;
    typedef typename K::Point_3 Point;
    typedef typename K::Segment_3 Segment;
    typedef typename Tree::Primitive Primitive;

    Point p((FT)-0.5, (FT)-0.5, (FT)-0.5);
    Point q((FT) 0.5, (FT) 0.5, (FT) 0.5);
    Ray ray(p,q);
    Segment segment(p,q);
    bool success = false;

    // do_intersect
    success = tree.do_intersect(ray);
    success = tree.do_intersect(segment);
    (void) success;

    // number_of_intersected_primitives
    tree.number_of_intersected_primitives(ray);
    tree.number_of_intersected_primitives(segment);

    // all_intersected_primitives
    std::list<typename Tree::Primitive_id> primitives;
    tree.all_intersected_primitives(ray,std::back_inserter(primitives));
    tree.all_intersected_primitives(segment,std::back_inserter(primitives));

    // any_intersection
    boost::optional<typename Tree::Intersection_and_primitive_id> r = tree.any_intersection(ray);
    boost::optional<typename Tree::Intersection_and_primitive_id> s = tree.any_intersection(segment);

    // any_intersected_primitive
    boost::optional<typename Tree::Primitive_id> optional_primitive;
    optional_primitive = tree.any_intersected_primitive(ray);
    optional_primitive = tree.any_intersected_primitive(segment);

    // all_intersections
    std::list< boost::optional<typename Tree::Intersection_and_primitive_id> > intersections_r;
    std::list< boost::optional<typename Tree::Intersection_and_primitive_id> > intersections_s;
    tree.all_intersections(ray,std::back_inserter(intersections_r));
    tree.all_intersections(segment,std::back_inserter(intersections_s));
}


template <class Tree, class K>
void test_all_distance_query_types(Tree& tree)
{
    typedef typename K::FT FT;
    typedef typename K::Point_3 Point;
    typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;

    Point query = random_point_in<K>(tree.bbox());

    FT sqd1 = tree.squared_distance(query);

    Point p1 = tree.closest_point(query);

    Point_and_primitive_id pp1 = tree.closest_point_and_primitive(query);
}


template <class Tree, class K>
void test_distance_speed(Tree& tree,
                         const double duration)
{
    typedef typename K::Point_3 Point;

    CGAL::Timer timer;
    timer.start();
    unsigned int nb = 0;
    while(timer.time() < duration)
    {
            // picks a random point in the tree bbox
            Point query = random_point_in<K>(tree.bbox());
            Point closest = tree.closest_point(query);
      (void) closest;
            nb++;
    }
    double speed = (double)nb / timer.time();
    std::cout << speed << " distance queries/s" << std::endl;
    timer.stop();
}


/**
 * Declaration only, implementation should be given in .cpp file
 */
template<class K, class Tree, class Polyhedron>
void test_impl(Tree& tree, Polyhedron& p, const double duration);


/**
 * Generic test method. Build AABB_tree and call test_impl()
 */
template <class K>
void test(const char *filename,
          const double duration)
{
    typedef CGAL::Polyhedron_3<K> Polyhedron;
    typedef CGAL::Embree::Triangle_mesh_geometry<Polyhedron, K> TriangleMesh;
    typedef CGAL::Embree::AABB_tree<TriangleMesh, K> Tree;

    Polyhedron polyhedron;
    std::ifstream ifs(filename);
    ifs >> polyhedron;

    Tree tree;
    tree.insert(polyhedron);

    // call all tests
    test_impl<K,Tree,Polyhedron>(tree,polyhedron,duration);
}


/**
 * Generic test_kernel method. call test<K> for various kernel K.
 */
void test_kernels(const char *filename,
                  const double duration)
{
    std::cout << std::endl;
    std::cout << "Polyhedron " << filename << std::endl;
    std::cout << "============================" << std::endl;

    std::cout << std::endl;
    std::cout << "Simple cartesian float kernel" << std::endl;
    test<CGAL::Simple_cartesian<float>>(filename,duration);

    std::cout << std::endl;
    std::cout << "Cartesian float kernel" << std::endl;
    test<CGAL::Cartesian<float>>(filename,duration);

    std::cout << std::endl;
    std::cout << "Simple cartesian double kernel" << std::endl;
    test<CGAL::Simple_cartesian<double>>(filename,duration);

    std::cout << std::endl;
    std::cout << "Cartesian double kernel" << std::endl;
    test<CGAL::Cartesian<double>>(filename,duration);

    std::cout << std::endl;
    std::cout << "Epic kernel" << std::endl;
    test<CGAL::Exact_predicates_inexact_constructions_kernel>(filename,duration);
}

#include <CGAL/enable_warnings.h>
