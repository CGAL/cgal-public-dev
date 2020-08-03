#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/Embree/AABB_tree.h>

// #include <CGAL/AABB_face_graph_triangle_primitive.h>

template <class K>
int test()
{
  // types
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point;
  typedef typename K::Segment_3 Segment;

  // load polyhedron
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  Polyhedron polyhedron;
  std::ifstream ifs("./data/tetrahedron.off");
  ifs >> polyhedron;

  // construct tree from facets
  typedef CGAL::Embree::Triangle_mesh_geometry<Polyhedron, K> TriangleMesh;
  typedef CGAL::Embree::AABB_tree<TriangleMesh, K> Tree;
  typedef typename Tree::Intersection_and_primitive_id Object_and_primitive_id;

  Tree tree;
  tree.insert(polyhedron);

  // segment intersection query
  Point p((FT)-0.25,  (FT)0.251, (FT)0.255);
  Point q((FT) 0.25,  (FT)0.253, (FT)0.256);
  Segment pq(p,q);

  if(!tree.do_intersect(pq))
  {
    std::cerr << "no intersection found" << std::endl;
    return EXIT_FAILURE;
  }

  if(tree.number_of_intersected_primitives(pq) != 1)
  {
    std::cerr << "number of intersections different than one" << std::endl;
    return EXIT_FAILURE;
  }

  boost::optional<Object_and_primitive_id> any;
  any = tree.any_intersection(pq);
  if(any)
  {
    Point p = any->first;
    std::cout << "Intersection point: " << p << std::endl;
  }
  if(!any)
  {
    std::cerr << "did not find any intersection" << std::endl;
    return EXIT_FAILURE;
  }

  // // closest point query
  // Point r((FT)0.0, (FT)0.0, (FT)3.0);
  // Point closest((FT)0.0, (FT)0.0, (FT)1.0);
  // Point result = tree.closest_point(r);
  // if(result != closest)
  // {
  //   std::cerr << "wrong closest point" << std::endl;
  //   return EXIT_FAILURE;
  // }

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

/***EMACS SETTINGS***/
/* Local Variables: */
/* tab-width: 2     */
/* End:             */

