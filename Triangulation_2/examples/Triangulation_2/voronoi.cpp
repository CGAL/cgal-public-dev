#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <boost/foreach.hpp>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_2<K>  Triangulation;
typedef Triangulation::Edge                Edge;
typedef Triangulation::Point               Point;

int main( )
{
  std::ifstream in("data/voronoi.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;
  Triangulation T;
  T.insert(begin, end);

  int ns = 0;
  int nr = 0;
  BOOST_FOREACH(Edge &e, T.edges()) {
    CGAL::Object o = T.dual(e);
    if (CGAL::object_cast<K::Segment_2>(&o)) {++ns;}
    else if (CGAL::object_cast<K::Ray_2>(&o)) {++nr;}
  }
  std::cout << "The Voronoi diagram has " << ns << " finite edges "
	    << " and " << nr << " rays" << std::endl;
  return 0;
}
