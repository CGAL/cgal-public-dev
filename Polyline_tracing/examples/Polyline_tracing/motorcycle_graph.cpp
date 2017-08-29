#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyline_tracing/Motorcycle_graph.h>

#include <vector>

int main()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

  typedef K::Point_2                                          Point_2;

  std::vector<Point_2> motorcycles;
  std::vector<Point_2> destinations;

  CGAL::Polyline_tracing::Motorcycle_graph<K> mg;
  mg.trace_motorcycle_graph(motorcycles.begin(), motorcycles.end(),
                            destinations.begin(), destinations.end());

  return 0;
}
