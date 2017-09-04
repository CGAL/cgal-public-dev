#define CGAL_MOTORCYCLE_GRAPH_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyline_tracing/Motorcycle_graph.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h> // for 'copy_n_unique()'

#include <vector>

int main()
{
  std::cout.precision(17);
  std::cerr.precision(17);

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

  typedef K::Point_2                                          Point_2;

  std::vector<Point_2> motorcycles;
  std::vector<Point_2> destinations;

#if 0
  // add some motorcycles
  motorcycles.push_back(Point_2(-1,0));
  motorcycles.push_back(Point_2(1,-5));
  motorcycles.push_back(Point_2(3,2));
  motorcycles.push_back(Point_2(2,-3));

  // and their respective destinations
  destinations.push_back(Point_2(3,0));
  destinations.push_back(Point_2(1,1));
  destinations.push_back(Point_2(-1,-1));
  destinations.push_back(Point_2(2,1));
#else
  typedef CGAL::Random_points_in_square_2<Point_2>            Generator;

  const int side = 10;
  const int size = 50;

#if 0 // if using specific seeds
  CGAL::Random rand_s(6), rand_d(7);
  std::cerr << "Seeds = " << rand_s.get_seed() << " " << rand_d.get_seed() << std::endl;
  Generator gen_s(side, rand_s), gen_d(side, rand_d);
#else
  Generator gen_s(side), gen_d(side);
#endif

  CGAL::copy_n_unique(gen_s, size, std::back_inserter(motorcycles));
  CGAL::copy_n_unique(gen_d, size, std::back_inserter(destinations));
#endif

  // trace the graph
  CGAL::Polyline_tracing::Motorcycle_graph<K> motorcycle_graph;
  motorcycle_graph.trace_motorcycle_graph(motorcycles.begin(), motorcycles.end(),
                                          destinations.begin(), destinations.end());

  // output
  motorcycle_graph.output_motorcycles_sources_and_destinations();
  motorcycle_graph.output_all_dictionary_points();
  for(std::size_t i=0; i<motorcycles.size(); ++i)
  {
    motorcycle_graph.motorcycle(i).output_intended_track();
    motorcycle_graph.motorcycle(i).output_track();
  }

  return 0;
}
