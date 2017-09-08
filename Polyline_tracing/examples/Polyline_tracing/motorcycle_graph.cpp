#define CGAL_MOTORCYCLE_GRAPH_VERBOSE

//#define CGAL_MOTORCYCLE_GRAPH_USE_MANUAL_POINTS
#define CGAL_MOTORCYCLE_GRAPH_RANDOM_POINTS_IN_SQUARE
//#define CGAL_MOTORCYCLE_GRAPH_RANDOM_POINTS_ON_SEGMENT

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

  bool is_loop_infinite = true;
  while(is_loop_infinite)
  {
//    is_loop_infinite = false;

    std::vector<Point_2> motorcycles;
    std::vector<Point_2> destinations;

#ifdef CGAL_MOTORCYCLE_GRAPH_USE_MANUAL_POINTS
    // add some motorcycles
    motorcycles.push_back(Point_2(0,0));
    motorcycles.push_back(Point_2(0,0));
    motorcycles.push_back(Point_2(0,0));

    // and their respective destinations
    destinations.push_back(Point_2(10,0));
    destinations.push_back(Point_2(5,0));
    destinations.push_back(Point_2(5,5));
#else // random stuff below
    const int size = 50; // number of random points

 #ifdef CGAL_MOTORCYCLE_GRAPH_RANDOM_POINTS_IN_SQUARE
    typedef CGAL::Random_points_in_square_2<Point_2>            Generator;
    const int side = 10;
  #if CGAL_MOTORCYCLE_GRAPH_USE_FIXED_SEEDS
    CGAL::Random rand_s(6), rand_d(7);
    std::cerr << "Seeds = " << rand_s.get_seed() << " " << rand_d.get_seed() << std::endl;
    Generator gen_s(side, rand_s), gen_d(side, rand_d);
  #else
    Generator gen_s(side), gen_d(side);
  #endif

    CGAL::copy_n_unique(gen_s, size, std::back_inserter(motorcycles));
    CGAL::copy_n_unique(gen_d, size, std::back_inserter(destinations));

 #elif defined(CGAL_MOTORCYCLE_GRAPH_RANDOM_POINTS_ON_SEGMENT)
    typedef CGAL::Random_points_on_segment_2<Point_2>            Generator;
    const Point_2 s0(0,0), t0(4,0), s1(CGAL_PI, -0.1), t1(-CGAL::sqrt(3.), std::cos(0.1));
    Generator gen_s0(s0, t0), gen_d0(s0, t0);
    Generator gen_s1(s1, t1), gen_d1(s1, t1);

    CGAL::copy_n_unique(gen_s0, size, std::back_inserter(motorcycles));
    CGAL::copy_n_unique(gen_d0, size, std::back_inserter(destinations));
    CGAL::copy_n_unique(gen_s1, size, std::back_inserter(motorcycles));
    CGAL::copy_n_unique(gen_d1, size, std::back_inserter(destinations));
 #endif

#endif // CGAL_MOTORCYCLE_GRAPH_USE_MANUAL_POINTS

    // trace the graph
    CGAL::Polyline_tracing::Motorcycle_graph<K> motorcycle_graph;
    motorcycle_graph.trace_motorcycle_graph(motorcycles.begin(), motorcycles.end(),
                                            destinations.begin(), destinations.end());
    CGAL_postcondition(motorcycle_graph.is_valid());

    // output
    motorcycle_graph.output_motorcycles_sources_and_destinations();
    motorcycle_graph.output_all_dictionary_points();
    for(std::size_t i=0; i<motorcycles.size(); ++i)
    {
  //    motorcycle_graph.motorcycle(i).output_intended_track();
      motorcycle_graph.motorcycle(i).output_track();
    }
  }

  return 0;
}
