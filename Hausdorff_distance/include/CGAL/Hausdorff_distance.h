#ifndef CGAL_HAUSDORFF_DISTANCE_H
#define CGAL_HAUSDORFF_DISTANCE_H

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Timer.h>
#include <boost/foreach.hpp>

#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>


namespace CGAL {

  namespace HD {


class Apply {

  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_2 Point_2;
  typedef CGAL::Search_traits_2<K> TreeTraits;
  typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;

  Tree& tree;
  Point_2 const * points;
  double* dist;

public:
  void operator()(const tbb::blocked_range<std::size_t>& r) const {
      Point_2 const * a = points; 
      for( std::size_t i = r.begin(); i != r.end(); ++i){
        Neighbor_search search(tree, a[i], 1);
        dist[i] = search.begin()->second;
      }
  }

  Apply(Tree& tree, const std::vector<Point_2>& points, std::vector<double>& dist)
    : tree(tree), points(&(points[0])), dist(&(dist[0]))
  {}
};

struct Max {
    double value;
    Max() : value(0) {}
  Max( Max& s, tbb::split ) {value = 0;}
  void operator()( const tbb::blocked_range<double*>& r ) {
        double temp = value;
        for( double* a=r.begin(); a!=r.end(); ++a ) {
          temp = (std::max)(temp,*a);
        }
        value = temp;
    }
  void join( Max& rhs ) {value = (std::max)(value,rhs.value);}
};

  }

template <typename PointRange>
double Hausdorff_distance_2(const PointRange& A, const PointRange& B)
{
  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_2 Point_2;
  typedef CGAL::Search_traits_2<K> TreeTraits;
  typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;
  
  std::vector<double> dist((std::max)(A.size(),B.size()));

  Timer tb, tq;
  tb.start();
  Tree treeB(B.begin(), B.end());
  treeB.build();
  tb.stop();

  tq.start();
  double res = 0;
#if 0  
  BOOST_FOREACH(const Point_2& a, A){
    Neighbor_search search(treeB, a, 1);
    double sd = search.begin()->second;
    res = (std::max)(sd, res);
  } 
#else
  {
    tbb::parallel_for(tbb::blocked_range<size_t>(0,A.size()),HD::Apply(treeB, A, dist));
    HD::Max m;
    tbb::parallel_reduce( tbb::blocked_range<double*>( &(dist[0]), &(dist[A.size()])), m);
    res = m.value;
  }
#endif

  tq.stop();
  
  tb.start();
  Tree treeA(A.begin(), A.end());
  treeA.build();
  tb.stop();

  tq.start();
#if 0  
  BOOST_FOREACH(const Point_2& b, B){
    Neighbor_search search(treeA, b, 1);
    double sd = search.begin()->second;
    res = (std::max)(sd, res);
  }
#else
  {
    tbb::parallel_for(tbb::blocked_range<size_t>(0,B.size()),HD::Apply(treeA, B, dist));
    HD::Max m;
    tbb::parallel_reduce( tbb::blocked_range<double*>( &(dist[0]), &(dist[B.size()])), m);
    res = (std::max)(res, m.value);
  }
#endif
  tq.stop();

  std::cerr << "build time: " << tb.time() << " sec." << std::endl; 
  std::cerr << "query time: " << tq.time() << " sec." << std::endl; 
  return sqrt(res);
}
  
} // namespace CGAL


#endif // CGAL_HAUSDORFF_DISTANCE_H
