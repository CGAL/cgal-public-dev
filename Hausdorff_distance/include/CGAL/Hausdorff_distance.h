#ifndef CGAL_HAUSDORFF_DISTANCE_H
#define CGAL_HAUSDORFF_DISTANCE_H

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Timer.h>
#include <CGAL/hilbert_sort.h>

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
    tbb::parallel_reduce( tbb::blocked_range<double*>( &(dist[0]), &(dist[A.size()])), m);
    res = (std::max)(res, m.value);
  }
#endif
  tq.stop();

  std::cerr << "build time: " << tb.time() << " sec." << std::endl; 
  std::cerr << "query time: " << tq.time() << " sec." << std::endl; 
  return sqrt(res);
}

  /*
p = first point of P
while(not at end of P)
  q = find closest point(p, T)
  D = dist(p,q);
  do {
    p = next point of P
    pair<point,bool> pa = circular_range_query(p,D,T)
    if (range empty) {
      break;  //the outer loop will compute a larger D
    }
    // D is still the maximum, so continue inner loop
  }

   */


template <typename PointRange>
double Hausdorff_distance_rq_2(const PointRange& A, const PointRange& B)
{
  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_2 Point_2;
  typedef CGAL::Search_traits_2<K> TreeTraits;
  typedef CGAL::Fuzzy_sphere<TreeTraits> Fuzzy_circle;
  typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;

  Tree treeB(B.begin(), B.end());
  treeB.build();

 Tree treeA(A.begin(), A.end());
  treeA.build();

  double res = 0;
  typename PointRange::const_iterator b = A.begin(), e = A.end();
  Point_2 a = *b; 
  while(b != e){
    Neighbor_search search(treeB, a, 1);
    Point_2 q = search.begin()->first;
    double sd = search.begin()->second;
    res = (std::max)(res,sd);
    while(true){
      ++b;
      if(b != e){
        a = *b;
        Fuzzy_circle circ(a,sqrt(res));
        std::list<Point_2> result;
        treeB.search(std::back_inserter( result ),circ);
        if(result.empty()){
          break;
        }
      } else {
        break;
      }
    }
  }

{
  typename PointRange::const_iterator b = B.begin(), e = B.end();
  Point_2 a = *b; 
  while(b != e){
    Neighbor_search search(treeA, a, 1);
    Point_2 q = search.begin()->first;
    double sd = search.begin()->second;
    res = (std::max)(res,sd);
    while(true){
      ++b;
      if(b != e){
        a = *b;
        Fuzzy_circle circ(a,sqrt(res));
        std::list<Point_2> result;
        treeB.search(std::back_inserter( result ),circ);
        if(result.empty()){
          break;
        }
      } else {
        break;
      }
    }
  }
}
  return sqrt(res);
 }



template <typename PointRange>
double Hausdorff_distance_hs_2(const PointRange& A, const PointRange& B)
{
  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_2 Point_2;
  typedef CGAL::Search_traits_2<K> TreeTraits;
  typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;
  std::vector<Point_2> AQ(A.begin(),A.end());
  
  hilbert_sort(AQ.begin(),AQ.end());

  Tree treeB(B.begin(), B.end());
  treeB.build();
  
  double res = 0;
  Point_2 q(0,0);

  for(int i=0; i < AQ.size(); i++){
    const Point_2& a = AQ[i];
    if((i==0) || (squared_distance(a,q) > res)){
      Neighbor_search search(treeB, a, 1);
      q = search.begin()->first;
      double sd = search.begin()->second;
      res = (std::max)(sd, res);
    }
  }

  std::vector<Point_2> BQ(B.begin(),B.end());
  
  hilbert_sort(BQ.begin(),BQ.end());

  Tree treeA(A.begin(), A.end());
  treeB.build();
  

  for(int i=0; i < BQ.size(); i++){
    const Point_2& a = BQ[i];
    if((i==0) || (squared_distance(a,q) > res)){
      Neighbor_search search(treeA, a, 1);
      q = search.begin()->first;
      double sd = search.begin()->second;
      res = (std::max)(sd, res);
    }
  }
  return sqrt(res);
 }
} // namespace CGAL


#endif // CGAL_HAUSDORFF_DISTANCE_H
