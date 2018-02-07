#ifndef ISLAND_TRIANGULATE_HOLE_POLYLINE_H
#define ISLAND_TRIANGULATE_HOLE_POLYLINE_H

#include <vector>
#include <tuple>
#include <stack>
#include <CGAL/Combination_enumerator.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>


namespace CGAL {
namespace internal {


// Domain structure //
// ---------------- //

template <typename PointRange>
struct Domain
{

  Domain() {}

  // boundary should always be given without the stupid last(first) one.
  // not used
  Domain(PointRange& boundary) : boundary(boundary)
  {
    std::size_t n_p = boundary.size();
    if(boundary.size() > 0)
      assert(boundary[n_p - 1] != boundary[0]);
  }

  // constructor with indices
  Domain(std::vector<int> ids) : b_ids(ids) {}

  ~Domain() {}

  void add_hole(std::vector<int>& ids)
  {
    holes_list.push_back(ids);

    for(int i=0; i<ids.size(); ++i)
      h_ids.push_back(ids[i]);
  }

  bool is_empty()
  {
    holes_list.empty() ? true : false;
  }

  std::pair<int, int> get_access_edge()
  {
    std::size_t number_of_points = b_ids.size();
    //assert(number_of_points > 0);
    assert(number_of_points >= 2);
    int i = b_ids[number_of_points - 1];
    int k = b_ids[0];

    return std::make_pair(i, k);
  }

  void print_boundary()
  {
    print(boundary);
  }

  void print_size()
  {
    std::cout << boundary.size() << std::endl;
  }

  //data
  PointRange boundary; // not used in the main algorithm
  std::vector<PointRange> holes; // not used

  std::vector<int> b_ids;
  std::vector<int> h_ids;
  std::vector<std::vector<int>> holes_list;

};


template <typename T>
void print(T &v)
{
  for(int i=0; i<v.size(); ++i)
  {
    std::cout << v[i] << " ";//<< std::endl;
  }
  std::cout << std::endl;
}

template <typename T>
void print(T &v, std::ofstream& out)
{
  out.open("data/domainV.dat", std::ofstream::app);

  for(int i=0; i<v.size(); ++i)
  {
    out << v[i] << " ";//<< std::endl;
  }
  out << std::endl;

  out.close();

}



// partition permutations // 
// ---------------------- //

struct Phi
{
  typedef std::pair<std::vector<int>, std::vector<int>> Sub_domains_pair;
  void put(std::vector<int>& left, std::vector<int>& right)
  {
    sub_domains_list.push_back(std::make_pair(left, right)); // preallocate list
  }

  std::size_t size()
  {
    return sub_domains_list.size();
  }

  bool empty()
  {
    return size() == 0 ? true : false;
  }

  std::vector<int> lsubset(const int i)
  {
    assert(i >= 0);
    assert(i < sub_domains_list.size());
    return sub_domains_list[i].first;
  }

  std::vector<int> rsubset(const int i)
  {
    assert(i >= 0);
    assert(i < sub_domains_list.size());
    return sub_domains_list[i].second;
  }

  std::vector<Sub_domains_pair> sub_domains_list;
};


void do_permutations(std::vector<std::vector<int>>& hole_list, Phi& subsets)
{
  if(hole_list.empty())
    return;

  std::vector<int> hs;

  for(int n = 0; n < hole_list.size(); ++n)
    hs.push_back(n);

  const int first = hs.front();
  const int last = hs.back();
  //std::sort(hs.begin(), hs.end()); // already sorted

  for(int s = 0; s <= hs.size(); ++s) // s = number of holes on one (left) side
  {
    std::vector<int> p1(s);
    std::vector<int> p2(hole_list.size() - s);

    if(s == 0)
    {
      subsets.put(p1, hs);

      print(p1); std::cout << "-- "; print(hs); std::cout << std::endl;
      continue;
    }

    CGAL::Combination_enumerator<int> permutations(s, first, last+1);

    int p = 0;
    while(!permutations.finished())
    {
      for(int i=0; i<s; ++i)
      {
        p1[i] = permutations[i];
      }

      ++permutations;

      std::sort(p1.begin(), p1.end());
      std::set_symmetric_difference(p1.begin(), p1.end(), hs.begin(), hs.end(),
                                    p2.begin());

      assert(p1.size() == s);
      assert(p2.size() == hs.size() - s);

      print(p1); std::cout << "-- "; print(p2); std::cout << std::endl;

      subsets.put(p1, p2);
    }

  }
}

// split //
// ----- //

template <typename PointRange>
void split_domain(const Domain<PointRange>& init_domain,
                    Domain<PointRange>& left_dom, Domain<PointRange>& right_dom,
                    const int& i, const int& pid, const int& k)
{

  typedef std::vector<int> Ids;
  Ids ids = init_domain.b_ids;
  Ids left;
  Ids right;

  // i, k indices of access edge = first and last

  // find position of pid
  Ids::iterator it;
  it = find(ids.begin(), ids.end(), pid);
  assert(it != ids.end());


  left.insert(left.begin(), ids.begin(), it + 1);
  right.insert(right.begin(), it, ids.end());


  //assert(left.front() == i);
  assert(left.back() == pid);
  assert(right.front() == pid);
  //assert(right.back() == k); // maybe switch i, and k

  left_dom.b_ids = left;
  right_dom.b_ids = right;

}


void reorder_island(std::vector<int>& h_ids, const int& v)
{

  std::vector<int>::iterator it = find(h_ids.begin(), h_ids.end(), v);
  assert(it != h_ids.end());

  // 2) rotate by the third vertex of t
  std::size_t dist = std::distance(h_ids.begin(), it); // std::size_t?
  std::rotate(h_ids.begin(), h_ids.begin() + dist, h_ids.end());

  // 3) add the first removed element
  h_ids.push_back(h_ids[0]);

  // 4) reverse. Todo: Check and do it iff reversal is needed.
  std::reverse(h_ids.begin(), h_ids.end());

}


void merge_id_sets(std::vector<int>& b_ids,
                   const int& i, const int& v, const int& k,
                   std::vector<int>& hole_ids)
{

  // temp solution: add last vertex so that it's closed.
  //b_ids.push_back(b_ids[0]);

  // mst take cases with & without reordering
  reorder_island(hole_ids, v);

  std::size_t initial_b_size = b_ids.size();

  // insertion position = just before k
  //typename std::vector<int>::iterator insertion_point = b_ids.begin() + k;

  // insertion position = just after i
  // i is at position n - 1 = last element. Can I just append at the end always?
  //typename std::vector<int>::iterator insertion_point = b_ids.begin() + (initial_b_size - 1);

  // just append at the end - assuming that i is the last point on b_ids
  // and k is the first. t triangle is (i, v, k)
  typename std::vector<int>::iterator insertion_point = b_ids.end();

  b_ids.insert(insertion_point, hole_ids.begin(), hole_ids.end());

  assert(*(b_ids.begin() + i) == b_ids[i] );
  assert(b_ids.size() == initial_b_size + hole_ids.size());

}

template<typename PointRange>
void join_domain(const Domain<PointRange>& domain, Domain<PointRange>& new_domain,
                 const int& i, const int& v, const int& k)
{
  typedef std::vector<int> Ids;
  Ids id_set = domain.b_ids;
  Ids hole_ids = domain.h_ids; // for now assume just one hole.

  merge_id_sets(id_set, i, v, k, hole_ids);
  new_domain.b_ids = id_set;
}


struct Tracer
{

  template <class LookupTable>
  void operator()(const LookupTable& lambda, int v0, int v1)
  {
    CGAL_assertion_code( const int n = lambda.n; )
    std::stack<std::pair<int, int> > ranges;
    ranges.push(std::make_pair(v0, v1));

    while(!ranges.empty()) {
      std::pair<std::size_t, std::size_t> r = ranges.top();
      ranges.pop();
      CGAL_assertion(r.first >= 0 && r.first < n);
      CGAL_assertion(r.second >= 0 && r.second < n);

      // if on border
      if(r.first + 1 == r.second) {
        continue; }
      if(r.first + 1 == n) {
        continue;
      }

      std::size_t la = lambda.get(r.first, r.second);
      if(la == -1) {
          std::cerr << "out hole" << std::endl;
          //*out_hole++ = std::make_pair(r.first, r.second);
        continue;
      }

      CGAL_assertion(la >= 0 && la < n);
     // CGAL_assertion(r.first < la && r.second > la); not with islands
      //auto triangle = std::make_tuple(r.first, la, r.second);

      std::vector<std::size_t> triangle = {r.first, la, r.second};
      collection.push_back(triangle);

      ranges.push(std::make_pair(r.first, la));
      ranges.push(std::make_pair(la, r.second));
    }
  }


  //std::vector<std::tuple<int, int, int>> collection;
  std::vector<std::vector<std::size_t>> collection;
};






template<typename PointRange, typename WeightCalculator,
         typename WeightTable, typename LambdaTable>
class Triangulate
{
  typedef typename WeightCalculator::Weight Weight;


public:

  Triangulate(Domain<PointRange> domain,
              PointRange allpoints,
              WeightTable& W,
              LambdaTable& l,
              const WeightCalculator & WC) :
              points(allpoints),
              W(W),
              lambda(l),
              domain(domain),
              WC(WC){}

  std::size_t do_triangulation(int& i, int& k, std::size_t& count)
  {

    init_triangulation();

    processDomain(domain, i, k, count);

    //ambda.print("data/lambda-rec.dat");
    //W.print("data/weight-rec.dat");
  }

  void collect_triangles(std::vector<std::vector<std::size_t>>& triplets,
                         int& i, int& k)
  {
    Tracer tracer;
    tracer(lambda, i, k);
    triplets = tracer.collection;
  }

  template <typename PolygonMesh>
  void visualize(PointRange& points, std::vector<std::vector<std::size_t>>& polygon_soup,
                 PolygonMesh& mesh)
  {
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygon_soup, mesh);
  }




private:


  void init_triangulation()
  {
    for(int i=0; i < domain.b_ids.size(); ++i)
    {
      init_b.push_back(domain.b_ids[i]);
    }

    for(int i=0; i < domain.h_ids.size(); ++i)
    {
      init_island.push_back(domain.h_ids[i]);
    }

  }


  // main loop //
  // --------- //
  void processDomain(Domain<PointRange> domain, int& i, int& k, std::size_t& count)
  {
    // (i, k) = acccess edge


    std::cout << "count: " << count << std::endl;

    //print(domain.b_ids, out_domain);



    // domains consisting of only one edge
    if(domain.b_ids.size() == 2)
      return;

    // base case
    if(domain.b_ids.size() == 3 && domain.holes_list.empty())
    {
      //assert(domain.b_ids[0] == i); // access edge source
      //assert(domain.b_ids[2] == k); // access edge target


      int m = domain.b_ids[1]; //third vertex



      ///////////////////////////////////////////////////////////////
      //std::cout<<"Evaluating t= ("<<i<<","<<m<<","<<k<<")"<<std::endl;
      calculate_weight(i, m, k);
      count++;

      return;
    }
    assert(domain.b_ids.size() >= 3);

    // pid : third vertex

    // CASE I - if there are islands, join until there are no islands.
    for(int pid : domain.h_ids)
    {
      //std::cout << "i= " << i << " k= " << k << std::endl;
      //std::cout << "pid= " << pid << std::endl;

      // avoid source & target of e_D
      if(pid == i || pid == k)
      {
        //std::cout << " point aborted" << std::endl;
        continue;
      }

      Domain<PointRange> D1;
      join_domain(domain, D1, i, pid, k);

      // get a new e_D
      std::pair<int, int> e_D1 = D1.get_access_edge();

      processDomain(D1, e_D1.first, e_D1.second, count);


      //////////////////////////////////////////////////////////////
      // calculate weight of triangle t - after the subdomains left and right have been checked
      int m = pid; //third vertex


      //std::cout<<"Evaluating t= ("<<i<<","<<m<<","<<k<<")"<<std::endl;
      calculate_weight(i, m, k);
      count++;

    }

    // CASE II
    for(int pid : domain.b_ids)
    {
      //std::cout << "i= " << i << " k= " << k << std::endl;
      //std::cout << "pid= " << pid << std::endl;

      // avoid source & target of e_D
      if(pid == i || pid == k)
      {
        //std::cout << " point aborted" << std::endl;
        continue;
      }

      // split to two sub-domains
      Domain<PointRange> D1;
      Domain<PointRange> D2;
      // essentially splitting boundaries
      split_domain(domain, D1, D2, i, pid, k);
      // D1, D2 have just new boundaries - no hole information.

      // get new access edges for each
      std::pair<int, int> e_D1 = D1.get_access_edge();
      std::pair<int, int> e_D2 = D2.get_access_edge();

      // assign all combination of holes to subdomains and process each pair
      Phi partition_space;
      do_permutations(domain.holes_list, partition_space);

      if(partition_space.empty())
      {
        // when the domain has been merged so that there is no holes inside
        processDomain(D1, e_D1.first, e_D1.second, count);
        processDomain(D2, e_D2.first, e_D2.second, count);
      }
      else
      {
        // when t is formed with a vertex on the boundary of a domain with holes
        for(std::size_t p = 0; p < partition_space.size(); ++p)
        {
          std::vector<int> lholes = partition_space.lsubset(p);
          std::vector<int> rholes = partition_space.rsubset(p);

          for(int lh : lholes)
            D1.add_hole(domain.holes_list[lh]);

          for(int rh : rholes)
            D2.add_hole(domain.holes_list[rh]);

          processDomain(D1, e_D1.first, e_D1.second, count);
          processDomain(D2, e_D2.first, e_D2.second, count);
        }

      }

      ///////////////////////////////////////////////////////////////
      // calculate weight of triangle t - after the subdomains left and right have been checked
      int m = pid; //third vertex
      //std::cout<<"Evaluating t= ("<<i<<","<<m<<","<<k<<")"<<std::endl;


      calculate_weight(i, m, k);
      count++;

    }
  }


  bool are_vertices_in_island(const int& i, const int& m, const int& k)
  {

    std::vector<int>::iterator it1, it2, it3;

    it1 = std::find(init_island.begin(), init_island.end(), i);
    it2 = std::find(init_island.begin(), init_island.end(), m);
    it3 = std::find(init_island.begin(), init_island.end(), k);

    return (it1 != init_island.end()) && (it2 != init_island.end()) && (it3 != init_island.end()) ?  true : false;

  }


  void calculate_weight(int& i, int& m, int& k)
  {

    if(are_vertices_in_island(i, m, k))
      return;


    // i, m, k are global indices
    assert(m != i);
    assert(m != k);

    if(i > k)
    {
      std::swap(i, k); // needed to store the edge (i,k) sorted. Maybe move this in the Look_up_map.
    }
    assert(i < k);


    PointRange Q;
    const Weight& w_imk = WC(points, Q, i, m, k, lambda);


    if(w_imk == Weight::NOT_VALID())
    {
      std::cerr << "non-manifold edge"  << std::endl;
      return;
    }

    const Weight& w = W.get(i,m) + W.get(m,k) + w_imk;

    if(lambda.get(i, k) == -1 || w < W.get(i, k)) {
      W.put(i, k, w);
      lambda.put(i, k, m);

      //W.print("data/weight.dat");
      //lambda.print("data/lambda.dat");
      std::cout << std::endl;
    }
  }



  // data
  PointRange points;

  //std::set<std::vector<std::size_t>> memoized;

  WeightTable W;
  LambdaTable lambda;

  Domain<PointRange> domain;
  const WeightCalculator& WC;


  // initial boundary & island
  std::vector<int> init_b;
  std::vector<int> init_island;


  std::ofstream out_domain;



};





} // namespace internal
} // namsepace CGAL





#endif // ISLAND_TRIANGULATE_HOLE_POLYLINE_H
