#ifndef ISLAND_TRIANGULATE_HOLE_POLYLINE_H
#define ISLAND_TRIANGULATE_HOLE_POLYLINE_H

#include <vector>
#include <CGAL/Combination_enumerator.h>

namespace CGAL {
namespace internal {


// Domain structure //
// ---------------- //

template <typename PointRange>
struct Domain
{

  Domain() {}

  // boundary should always be given without the stupid last(first) one.
  Domain(PointRange& boundary) : boundary(boundary)
  {
    std::size_t n_p = boundary.size();
    if(boundary.size() > 0)
      assert(boundary[n_p - 1] != boundary[0]);
  }

  // constructor with indices
  Domain(std::vector<int> ids) : b_ids(ids) {}

  ~Domain() {}

  void add_hole(PointRange& hole)
  {
    holes.push_back(hole);

    // todo: deal with hanging point when merging multiple holes
    for(int i=0; i<hole.size(); ++i)
    {
      // for now holeVertices is a concatenation of all hole points
      holeVertices.push_back(hole[i]);
    }
  }

  void add_hole(std::vector<int>& ids)
  {
    holes_list.push_back(ids);

    for(int i=0; i<ids.size(); ++i)
      h_ids.push_back(ids[i]);
  }

  bool is_empty()
  {
    holes.empty() ? true : false;
  }

  std::pair<int, int> get_access_edge()
  {
    std::size_t number_of_points = boundary.size();
    assert(number_of_points > 0);
    int i = number_of_points - 1; // the one before the last one.
    int k = 0;


    return std::make_pair(i, k);
  }

  std::pair<int, int> get_access_edge_ids()
  {
    std::size_t number_of_points = b_ids.size();
    assert(number_of_points > 0);
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
  PointRange boundary;
  std::vector<PointRange> holes;
  PointRange holeVertices;

  std::vector<int> b_ids;
  std::vector<int> h_ids;
  std::vector<std::vector<int>> holes_list;

};


template <typename PointRange>
void print(PointRange &v)
{
  for(int i=0; i<v.size(); ++i)
  {
    std::cout << v[i] << " ";//<< std::endl;
  }
  //std::cout << std::endl;
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


template<typename PointRange>
void do_permutations(std::vector<PointRange>& holes, Phi& subsets)
{
  if(holes.empty())
    return;

  // enumerate holes
  std::vector<int> hs;
  for(int n = 0; n < holes.size(); ++n)
    hs.push_back(n);

  const int first = hs.front();
  const int last = hs.back();
  //std::sort(hs.begin(), hs.end()); // already sorted

  for(int s = 0; s <= hs.size(); ++s) // s = number of holes on one (left) side
  {
    std::vector<int> p1(s);
    std::vector<int> p2(holes.size() - s);

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
void split_domain(Domain<PointRange>& init_domain, Domain<PointRange>& left_dom, Domain<PointRange>& right_dom,
                  const int& i, const int& v, const int& k)
{
  PointRange b_points = init_domain.boundary;
  PointRange left;// = left_dom.boundary;
  PointRange right;// = right_dom.boundary;

  // i, k indices of access edge
  // v index of third vertex forming triangle t.

  int next_v = v;

  // left subset
  left.push_back(b_points[next_v]);
  while (next_v != i) {

    if(next_v == b_points.size() - 1) // if we reached the last element
      next_v = 0;
    else
      next_v++;

    left.push_back(b_points[next_v]);
  }

  // add the new last(=first)
  //left.push_back(b_points[v]);

  // right subset
  next_v = k;
  right.push_back(b_points[next_v]);
  while (next_v != v) {

    if(next_v == b_points.size() - 1)
      next_v = 0;
    else
      next_v++;

    right.push_back(b_points[next_v]);
  }

  // add the new last(=first)
  //right.push_back(b_points[k]);

  assert(left.front() == b_points[v]);
  assert(left.back() == b_points[i]);
  assert(right.front() == b_points[k]);
  assert(right.back() == b_points[v]);

  left_dom.boundary = left;
  right_dom.boundary = right;

}

template <typename PointRange>
void split_domain_i(const Domain<PointRange>& init_domain,
                    Domain<PointRange>& left_dom, Domain<PointRange>& right_dom,
                    const int& i, const int& pid, const int& k)
{
  typedef std::vector<int> Ids;
  Ids ids = init_domain.b_ids;
  Ids left;
  Ids right;
  const int n = ids.size();

  // i, k indices of access edge

  Ids::iterator it;
  it = find(ids.begin(), ids.end(), pid);
  assert(it != ids.end());

  // left subset
  left.push_back(*it);

  while (*it != i) {

    if(it == ids.end()-1)
      it = ids.begin();
    else
      ++it;

    left.push_back(*it);
  }

  // right subset

  it = find(ids.begin(), ids.end(), k);
  assert(it != ids.end());

  right.push_back(*it);

  while (*it != pid) {

    if(it == ids.end()-1)
      it = ids.begin();
    else
      ++it;

    right.push_back(*it);
  }


  assert(left.front() == pid);
  assert(left.back() == i);
  assert(right.front() == k);
  assert(right.back() == pid);

  left_dom.b_ids = left;
  right_dom.b_ids = right;

}


// join //
// ---- //
/*
template<typename PointRange>
void join_domain(const Domain<PointRange>& domain, Domain<PointRange>& new_domain,
                 const int& i, const int& v, const int& k)
{
  PointRange boundary_point_set = domain.boundary;
  PointRange hole = domain.holeVertices; // for now assume just one island. todo for more.

  merge_point_sets(boundary_point_set, i, v, k, hole);

  new_domain.boundary = boundary_point_set;

}

template <typename PointRange>
void merge_point_sets(PointRange& boundary,
                      const int& i, const int& v, const int& k,
                      PointRange& hole)
{
  // insert hole on the boundary list

  reorder_island(hole, v);

  std::size_t initial_b_size = boundary.size();

  // insertion point = just before k
  typename PointRange::iterator insertion_point = boundary.begin() + k;

  boundary.insert(insertion_point, hole.begin(), hole.end());

  assert(*(boundary.begin() + i) == boundary[i] );
  assert(boundary.size() == initial_b_size + hole.size());

}

template <typename PointRange>
void reorder_island(PointRange& hole, const int& v)
{
  assert(v >= 0);
  assert(v < ho.size());

  // 1) take the last(=first) out
  //hole.pop_back();

  // 2) rotate by the third vertex of t
  std::rotate(hole.begin(), hole.begin() + v, hole.end());

  // 3) add the first removed element
  hole.push_back(hole[0]);

  // 4) reverse. Todo: Check and do it iff reversal is needed.
  std::reverse(hole.begin(), hole.end());

}
*/


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
  reorder_island(hole_ids, v);

  std::size_t initial_b_size = b_ids.size();

  // insertion point = just before k
  typename std::vector<int>::iterator insertion_point = b_ids.begin() + k;

  b_ids.insert(insertion_point, hole_ids.begin(), hole_ids.end());

  assert(*(b_ids.begin() + i) == b_ids[i] );
  assert(b_ids.size() == initial_b_size + hole_ids.size());

}

template<typename PointRange>
void join_domain_i(const Domain<PointRange>& domain, Domain<PointRange>& new_domain,
                  const int& i, const int& v, const int& k)
{
  PointRange test;
  typedef std::vector<int> Ids;
  Ids id_set = domain.b_ids;
  Ids hole_ids = domain.h_ids; // for now assume just one.

  merge_id_sets(id_set, i, v, k, hole_ids);

  new_domain.b_ids = id_set;

}




// main loop //
// --------- //
template <typename PointRange, typename WeightMap, typename LambdaMap>
void processDomain(Domain<PointRange>& domain, const int& i, const int& k, std::size_t& count,
                   WeightMap& w_map, LambdaMap& l_map)
{
  // (i, k) = acccess edge
  // v = trird vertex
  int v = 0;

  // domains consisting of only one edge
  if(domain.b_ids.size() == 2)
    return;

  // base case
  if(domain.b_ids.size() == 3 && domain.holes_list.empty())
  {
    //calc weight
    count++;
    return;
  }
  assert(domain.b_ids.size() >= 3);


  //CASE I - if there are islands, join until there are no islands.
  //for (auto point_3 : domain.holeVertices) // todo: check holeVertices
  for(int pid : domain.h_ids)
  {
    //std::cout << "i= " << domain.boundary[i] << " k= " << domain.boundary[k] << std::endl;
    //std::cout << "v = " << point_3 << std::endl;

    std::cout << "i= " << i << " k= " << k << std::endl;
    std::cout << "pid= " << pid << std::endl;

    Domain<PointRange> D1;
    join_domain_i(domain, D1, i, pid, k);

    // get a new e_D
    std::pair<int, int> e_D1 = D1.get_access_edge_ids();

    processDomain(D1, e_D1.first, e_D1.second, count, w_map, l_map);
    v++;

  }


  // CASE II
  v = 0; // temp: index to boundary vertices
  //for(auto point_3 : domain.boundary)
  for(int pid : domain.b_ids)
  {

    //std::cout << "i= " << domain.boundary[i] << " k= " << domain.boundary[k] << std::endl;
    //std::cout << "v = " << point_3 << std::endl;

    std::cout << "i= " << i << " k= " << k << std::endl;
    std::cout << "pid= " << pid << std::endl;

    // avoid source, target of e_D
    if(pid == i || pid == k)
    {
      //++v;
      //std::cout << " point aborted" << std::endl;
      continue;
    }





    // split to two sub-domains
    Domain<PointRange> D1;
    Domain<PointRange> D2;

    // essentially splitting boundaries - change that maybe to work on boundaries directly
    split_domain_i(domain, D1, D2, i, pid, k);
    // D1, D2 have just new boundaries - no hole information.

    // get new access edges for each
    std::pair<int, int> e_D1 = D1.get_access_edge_ids();
    std::pair<int, int> e_D2 = D2.get_access_edge_ids();


    // assign all combination of holes to subdomains and process each pair
    Phi partition_space;
    do_permutations(domain.holes_list, partition_space);
    
    if(partition_space.empty())
    {
      // when the domain has been joined so that there is no holes inside
      processDomain(D1, e_D1.first, e_D1.second, count, w_map, l_map);
      processDomain(D2, e_D2.first, e_D2.second, count, w_map, l_map);
    }
    else
    {
      // when form a t with a vertex on the boundary of a domain with holes
      for(std::size_t p = 0; p < partition_space.size(); ++p)
      {
        std::vector<int> lholes = partition_space.lsubset(p); // vector<int>
        std::vector<int> rholes = partition_space.rsubset(p);

        for(int lh : lholes)
          D1.add_hole(domain.holes_list[lh]);

        for(int rh : rholes)
          D2.add_hole(domain.holes_list[rh]);

        processDomain(D1, e_D1.first, e_D1.second, count, w_map, l_map);
        processDomain(D2, e_D2.first, e_D2.second, count, w_map, l_map);
      }

    }

    ++v; // take next point

  }


}








} // namespace internal
} // namsepace CGAL





#endif // ISLAND_TRIANGULATE_HOLE_POLYLINE_H
