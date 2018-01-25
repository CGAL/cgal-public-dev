#ifndef ISLAND_TRIANGULATE_HOLE_POLYLINE_H
#define ISLAND_TRIANGULATE_HOLE_POLYLINE_H

#include <vector>
#include <CGAL/Combination_enumerator.h>

namespace CGAL {
namespace internal {


struct Domain
{


  //data
  //1) boundary
  //2) holes

};

struct Phi
{

  typedef std::pair<std::vector<int>, std::vector<int>> Sub_domains_pair;

  void put(std::vector<int>& left, std::vector<int>& right)
  {
    sub_domains_list.push_back(std::make_pair(left, right)); // preallocate list
  }

  std::vector<Sub_domains_pair> sub_domains_list;
};

template <typename PointRange>
void print(PointRange &v)
{
  for(int i=0; i<v.size(); ++i)
  {
    std::cout << v[i] << std::endl;
  }
  std::cout << std::endl;
}

void do_permutations(const int s, std::vector<int>& hs, Phi& subsets)
{
  const int first = hs.front();
  const int last = hs.back();
  std::sort(hs.begin(), hs.end());

  std::vector<int> p1(s);
  std::vector<int> p2(hs.size() - s);

  if(s == 0)
  {
    subsets.put(p1, hs);

    print(p1); std::cout << "-- "; print(hs); std::cout << std::endl;
    return;
  }

  CGAL::Combination_enumerator<int> permutations(s, first, last + 1);

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

template <typename PointRange>
void split_domain(PointRange& boundary, PointRange& left, PointRange& right,
                  const int& i, const int& v, const int& k)
{

  // take out last(=first)
  boundary.pop_back();

  int next_v = v;

  // left subset
  left.push_back(boundary[next_v]);
  while (next_v != i) {

    if(next_v == boundary.size() - 1) // if we reached the last element
      next_v = 0;
    else
      next_v++;

    left.push_back(boundary[next_v]);
  }
  // add the new last(=first)
  left.push_back(boundary[v]);

  // right subset
  next_v = k;
  right.push_back(boundary[next_v]);
  while (next_v != v) {

    if(next_v == boundary.size() - 1)
      next_v = 0;
    else
      next_v++;

    right.push_back(boundary[next_v]);
  }
  // add the new last(=first)
  right.push_back(boundary[k]);

  assert(left.front() == boundary[v]);
  assert(left.back() == boundary[v]);
  assert(right.front() == boundary[k]);
  assert(right.back() == boundary[k]);
}


template <typename PointRange>
void reorder_island(PointRange& hole, const int& v)
{
  assert(v >= 0);
  assert(v < hole.size());

  // 1) take the last(=first) out
  hole.pop_back();

  // 2) rotate by the third vertex of t
  std::rotate(hole.begin(), hole.begin() + v, hole.end());

  // 3) add the first removed element
  hole.push_back(hole[0]);

  // 4) reverse. Todo: Check and do it iff reversal is needed.
  std::reverse(hole.begin(), hole.end());

}


template <typename PointRange>
void join_domain(PointRange& boundary,
                 const int& i, const int& v, const int& k,
                 PointRange& hole)
{
  // assuming the old boundary is not needed any more, insert on this one.
  //new_domain.resize(boundary.size() + hole.size() - 1); // there are 2 extra points repeated, one on the boundary and one on the hole.

  reorder_island(hole, v);

  std::size_t initial_b_size = boundary.size();

  // insertion point = third vertex of t
  typename PointRange::iterator insertion_point = boundary.begin() + k;

  boundary.insert(insertion_point, hole.begin(), hole.end());

  assert(*(boundary.begin() + i) == boundary[i] );
  assert(boundary.size() == initial_b_size + hole.size());

}

template <typename PointRange>
void test_split_domain(PointRange& boundary)
{
  // e_D (i, k)
  const int i = 1;
  const int k = 2;
  // trird vertex - on the boundary
  const int v = 4;

  PointRange left;
  PointRange right;

  split_domain(boundary, left, right, i, v, k);
  std::cout << "left: \n";
  print(left);
  std::cout << "right: \n";
  print(right);

}

template <typename PointRange>
void test_join_domain(PointRange& boundary, PointRange& hole)
{
  // e_D (i, k)
  const int i = 1;
  const int k = 2;
  // trird vertex - index of hole vertices
  const int v = 1;

  join_domain(boundary, i, v, k, hole);


}



template <typename PointRange>
void create_subsets(PointRange& boundary, PointRange& hole)
{

  test_split_domain(boundary);

  test_join_domain(boundary, hole);

}












} // namespace internal
} // namsepace CGAL





#endif // ISLAND_TRIANGULATE_HOLE_POLYLINE_H