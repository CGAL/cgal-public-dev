#ifndef ISLAND_TRIANGULATE_HOLE_POLYLINE_H
#define ISLAND_TRIANGULATE_HOLE_POLYLINE_H

#include <vector>
#include <CGAL/Combination_enumerator.h>

namespace CGAL {
namespace internal {



struct Phi
{

  typedef std::pair<std::vector<int>, std::vector<int>> Sub_domains_pair;

  void put(std::vector<int>& left, std::vector<int>& right)
  {
    sub_domains_list.push_back(std::make_pair(left, right)); // preallocate list
  }

  std::vector<Sub_domains_pair> sub_domains_list;
};

void print(std::vector<int> &v)
{
  for(int i=0; i<v.size(); ++i)
  {
    std::cout << v[i] <<" ";
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

void add_subsets(std::vector<int>& D1, std::vector<int>& D2, std::vector<int>& sumD)
{
  // concatenate
  sumD.resize(D1.size() + D2.size());
  sumD.insert(sumD.begin(), D1.begin(), D1.end());
  sumD.insert(sumD.end(), D2.begin(), D2.end());
}

void split_domain(const std::pair<int, int>& range, std::vector<int>& left, std::vector<int>& right,
                  const int& i, const int& v, const int& k)
{

  int next_v = v;

  // left subset
  left.push_back(next_v);
  while (next_v != i) {

    if(next_v == range.second)
      next_v = range.first;
    else
      next_v++;

    left.push_back(next_v);
  }
  left.push_back(i);

  // right subset
  next_v = k;
  right.push_back(next_v);
  while (next_v != v) {

    if(next_v == range.second)
      next_v = range.first; // =0
    else
      next_v++;

    right.push_back(next_v);
  }
  left.push_back(v); // convenction: add the first one
}




template <typename PointRange>
void create_subsets(const PointRange& boundary, PointRange& hole)
{

  // e_D (i, k)
  const int i = 1;
  const int k = 2;

  // trird vertex
  const int v = 4;

  std::pair<int, int> range(0, boundary.size() - 1 - 1); // last = first

  std::vector<int> left;
  std::vector<int> right;


  split_domain(range, left, right, i, v, k);

  std::cout << "left: \n";
  print(left);
  std::cout << "right: \n";
  print(right);



}









} // namespace internal
} // namsepace CGAL





#endif // ISLAND_TRIANGULATE_HOLE_POLYLINE_H
