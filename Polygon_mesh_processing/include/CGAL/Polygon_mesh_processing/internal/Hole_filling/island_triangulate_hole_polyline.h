// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Konstantinos Katrioplas

#ifndef CGAL_PMP_INTERNAL_HOLE_FILLING_ISLAND_TRIANGULATE_HOLE_POLYLINE_H
#define CGAL_PMP_INTERNAL_HOLE_FILLING_ISLAND_TRIANGULATE_HOLE_POLYLINE_H

#include <vector>
#include <limits>
#include <CGAL/Combination_enumerator.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>


namespace CGAL {
namespace internal {


// Domain structure //
// ---------------- //
// todo: take PointRange out
template <typename PointRange>
struct Domain
{

  Domain() {}

  // boundary should always be given without the stupid last(first) one.
  // used only in tests
  Domain(const PointRange& boundary) : boundary(boundary)
  {
    std::size_t n_p = boundary.size();
    if(boundary.size() > 0)
      CGAL_assertion(boundary[n_p - 1] != boundary[0]);
  }


  // constructor with indices
  Domain(const std::vector<int>& ids) : b_ids(ids) {}

  void clear_islands()
  {
    islands_list.clear();
  }

  void add_hole(const std::vector<int>& ids) // to change
  {
    islands_list.push_back(ids);
    //h_ids.insert(h_ids.end(), ids.begin(), ids.end());
  }

  bool is_empty()
  {
    return b_ids.size() == 2;
  }

  bool has_islands()
  {
    return !islands_list.empty();
  }

  void add_islands(const std::vector<std::vector<int>> islands)
  {
    assert(this->islands_list.empty());
    islands_list = islands; // to make data private
  }

  void add_islands(const Domain<PointRange> domain,
                   const std::vector<int> island_ids)
  {
    assert(this->islands_list.empty());
    for(int id : island_ids)
    {
      islands_list.push_back(domain.islands_list[id]);
    }
  }

  std::pair<int, int> get_access_edge() const
  {
    CGAL_assertion(b_ids.size() >= 2);
    const int i = b_ids.front();
    const int k = b_ids.back();
    return std::make_pair(i, k);
  }


  //data
  PointRange boundary; // not used in the main algorithm
  std::vector<int> b_ids;
  std::vector<std::vector<int> > islands_list;
};


template <typename T>
void print(T &v)
{
  for(int i=0; i<v.size(); ++i)
  {
    std::cout << v[i] << " ";//<< std::endl;
  }
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
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < sub_domains_list.size());
    return sub_domains_list[i].first;
  }

  std::vector<int> rsubset(const int i)
  {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < sub_domains_list.size());
    return sub_domains_list[i].second;
  }

  std::vector<Sub_domains_pair> sub_domains_list;
};


void do_permutations(std::vector<std::vector<int>>& island_list, Phi& subsets)
{
  if(island_list.empty())
    return;

  std::vector<int> hs;

  for(int n = 0; n < island_list.size(); ++n)
    hs.push_back(n);

  const int first = hs.front();
  const int last = hs.back();
  //std::sort(hs.begin(), hs.end()); // already sorted

  for(int s = 0; s <= hs.size(); ++s) // s = number of islands on one (left) side
  {
    std::vector<int> p1(s);
    std::vector<int> p2(island_list.size() - s);

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

      CGAL_assertion(p1.size() == s);
      CGAL_assertion(p2.size() == hs.size() - s);

      print(p1); std::cout << "-- "; print(p2); std::cout << std::endl;

      subsets.put(p1, p2);
    }

  }
}

// split //
// ----- //

template <typename PointRange>
void split_domain_case_2(const Domain<PointRange>& init_domain,
                         Domain<PointRange>& left_dom, Domain<PointRange>& right_dom,
                         const int i, std::vector<int>::const_iterator it, const int k)
{
  typedef std::vector<int> Ids;
  const Ids& ids = init_domain.b_ids;
  const int pid = *it;

  // i, k indices of access edge = first and last
  // FIXME: as soon as there is a duplicate vertex on the boundary (due to a case I split) only one copy of the attached will be considered
  // fixed: passing iterator to the function to avoid confusion between duplicates

  left_dom.b_ids.assign(ids.begin(), it + 1);
  right_dom.b_ids.assign(it, ids.end());

  CGAL_assertion(left_dom.b_ids.front() == i);
  CGAL_assertion(left_dom.b_ids.back() == pid);
  CGAL_assertion(right_dom.b_ids.front() == pid);
  CGAL_assertion(right_dom.b_ids.back() == k);
}

void rotate_island_vertices(std::vector<int>& ids, const int& position)
{
  // 1) find v's position
  //std::vector<int>::iterator it = find(ids.begin(), ids.end(), v);
  //CGAL_assertion(it != ids.end());

  // get an iterator at the position on the island_ids
  CGAL_assertion(position >=0);
  CGAL_assertion(position < ids.size());

  std::vector<int>::iterator it = ids.begin() + position;

  // 2) rotate by the third vertex of t
  const int v = *it;
  std::rotate(ids.begin(), it, ids.end()); // not working
  CGAL_assertion(ids.front() == v);

  // 3) add the first removed element
  ids.push_back(v);
}

void merge_island_and_boundary(std::vector<int>& b_ids,
                               const int i, const int k,
                               std::vector<int>& island_ids, const int& position)
{
  std::size_t initial_b_size = b_ids.size();
  const int v = island_ids[position];

  rotate_island_vertices(island_ids, position);

  // insertion position = just after k
  // k is at position n - 1 = last element.
  // just append at the end - i is the first point on b_ids
  // and k is the last. t triangle is (i, v, k)
  typename std::vector<int>::iterator insertion_point = b_ids.end();
  b_ids.insert(insertion_point, island_ids.begin(), island_ids.end());

  CGAL_assertion(b_ids[initial_b_size - 1] == k);
  CGAL_assertion(b_ids[0] == i);
  CGAL_assertion(b_ids[initial_b_size] == v);
  CGAL_assertion(b_ids[b_ids.size() - 1] == v);
  CGAL_assertion(b_ids.size() == initial_b_size + island_ids.size());
}

template<typename PointRange>
void split_domain_case_1(const Domain<PointRange>& domain, Domain<PointRange>& D1, Domain<PointRange>& D2,
                         const int i, const int k, std::vector<int> island_ids, const int& position)
{
  typedef std::vector<int> Ids;

  // get boundary ids
  Ids id_set1(domain.b_ids.begin(), domain.b_ids.end());
  Ids id_set2(id_set1);

  // get island ids from the domains's list
  //CGAL_assertion(id < domain.islands_list.size());
  //std::vector<int> ids(domain.islands_list[id]);

  // create two sets - one with reversed orientation
  Ids island_ids1(island_ids.begin(), island_ids.end());
  Ids island_ids2(island_ids.rbegin(), island_ids.rend()); // reversed orientation

  // position points to the correct element of the list
  // merge once with input island
  merge_island_and_boundary(id_set1, i, k, island_ids1, position);
  D1.b_ids = id_set1;

  // merge again with island with reversed orientation
  merge_island_and_boundary(id_set2, i, k, island_ids2, position);
  D2.b_ids = id_set2;
}

// overload + for pairs of doubles
const std::pair<double, double> add_weights(const std::pair<double, double>& p1,
                                            const std::pair<double, double>& p2)
{

  const double angle = p1.first < p2.first ? p1.first : p2.first;
  const double area = p1.second + p2.second;

  return {angle, area};
}


template<typename PointRange, typename WeightCalculator,
         typename WeightTable, typename LambdaTable>
class Triangulate_hole_with_islands
{
  typedef typename WeightCalculator::Weight Weight;
  typedef std::vector<std::size_t> Triangle;
  typedef std::pair<double, double> Wpair;


public:

  Triangulate_hole_with_islands(const Domain<PointRange>& domain,
                                const PointRange& allpoints,
                                WeightTable& W,
                                LambdaTable& l,
                                const WeightCalculator & WC)
    : points(allpoints)
    , W(W)
    , lambda(l)
    , domain(domain)
    , WC(WC)
  {}

  std::size_t do_triangulation(const int i, const int k, std::vector<Triangle>& triangles,
                               std::size_t& count)
  {
    init_triangulation();


    process_domain(domain, std::make_pair(i, k), triangles, count);



    std::cout << std::endl;
    std::cout << "Number of triangles collected: " << triangles.size() << std::endl;

    // will be removed, keep for testing
    // different number before and after indicates bug
    sort(triangles.begin(), triangles.end());
    triangles.erase(unique(triangles.begin(), triangles.end()), triangles.end());

    std::cout << "Number of unique triangles: " << triangles.size() << std::endl;

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
    // will have to include all ids on islands
    for(auto island : domain.islands_list)
    {
      init_island.insert(init_island.end(), island.begin(), island.end());
    }
  }



  // todo: pass Wpair as a reference - maybe
  const Wpair process_domain(Domain<PointRange> domain, const std::pair<int, int> e_D,
                                   std::vector<Triangle>& triangles,
                                   std::size_t& count)
  {

    std::pair<double, double> best_weight = std::make_pair( // todo: use an alias for this
                                            std::numeric_limits<double>::max(),
                                            std::numeric_limits<double>::max());

    int i = e_D.first;
    int k = e_D.second;


    // avoid non-manifold access edges, return invalid triangulation
    if (i == k)
    {
      #ifdef PMP_ISLANDS_DEBUG
      std::cout << "on domain: ";
      for(int j=0; j<domain.b_ids.size(); ++j)
        std::cout << domain.b_ids[j] << " ";
      std::cout << std::endl;
      std::cout <<"i == k: " << i << "=" << k << " returning invalid triangualtion..." <<std::endl;
      // because assess edge would be degenerate
      #endif

      return std::make_pair( // todo: use an alias for this
                             std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::max());
    }

    // empty domain: adds nothing and is not invalid
    if(domain.b_ids.size() == 2)
      return std::make_pair(0, 0);

    // base case triangle evaluation
    if(domain.b_ids.size() == 3 && !domain.has_islands())
    {
      CGAL_assertion(domain.b_ids[0] == i); // access edge source
      CGAL_assertion(domain.b_ids[2] == k); // access edge target

      int m = domain.b_ids[1]; //third vertex
      const Wpair weight = calc_weight(i, m, k);
      assert(weight.first >= 0);
      assert(weight.second >= 0);

      // return the triangle and its weight
      ++count;
      triangles = {{i, m, k}};
      return weight;
    }

    CGAL_assertion(domain.b_ids.size() >= 3);



    // case I

    // merge each island
    for(std::size_t island_id = 0; island_id < domain.islands_list.size(); ++island_id)
    {


      // local islands are the islands left without the one that is being merged with case I below.
      std::vector<std::vector<int>> local_islands(domain.islands_list);
      local_islands.erase(local_islands.begin() + island_id);



      // take each vertex of this island
      for(int j = 0; j < domain.islands_list[island_id].size(); ++j)
      {

        // point that is being connected to the boundary
        const int pid = domain.islands_list[island_id][j];


        //std::cout << "pid = " << pid << std::endl;
        //std::cin.get();

        assert(std::find(domain.islands_list[island_id].begin(), domain.islands_list[island_id].begin(), pid) !=
               domain.islands_list[island_id].end());


        Domain<PointRange> D1, D2;
        // assign the remaining islands the domain(both orientations) that are produced
        D1.add_islands(local_islands);
        D2.add_islands(local_islands);
        // split_domain_case_1 joins the island that pid belongs to
        split_domain_case_1(domain, D1, D2, i, k, domain.islands_list[island_id], j);

        std::pair<int, int> e_D1 = D1.get_access_edge();
        std::pair<int, int> e_D2 = D2.get_access_edge();
        std::vector<Triangle> triangles_D1, triangles_D2;

        const Wpair w_D1 = process_domain(D1, e_D1, triangles_D1, count);
        const Wpair w_D2 = process_domain(D2, e_D2, triangles_D2, count);

        CGAL_assertion(w_D1.first <= 180);
        CGAL_assertion(w_D2.first <= 180);


        // evaluate triangulations

        // choose the best orientation
        if(w_D1 < w_D2)
        {
          // calculate w(t) & add w(t) to w_D1
          const Wpair weight_t = calc_weight(i, pid, k);
          const Wpair w = add_weights(w_D1, weight_t);

          if(w < best_weight)
          {
            // update the best weight
            best_weight = w;

            // keep only the best
            triangles.clear();

            // add t to triangles_D2 and return them
            Triangle t = {i, pid, k};
            triangles.insert(triangles.begin(), triangles_D1.begin(), triangles_D1.end());
            triangles.insert(triangles.end(), t);
         }
        }
        else
        {
          // calculate w(t) & add w(t) to w_D2
          const Wpair weight_t = calc_weight(i, pid, k);
          const Wpair w = add_weights(w_D2, weight_t);

          if(w < best_weight)
          {
            // update the best weight
            best_weight = w;

            // keep only the best
            triangles.clear();

            // add t to triangles_D2 and return them
            Triangle t = {i, pid, k};
            triangles.insert(triangles.begin(), triangles_D2.begin(), triangles_D2.end());
            triangles.insert(triangles.end(), t);
         }
        }

        // does not return: need to evaluate permutations of islands for case 2 splits

        // triangles and best weight for this level have been calculated and
        // compared with those from the case II splitting below, which occurs without
        // case I before.

      } // pid : domains.all_h_ids - case 1 split


      std::cout << "reached end of islands for case 1 " << std::endl;

    } // end list of islands




    // case 2

    // avoid begin and end of the container which is the source and target of the access edge
    for(std::vector<int>::iterator pid_it = domain.b_ids.begin() + 1; pid_it != domain.b_ids.end() - 1; ++pid_it)
    {

      const int pid = *pid_it;

#ifdef PMP_ISLANDS_DEBUG
      std::cout << "on domain: ";
      for(int j=0; j<domain.b_ids.size(); ++j)
        std::cout << domain.b_ids[j] << " ";
      std:: cout <<", pid: " << pid << ", splitting..." <<std::endl;
#endif

      // a triangle that has islands is considered
      // an invalid triangulation because it disconnects boundary and island
      if(domain.b_ids.size() == 3 && domain.has_islands())
        return std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());


      Domain<PointRange> D1, D2;
      // split_domain_case_2 splits the domain to the boundary by creating 2 subdomains
      // D1, D2 have just new boundaries - no island information.
      split_domain_case_2(domain, D1, D2, i, pid_it, k);

      CGAL_assertion(D1.b_ids[0] == i);
      CGAL_assertion(D2.b_ids[D2.b_ids.size() - 1] == k);


      // get access edge: it used in any of the following cases at this level
      std::pair<int, int> e_D1 = D1.get_access_edge();
      std::pair<int, int> e_D2 = D2.get_access_edge();


      // --- assign islands to each subdomain ---
      Wpair w_D12 = std::make_pair(std::numeric_limits<double>::max(),
                                  std::numeric_limits<double>::max());
      std::vector<Triangle> triangles_D1, triangles_D2;

      CGAL_assertion(!D1.has_islands() && !D2.has_islands());

      // if domain does not have islands, then just the main loop is called for each.
      // This is evaluated after case I, since case I stops happening
      // only after there are no islands left.
      if(!domain.has_islands())
      {
        // no islands have been assigned to D1 and D2 after the case_2
        w_D12 = add_weights(process_domain(D1, e_D1, triangles_D1, count),
                            process_domain(D2, e_D2, triangles_D2, count) );
      }

      // if domain does have islands, then we don't need to parition if
      // one of D1,D2 is empty: all islands go to the not empty
      else
      {
        if(D1.is_empty())
        {
          // assign all left: local_islands
          D2.add_islands(domain.islands_list);
          w_D12 = process_domain(D2, e_D2, triangles_D2, count);
        }
        else
        {
          if(D2.is_empty())
          {
            // assign all left: local_islands
            D1.add_islands(domain.islands_list);
            w_D12 = process_domain(D2, e_D2, triangles_D1, count);
          }
          // if there are islands in domain, then take all combinations of them on
          // each domain
          else
          {
            assert(!D1.is_empty());
            assert(!D2.is_empty());
            assert(domain.has_islands());

            Phi partition_space;
            do_permutations(domain.islands_list, partition_space);

            for(int p = 0; p < partition_space.size(); ++p)
            {
              // for each permutation, start over
              D1.clear_islands();
              D2.clear_islands();

              // get indices to islands
              std::vector<int> islands_D1 = partition_space.lsubset(p);
              std::vector<int> islands_D2 = partition_space.rsubset(p);

              // assign
              D1.add_islands(domain, islands_D1);
              D2.add_islands(domain, islands_D2);


              std::vector<Triangle> local_triangles_D1, local_triangles_D2;
              const Wpair local_w_D1 = process_domain(D1, e_D1, local_triangles_D1, count);
              const Wpair local_w_D2 = process_domain(D2, e_D2, local_triangles_D2, count);

              if (add_weights(local_w_D1,local_w_D2) < w_D12)
              {
                w_D12=add_weights(local_w_D1, local_w_D2);
                triangles_D1.swap(local_triangles_D1);
                triangles_D2.swap(local_triangles_D2);
              }
            }
          }
        }
      }


      // calculate w(t)
      const Wpair weight_t = calc_weight(i, pid, k);
      ++count;
      // add it to its subdomains
      const Wpair w = add_weights(w_D12, weight_t);

      // those from case I splits
      if(w < best_weight)
      {
        // update the best weight
        best_weight = w;

        // since this triangulation is better, get rid of the one collected so far
        // at this level before adding the better one.
        triangles.clear();

        // joint subdomains with t and return them
        Triangle t = {i, pid, k};
        triangles.insert(triangles.begin(), triangles_D1.begin(), triangles_D1.end());
        triangles.insert(triangles.end(), triangles_D2.begin(), triangles_D2.end());
        triangles.insert(triangles.end(), t);

        #ifdef PMP_ISLANDS_DEBUG
        std::cout << "-->triangles" << std::endl;
        for(int t = 0; t < triangles.size(); ++t)
        {
          Triangle tr = triangles[t];
          std::cout << "-->" << tr[0] << " " << tr[1] << " " << tr[2] << std::endl;
        }
        #endif
      } // w < best weight


    } // case 2 splits






    // useful when return for case II splits that have followed case I,
    // so as to compare different case I splits.
    return best_weight;
  }


  bool are_vertices_on_island(const int i, const int m, const int k)
  {
    std::vector<int>::iterator it1, it2, it3;
    it1 = std::find(init_island.begin(), init_island.end(), i);
    it2 = std::find(init_island.begin(), init_island.end(), m);
    it3 = std::find(init_island.begin(), init_island.end(), k);
    return (it1 != init_island.end()) && (it2 != init_island.end()) && (it3 != init_island.end()) ?  true : false;
  }

  bool are_vertices_on_island(std::vector<int> ids)
  {
    // this works for only 3 points
    CGAL_assertion(ids.size() == 3);

    const int i = ids[0];
    const int m = ids[1];
    const int k = ids[2];
    std::vector<int>::iterator it1, it2, it3;
    it1 = std::find(init_island.begin(), init_island.end(), i);
    it2 = std::find(init_island.begin(), init_island.end(), m);
    it3 = std::find(init_island.begin(), init_island.end(), k);
    return (it1 != init_island.end()) && (it2 != init_island.end()) && (it3 != init_island.end()) ?  true : false;
  }


  const Wpair calc_weight(const int i, const int m, const int k)
  {

    if(are_vertices_on_island(i, m, k))
    {

#ifdef PMP_ISLANDS_DEBUG
      std::cout << "vertices on island, invalid triangulation" << std::endl;
#endif

      return std::make_pair( // todo: use an alias for this
                             std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::max());
    }

    // to remove this and use a new function object
    const Weight& w_t = WC(points, Q, i, m, k, lambda);

    double angle = w_t.w.first;
    double area = w_t.w.second;

    // temp: handle degenerate edges - will be taken care with a new structure for the weight
    // which will produce the numeric limit instead of -1
    if(angle == -1)
      angle = std::numeric_limits<double>::max();
    if(area == -1)
      area = std::numeric_limits<double>::max();

    assert(angle >= 0);
    assert(area >= 0);

    return std::make_pair(angle, area);
  }






  // data
  const PointRange& points;
  const PointRange Q; // empty - to be optimized out

  // todo: use weight tables to avoid calc same weight again
  //std::set<std::vector<std::size_t>> memoized;

  WeightTable W; // to be removed
  LambdaTable lambda; // to be removed

  const Domain<PointRange>& domain;
  const WeightCalculator& WC; // a new object will be used

  // initial island vertices
  std::vector<int> init_island;

};


} // namespace internal
} // namsepace CGAL





#endif // CGAL_PMP_INTERNAL_HOLE_FILLING_ISLAND_TRIANGULATE_HOLE_POLYLINE_H
