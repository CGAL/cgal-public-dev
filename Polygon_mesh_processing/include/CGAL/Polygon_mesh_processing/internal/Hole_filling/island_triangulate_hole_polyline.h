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

//#define LOCAL_DT_FACES
//#define LOCAL_DT_EDGES

#include <CGAL/Timer.h>

#include <vector>
#include <limits>
#include <unordered_map>
#include <CGAL/Combination_enumerator.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polyline.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <boost/container/flat_set.hpp>

namespace CGAL {
namespace internal {


struct Pair_set
{
  typedef std::pair<void*, bool> return_type;

  Pair_set(std::size_t n)
    : nkeys(n)
    , bs(n*n, false)
  {}

  return_type
  insert(const std::pair<std::size_t, std::size_t>& p)
  {
    std::size_t key = p.first + nkeys * p.second;
    if (bs[key])
      return return_type(NULL, false);
    bs[key]=true;
    return return_type(NULL, true);
  }

  void swap(Pair_set& set2)
  {
    nkeys = set2.nkeys;
    bs = set2.bs;
  }

  typedef std::vector<int> Triangle;

  void add_from_triangles(std::vector<Triangle>& triangles)
  {
    for(const auto t : triangles)
    {
      int i = t[0];
      int v = t[1];
      int k = t[2];
      insert(std::make_pair(k, i));
      insert(std::make_pair(i, v));
      insert(std::make_pair(v, k));
    }
  }

  std::size_t nkeys;
  std::vector<bool> bs;
};

// Domain structure //
// ---------------- //
struct Domain
{

  Domain() {}

  // constructor with indices
  // boundary should always be given without the stupid last(first) one.
  Domain(const std::vector<int>& ids) : b_ids(ids) {}

  void clear_islands()
  {
    islands_list.clear();
  }

  void add_island(const std::vector<int>& ids)
  {
    islands_list.push_back(ids);
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
    CGAL_assertion(this->islands_list.empty());
    islands_list = islands;
  }

  void add_islands(const Domain domain,
                   const std::vector<int> island_ids)
  {
    CGAL_assertion(this->islands_list.empty());
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

  // maybe do this in a free function
  std::vector<std::pair<int, int> >  edges() const
  {
    std::vector<std::pair<int, int> > edges;

    CGAL_assertion(b_ids.size() > 0);

    // boundary edges
    for(auto b_it = b_ids.begin() + 1; b_it != b_ids.end(); ++b_it)
      edges.push_back(std::make_pair(*b_it, *(b_it-1)));
    edges.push_back(std::make_pair(*b_ids.begin(), *(b_ids.end()-1)));

    // island edges
    for(auto island : islands_list)
    {
      for(auto i_it = island.begin() + 1; i_it != island.end(); ++i_it)
        edges.push_back(std::make_pair(*i_it, *(i_it-1)));
      edges.push_back(std::make_pair(*island.begin(), *(island.end()-1)));
    }

    return edges;
  }

  // data
  std::vector<int> b_ids;
  std::vector<std::vector<int> > islands_list;
};


template<class T>
class Weights_cache_map {
public:
  Weights_cache_map(int n, const T& default_) : n(n), default_(default_) { }

  void put(std::pair<int, int> e_D, int v, const T& t) {
    //CGAL_assertion(bound_check(e_D.first, e_D.second, v));

    table[e_D] = std::make_pair(v, t);
  }

  const T& get_best_weight(std::pair<int, int> e_D) const {
    CGAL_assertion(bound_check(e_D.first, e_D.second));

    typename Map::const_iterator it = table.find(e_D);
    if(it != table.end())
      return it->second.second;
    else
    {
      // is the default useful?
      //assert(false);
      return default_;
    }
  }

  const int get_access_vertex(const std::pair<int, int>& e_D) const
  {
    CGAL_assertion(bound_check(e_D.first, e_D.second));

    typename Map::const_iterator it = table.find(e_D);
    if(it != table.end())
      return it->second.first;
    else
    {
      // is the default useful?
      //assert(false);
      return -1;
    }
  }

  bool exists(const std::pair<int, int>& e_D)
  {
    return table.find(e_D) != table.end() ? true : false;
  }

  int n;
private:
  // map stores <e_D, <access_vertex, best_weight> >
  typedef std::map< std::pair<int, int>, std::pair<int, T> > Map;
  bool bound_check(int i, int m, int j) const {
    CGAL_assertion(i >= 0 && i < n);
    CGAL_assertion(j >= 0 && j < n);
    CGAL_assertion(m >= 0 && m < n);
    return true;
  }
  bool bound_check(int i, int j) const {
    CGAL_assertion(i >= 0 && i < n);
    CGAL_assertion(j >= 0 && j < n);
    return true;
  }
  Map table;
  T default_;
};


template<class T>
class Best_triangles_table_map {
public:
  Best_triangles_table_map(int n, const T& default_) : n(n), default_(default_) { }

  void put(std::pair<int, int> e_D, const T& t) {
    CGAL_assertion(bound_check(e_D.first, e_D.second));

    table[e_D] = t;
  }

  const T& get_best_triangles(std::pair<int, int> e_D) const {
    CGAL_assertion(bound_check(e_D.first, e_D.second));

    typename Map::const_iterator it = table.find(e_D);
    if(it != table.end())
      return it->second;
    else
    {
      // is the default useful?
      assert(false);
      return default_;
    }
  }

  bool exists(const std::pair<int, int>& e_D)
  {
    return table.find(e_D) != table.end() ? true : false;
  }

  int n;
private:
  // map stores <e_D, best_triangles >
  typedef std::map< std::pair<int, int>, T> Map;
  bool bound_check(int i, int j) const {
    CGAL_assertion(i >= 0 && i < n);
    CGAL_assertion(j >= 0 && j < n);
    return true;
  }
  Map table;
  T default_;
};





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

  std::vector<int> lsubset(const std::size_t i)
  {
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < sub_domains_list.size());
    return sub_domains_list[i].first;
  }

  std::vector<int> rsubset(const std::size_t i)
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

  for(std::size_t n = 0; n < island_list.size(); ++n)
    hs.push_back(n);

  const int first = hs.front();
  const int last = hs.back();
  //std::sort(hs.begin(), hs.end()); // already sorted

  for(std::size_t s = 0; s <= hs.size(); ++s) // s = number of islands on one (left) side
  {
    std::vector<int> p1(s);
    std::vector<int> p2(island_list.size() - s);

    if(s == 0)
    {
      subsets.put(p1, hs);
      continue;
    }

    CGAL::Combination_enumerator<int> permutations(s, first, last+1);
    while(!permutations.finished())
    {
      for(std::size_t i=0; i<s; ++i)
      {
        p1[i] = permutations[i];
      }

      ++permutations;

      std::sort(p1.begin(), p1.end());
      std::set_symmetric_difference(p1.begin(), p1.end(), hs.begin(), hs.end(),
                                    p2.begin());
      CGAL_assertion(p1.size() == s);
      CGAL_assertion(p2.size() == hs.size() - s);
      subsets.put(p1, p2);
    }

  }
}

// split //
// ----- //

void split_domain_case_2(const Domain& init_domain,
                         Domain& left_dom, Domain& right_dom,
                         const int i, std::vector<int>::const_iterator it, const int k)
{
  typedef std::vector<int> Ids;
  const Ids& ids = init_domain.b_ids;
  const int pid = *it;

  // i, k indices of access edge = first and last
  // passing iterator to avoid confusion between duplicates

  left_dom.b_ids.assign(ids.begin(), it + 1);
  right_dom.b_ids.assign(it, ids.end());

  CGAL_assertion(left_dom.b_ids.front() == i);
  CGAL_assertion(left_dom.b_ids.back() == pid);
  CGAL_assertion(right_dom.b_ids.front() == pid);
  CGAL_assertion(right_dom.b_ids.back() == k);
}

void merge_island_and_boundary(std::vector<int>& b_ids,
                               const int i, const int k,
                               std::vector<int>& island_ids)
{
  CGAL_assertion_code( std::size_t initial_b_size = b_ids.size() );

  // insertion position = just after k
  // k is at position n - 1 = last element.
  // just append at the end - i is the first point on b_ids
  // and k is the last. t triangle is (i, v, k)
  std::vector<int>::iterator insertion_point = b_ids.end();
  b_ids.insert(insertion_point, island_ids.begin(), island_ids.end());

  CGAL_assertion(b_ids[initial_b_size - 1] == k);
  CGAL_assertion(b_ids[0] == i);
  CGAL_assertion(b_ids[initial_b_size] == island_ids[0]);
  CGAL_assertion(b_ids[b_ids.size() - 1] == island_ids[island_ids.size() - 1]);
  CGAL_assertion(b_ids.size() == initial_b_size + island_ids.size());
}

void split_domain_case_1(const Domain& domain, Domain& D1, Domain& D2,
                         const int i, const int k, std::vector<int> island_ids, const int& position)
{
  // position points to the vertex on the island that is being joined to the boundary

  typedef std::vector<int> Ids;

  // get boundary ids
  Ids id_set1(domain.b_ids.begin(), domain.b_ids.end());
  Ids id_set2(id_set1);

  // rotate by the position
  std::vector<int> local_island_ids(island_ids.begin() + position, island_ids.end());
  local_island_ids.insert(local_island_ids.end(),
                          island_ids.begin(), island_ids.begin() + position);

  // add the first to the end: island is a closed
  local_island_ids.push_back(local_island_ids[0]);

  // create two sets - one with reversed orientation
  Ids island_ids1(local_island_ids.begin(), local_island_ids.end());
  Ids island_ids2(local_island_ids.rbegin(), local_island_ids.rend());

  // merge once with input island
  merge_island_and_boundary(id_set1, i, k, island_ids1);
  // merge again with island with reversed orientation
  merge_island_and_boundary(id_set2, i, k, island_ids2);

  D1.b_ids = id_set1;
  D2.b_ids = id_set2;

}

const std::pair<double, double> add_weights(const std::pair<double, double>& p1,
                                            const std::pair<double, double>& p2)
{

  const double angle = p1.first < p2.first ? p1.first : p2.first;
  const double area = p1.second + p2.second;

  return {angle, area};
}


struct Compare_triangles
{
  typedef std::vector<int> Triangle;

  bool operator() (const Triangle& t1, const Triangle& t2)
  {
    CGAL_assertion(t1.size() == 3);
    CGAL_assertion(t2.size() == 3);

    // t1, t2 are sorted

    if(t1[0] == t2[0] && t1[1] == t2[1])
      return t1[2] < t2[2];

    if(t1[0] == t2[0])
      return t1[1] < t2[1];

    return t1[0] < t2[0];

  }
};





template<typename PointRange, typename LambdaTable, typename WeightCalculator>
class Triangulate_hole_with_islands
{
  typedef typename WeightCalculator::Weight Weight; // min_max_angle_area
  typedef std::vector<int> Triangle;
  typedef std::pair<double, double> Wpair;
  typedef typename Kernel_traits<
    typename std::iterator_traits<
      typename PointRange::const_iterator>::value_type>::Kernel Kernel;
  typedef Triangulation_vertex_base_with_info_3<int, Kernel> Vb;
  typedef Triangulation_data_structure_3<Vb> TDS;

  typedef Delaunay_triangulation_3<Kernel, TDS> DT3;

  //typedef CGAL::internal::Lookup_table_map<int> LambdaTable;
  //typedef CGAL::internal::Lookup_table_map<Weight> WeightsTable;
  //typedef Weights_cache_map<Weight> WeightsTable;

  // use my weight
  typedef Weights_cache_map<Wpair> WeightsTable;
  typedef Best_triangles_table_map<std::vector<Triangle>> TrianglesTable;

public:

  Triangulate_hole_with_islands(const Domain& domain,
                                const PointRange& allpoints,
                                LambdaTable& lambda,
                                const WeightCalculator & WC,
                                const int n,
                                bool use_dt3 = false,
                                bool ci_orientation = false)
    : points(allpoints)
    , domain(domain)
    , lambda(lambda)
    , WC(WC)
    , n(n)
    , correct_island_orientation(ci_orientation)
  {
    if (use_dt3) build_dt3();
    weights_cache = NULL;
  }

  void build_dt3()
  {
    // collect initial boundary edges
    std::vector<std::pair<typename Kernel::Point_3, int> > points_and_indices;
    points_and_indices.reserve(points.size());
    int i=0;
    for(const typename Kernel::Point_3& p : points)
      points_and_indices.push_back( std::make_pair(p, i++) );

    // get Delaunay
    dt3.insert(points_and_indices.begin(), points_and_indices.end());
    dt3_vertices.resize(points.size());
    for(typename DT3::Finite_vertices_iterator vit=dt3.finite_vertices_begin(),
                                               end=dt3.finite_vertices_end(); vit!=end; ++vit)
    {
      dt3_vertices[vit->info()]=vit;
    }

    // collect dt3_edges
    #ifdef LOCAL_DT_EDGES
    for(typename DT3::Finite_edges_iterator eit=dt3.finite_edges_begin(),
                                            end=dt3.finite_edges_end(); eit!=end; ++eit)
    {
      int v0 = eit->first->vertex(eit->second)->info();
      int v1 = eit->first->vertex(eit->third)->info();
      dt3_edges.insert(std::make_pair(v0, v1));
    }
    #endif

    #ifdef LOCAL_DT_FACES
    // collect dt3 faces
    for(typename DT3::Finite_facets_iterator fit=dt3.finite_facets_begin(),
                                                end=dt3.finite_facets_end(); fit!=end; ++fit)
     {
        Triangle triangle(3);
        for(int i = 1; i <= 3; ++i)
        {
          int indx = fit->first->vertex((fit->second + i) % 4)->info();
          triangle[i-1] = indx;
        }

        std::sort(triangle.begin(), triangle.end());
        dt3_faces.insert(triangle);
      }

    /*
    std::ofstream outf("data/facesDT.txt");
    for(auto it = dt3_faces.begin(); it != dt3_faces.end(); ++it)
    {
      auto t = *it;
      outf << t[0] << " " << t[1] << " " << t[2] << "\n";
    }
    outf.close();
    */

    #endif

    //edges are not embedded in the DT3, clearing it
   #ifndef LOCAL_DT_EDGES
   if (!can_use_dt3())
   #else
   if (!can_use_dt3())
   #endif
    {
      dt3_vertices.clear();
      dt3.clear();
    }
  }

  bool skip_facet(int i, int j, int k) const
  {
    if (dt3_vertices.empty()) return false;

    #ifndef LOCAL_DT_FACES
    typename DT3::Cell_handle ch;
    int tmp_i, tmp_j, tmp_k;
    return !dt3.is_facet(dt3_vertices[i],
                         dt3_vertices[j],
                         dt3_vertices[k],
                         ch, tmp_i, tmp_j, tmp_k);
    #endif

    #ifdef LOCAL_DT_FACES
    Triangle t = { i, j, k };
    std::sort(t.begin(), t.end());
    return !std::binary_search(dt3_faces.begin(), dt3_faces.end(), t);
    #endif
  }

  bool can_use_dt3() const
  {
    #ifndef LOCAL_DT_EDGES
    typename DT3::Cell_handle ch;
    int tmp_i, tmp_j;

    // for (i, j index of edge vertices on the boundary of the domain and on island)
    std::vector<std::pair<int, int> > edges = domain.edges(); //crappy copying
    for(std::pair<int, int> e : edges)
    {
      int i = e.first;
      int j = e.second;
      if (!dt3.is_edge (dt3_vertices[i], dt3_vertices[j], ch, tmp_i, tmp_j))
        return false;
    }
    return true;

    #else
    // for (i, j index of edge vertices on the boundary of the domain and on island)
    std::vector<std::pair<int, int> > edges = domain.edges(); //crappy copying
    for(std::pair<int, int> e : edges)
    {
      std::pair<int, int> e_opp(std::make_pair(e.second, e.first));

      if(!std::binary_search(dt3_edges.begin(), dt3_edges.end(), e) &&
         !std::binary_search(dt3_edges.begin(), dt3_edges.end(), e_opp) )
          return false;
    }
    return true;
    #endif
  }

  void do_triangulation(const int i, const int k, std::vector<Triangle>& triangles)
  {    
    //boost::container::flat_set< std::pair<int,int> > boundary_edges_picked;
    //std::set< std::pair<int,int> > boundary_edges_picked;
    Pair_set boundary_edges_picked(n);

    // INSERT BOUNDARY EDGES

    // adds all boundary edges to the bep.
    // loop on b_ids + add in boundary_edges_picked  make_pair(b_ids[k],b_ids[k-1])
    for(auto b_it = domain.b_ids.begin() + 1; b_it != domain.b_ids.end(); ++b_it)
      boundary_edges_picked.insert(std::make_pair(*b_it, *b_it-1)); // *b_it - 1 == *(b_it-1) ?
    boundary_edges_picked.insert(std::make_pair(*domain.b_ids.begin(), *(domain.b_ids.end()-1)));

    count_DT_skips = 0;
    count_triangles = 0;
    count_table_skips = 0;
    count_avoiding_beps = 0;

    // boundary edges here contain boundary edges (halfedges)
    process_domain(domain, std::make_pair(i, k), triangles, boundary_edges_picked);

    std::cout << std::endl;
    std::cout << "Number of triangulations avoided with bep: " << count_avoiding_beps << std::endl;
    std::cout << "Number of triangulations not in DT skipped: " << count_DT_skips << std::endl;
    std::cout << "Number of triangulations in Tables skipped: " << count_table_skips << std::endl;
    std::cout << "Possible triangles tested: " << count_triangles << std::endl;
    std::cout << "Number of triangles collected: " << triangles.size() << std::endl;

    // generate offs
    /*
    std::ofstream out("data/output.off");
    out << "OFF\n" << points.size() << " " << triangles.size() << " 0\n";
    for(auto p : points)
      out << p <<"\n";
    for(auto t : triangles)
      out << "3 " << t[0] << " " << t[1] << " " << t[2] << std::endl;
    out.close();
    */

  }

  template <typename PolygonMesh>
  void visualize(PointRange& points, std::vector<std::vector<int>>& polygon_soup,
                 PolygonMesh& mesh)
  {
    //CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygon_soup);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygon_soup, mesh);
  }


private:

  void construct_triangulation(const std::pair<int, int>& e_D,
                               std::vector<Triangle>& triangles)
  {
    CGAL_assertion_code( const int n = weights_cache->n; )
    std::stack<std::pair<int, int> > edges;
    std::set<std::pair<int, int> > used_edges;

    edges.push(e_D);
    while(!edges.empty()) {
      std::pair<int, int> e = edges.top();
      edges.pop();
      CGAL_assertion(e.first >= 0 && e.first < n);
      CGAL_assertion(e.second >= 0 && e.second < n);

      if(used_edges.count(e) != 0)
      {
        continue;
      }

      const int v_D = weights_cache->get_access_vertex(e);
      const auto weight = weights_cache->get_best_weight(e);
      const double half_weight = weight.second;

      if(v_D == -1 || half_weight == std::numeric_limits<double>::max())
      {
        continue;
      }

      CGAL_assertion(v_D >= 0 && v_D < n);
      triangles.push_back({e.first, v_D, e.second});

      // cache the triangles that are being reconstructed
      used_edges.insert(e);

      edges.push(std::make_pair(e.first, v_D));
      edges.push(std::make_pair(v_D, e.second));
    }
  }

  // used in reconstruction with weights
  void gather_boundary_edges(const std::vector<Triangle>& triangles,
                             //std::set< std::pair<int,int> >& boundary_edges_picked)
                             Pair_set& boundary_edges_picked)
  {
    // triangles are stored {i, pid, k}
    for(const Triangle& t : triangles)
    {
      int i = t[0];
      int v = t[1];
      int k = t[2];

      boundary_edges_picked.insert(std::make_pair(k, i));
      boundary_edges_picked.insert(std::make_pair(i, v));
      boundary_edges_picked.insert(std::make_pair(v, k));

      /* // why should there be no conflicts here?
      bool ok = true;
      ok &= boundary_edges_picked.insert(std::make_pair(k, i)).second;
      ok &= boundary_edges_picked.insert(std::make_pair(i, v)).second;
      ok &= boundary_edges_picked.insert(std::make_pair(v, k)).second;
      CGAL_assertion(ok);
      */
    }
  }

  const Wpair process_domain(Domain domain, const std::pair<int, int> e_D,
                             std::vector<Triangle>& triangles,
                             //boost::container::flat_set< std::pair<int,int> >& boundary_edges_picked
                             //std::set< std::pair<int,int> >& boundary_edges_picked
                             Pair_set& boundary_edges_picked
                             )
  {
    // empty domain: adds nothing and is not invalid
    if(domain.b_ids.size() == 2)
      return std::make_pair(0, 0);

    // use cache in islandless domains
    if(weights_cache != NULL && weights_cache->exists(e_D))
    {
      ++count_table_skips;

      // two approaches to restore the triangulation:
      // 1. use function construct_triangulation to create one triangle after the other
      //    using only the access edge and the opposite vertex
      // 2. cache it in table triangulation_cache - not the best idea

      // approach 1
      // reconstruct
      std::vector<Triangle> reconstructed_triangulation;
      construct_triangulation(e_D, reconstructed_triangulation);

      // approach 2
      // restores triangulation from a cache table
      //std::vector<Triangle> cached_triangles = triangulation_cache->get_best_triangles(e_D);

      // triangles is always empty
      CGAL_assertion(triangles.empty()); // todo: swap instead of insert

      // use approach 1
      triangles.insert(triangles.end(), reconstructed_triangulation.begin(), reconstructed_triangulation.end());
      // or use apprach 2
      //triangles.insert(triangles.end(), cached_triangles.begin(), cached_triangles.end());

      // compare if both do the same (they don't always)
      /*
      CGAL_assertion_code(
            std::sort(triangles.begin(), triangles.end());
            std::sort(cached_triangles.begin(), cached_triangles.end());  )
      CGAL_assertion(triangles == cached_triangles);
      */

      gather_boundary_edges(triangles, boundary_edges_picked);

      return weights_cache->get_best_weight(e_D);
    }

    // construct weight, triangles and bep for this level. Also i, k from e_D
    std::pair<double, double> best_weight = std::make_pair(
                                            std::numeric_limits<double>::max(),
                                            std::numeric_limits<double>::max());
    std::vector<Triangle> best_triangles;
    //boost::container::flat_set<std::pair<int,int> > best_bep;
    //std::set<std::pair<int,int> > best_bep;
    Pair_set best_bep(n);

    int i = e_D.first;
    int k = e_D.second;

    // base case triangle evaluation
    if(domain.b_ids.size() == 3 && !domain.has_islands())
    {
      CGAL_assertion(domain.b_ids[0] == i); // access edge source
      CGAL_assertion(domain.b_ids[2] == k); // access edge target

      int m = domain.b_ids[1]; //third vertex
      const Wpair weight = calculate_weight(i, m, k);
      CGAL_assertion(weight.first >= 0);
      CGAL_assertion(weight.second >= 0);

      // avoid it if it creates non-manifold edges
      if (!boundary_edges_picked.insert ( std::make_pair(k, i) ).second ||
          !boundary_edges_picked.insert ( std::make_pair(i, m) ).second ||
          !boundary_edges_picked.insert ( std::make_pair(m, k) ).second )
      {
        ++count_avoiding_beps;
        return std::make_pair(std::numeric_limits<double>::max(),
                              std::numeric_limits<double>::max());;
      }

      // return the triangle and its weight
      if (weights_cache != NULL)
      {
        weights_cache->put(e_D, m, weight);
        triangulation_cache->put(e_D, {{i, m, k}});
      }
      ++count_triangles;
      triangles.push_back( {{i, m, k}} );
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
      for(std::size_t j = 0; j < domain.islands_list[island_id].size(); ++j)
      {
        // point that is being connected to the boundary
        const int pid = domain.islands_list[island_id][j];

        CGAL_assertion(std::find(domain.islands_list[island_id].begin(), domain.islands_list[island_id].begin(), pid) !=
               domain.islands_list[island_id].end());

        if(skip_facet(i, pid, k))
        {
          count_DT_skips++;
          continue;
        }

        Domain D1, D2;
        // assign the remaining islands the domain(both orientations) that are produced
        D1.add_islands(local_islands);
        D2.add_islands(local_islands);
        // split_domain_case_1 joins the island that pid belongs to
        split_domain_case_1(domain, D1, D2, i, k, domain.islands_list[island_id], j); //to think about getting D2 when !correct orientation

        std::pair<int, int> e_D1 = D1.get_access_edge();
        std::pair<int, int> e_D2 = D2.get_access_edge();
        std::vector<Triangle> triangles_D1, triangles_D2;

        // todo: use bep2 only if !correct_island_orientation
        //boost::container::flat_set< std::pair<int,int> > bep1 = boundary_edges_picked, bep2=bep1;
        //std::set< std::pair<int,int> > bep1 = boundary_edges_picked, bep2=bep1;
        Pair_set bep1 = boundary_edges_picked, bep2=bep1;

        // insert island edges to the bep - todo: refactor this
        for(auto i_it = domain.islands_list[island_id].begin() + 1; i_it != domain.islands_list[island_id].end(); ++i_it)
        {
          bep1.insert(std::make_pair(*i_it, *i_it-1));
          bep2.insert(std::make_pair(*i_it-1, *i_it));
        }
        bep1.insert(std::make_pair(*domain.islands_list[island_id].begin(), *(domain.islands_list[island_id].end()-1)));
        bep2.insert(std::make_pair(*(domain.islands_list[island_id].end()-1), *domain.islands_list[island_id].begin()));

        // bep1 and bep2 here contain initial boundary edges and island edges
        const Wpair w_D1 = process_domain(D1, e_D1, triangles_D1, bep1);

        if(!correct_island_orientation)
        {
          const Wpair w_D2 = process_domain(D2, e_D2, triangles_D2, bep2);

          // is it guaranteed that there will be a valid triangulation after a case I?
          //CGAL_assertion(w_D1.first <= 180);
          //CGAL_assertion(w_D2.first <= 180);

          // choose the best orientation
          if(w_D1 < w_D2)
          {
            // calculate w(t) & add w(t) to w_D1
            const Wpair weight_t = calculate_weight(i, pid, k);
            const Wpair w = add_weights(w_D1, weight_t);
            if(w < best_weight)
            {
              // update the best weight
              best_weight = w;
              // add t to triangles_D2 and return them
              Triangle t = {i, pid, k};
              best_triangles.swap(triangles_D1);
              best_triangles.push_back(t);
              bep1.insert(std::make_pair(k, i));
              bep1.insert(std::make_pair(i, pid));
              bep1.insert(std::make_pair(pid, k));
              best_bep.swap(bep1);
            }
          }
          else
          {
            // calculate w(t) & add w(t) to w_D2
            const Wpair weight_t = calculate_weight(i, pid, k);
            const Wpair w = add_weights(w_D2, weight_t);
            if(w < best_weight)
            {
              // update the best weight
              best_weight = w;
              // add t to triangles_D2 and return them
              Triangle t = {i, pid, k};
              best_triangles.swap(triangles_D2);
              best_triangles.push_back(t);
              bep2.insert(std::make_pair(k, i));
              bep2.insert(std::make_pair(i, pid));
              bep2.insert(std::make_pair(pid, k));
              best_bep.swap(bep2);
           }
          }
        }
        else // if considering only one orientation
        {
          // calculate w(t) & add w(t) to w_D1
          const Wpair weight_t = calculate_weight(i, pid, k);
          const Wpair w = add_weights(w_D1, weight_t);
          if(w < best_weight)
          {
            // update the best weight
            best_weight = w;
            // add t to triangles_D2 and return them
            Triangle t = {i, pid, k};
            best_triangles.swap(triangles_D1);
            best_triangles.push_back(t);
            bep1.insert(std::make_pair(k, i));
            bep1.insert(std::make_pair(i, pid));
            bep1.insert(std::make_pair(pid, k));
            best_bep.swap(bep1);
          }
        }
      } // pid : domains.all_h_ids - case 1 split
    } // end list of islands


    bool table_created_here = false;
    int best_pid = -1;
    if(weights_cache == NULL && !domain.has_islands())
    {
      // create if it has been deleted
      table_created_here = true;
      weights_cache = new WeightsTable(n, std::make_pair(std::numeric_limits<double>::max(),
                                                         std::numeric_limits<double>::max()));

      // just to make sure & compare
      triangulation_cache = new TrianglesTable(n, std::vector<Triangle>());

      // maybe boundary edges picked should be cleared here and reset to
      // the initial bounrady + island edges
    }


    // case II

    // avoid begin and end of the container which is the source and target of the access edge
    for(std::vector<int>::iterator pid_it = domain.b_ids.begin() + 1; pid_it != domain.b_ids.end() - 1; ++pid_it)
    {
      // a triangle that has islands is considered
      // any case split 2 would produce an invalid triangulation because it disconnects boundary and island
      if(domain.b_ids.size() == 3 && domain.has_islands())
        break;

      const int pid = *pid_it; // todo: avoid this maybe

      if(skip_facet(i, pid, k))
      {
        //std::cout<< "Triangle t not in DT!\n" ;
        count_DT_skips++;
        continue;
      }

      // boundary_edges_picked here contain all boundary, island, and the edges that were added at the previous level.
      //boost::container::flat_set<std::pair<int,int> > bep_D1D2 = boundary_edges_picked;
      //std::set<std::pair<int,int> > bep_D1D2 = boundary_edges_picked;
      Pair_set bep_D1D2 = boundary_edges_picked;

      // collect and refuse split to avoid creation of non-manifold edges
      // adds weak edges (triangle t)
      if ( !bep_D1D2.insert ( std::make_pair(k, i) ).second   ||
          !bep_D1D2.insert ( std::make_pair(i, pid) ).second ||
          !bep_D1D2.insert ( std::make_pair(pid, k) ).second )
      {
        //std::cout<< "Triangle t in beps!\n" ;
        ++count_avoiding_beps;
        continue;
      }

      Domain D1, D2;
      // split_domain_case_2 splits the domain to the boundary by creating 2 subdomains
      // D1, D2 have just new boundaries - no island information yet.
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

      // This is evaluated after case I, since case I stops happening
      // only after there are no islands left.
      if(!domain.has_islands())
      {
        w_D12 = add_weights(process_domain(D1, e_D1, triangles_D1, bep_D1D2),
                            process_domain(D2, e_D2, triangles_D2, bep_D1D2) );
      }

      // if domain does have islands, then we don't need to parition if
      // one of D1,D2 is empty: all islands go to the not empty
      else
      {
        if(D1.is_empty())
        {
          // assign all left: local_islands
          D2.add_islands(domain.islands_list);
          w_D12 = process_domain(D2, e_D2, triangles_D2, bep_D1D2);
        }
        else
        {
          if(D2.is_empty())
          {
            // assign all left: local_islands
            D1.add_islands(domain.islands_list);
            w_D12 = process_domain(D1, e_D1, triangles_D1, bep_D1D2);
          }
          // if there are islands in domain, then take all combinations of them on
          // each domain
          else
          {
            CGAL_assertion(!D1.is_empty());
            CGAL_assertion(!D2.is_empty());
            CGAL_assertion(domain.has_islands());

            Phi partition_space;
            do_permutations(domain.islands_list, partition_space);

            for(std::size_t p = 0; p < partition_space.size(); ++p)
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
              //boost::container::flat_set<std::pair<int,int>> local_bep12 = bep_D1D2;
              //std::set<std::pair<int, int> > local_bep12 = bep_D1D2;
              Pair_set local_bep12 = bep_D1D2;

              const Wpair local_w_D1 = process_domain(D1, e_D1, local_triangles_D1,  local_bep12);
              const Wpair local_w_D2 = process_domain(D2, e_D2, local_triangles_D2,  local_bep12);

              if (add_weights(local_w_D1,local_w_D2) < w_D12)
              {
                w_D12=add_weights(local_w_D1, local_w_D2);
                triangles_D1.swap(local_triangles_D1);
                triangles_D2.swap(local_triangles_D2);

                // collect the best bep
                bep_D1D2.swap(local_bep12);
              }
            }
          }
        }
      }

      // calculate w(t)
      const Wpair weight_t = calculate_weight(i, pid, k);
      ++count_triangles;
      // add it to its subdomains
      const Wpair w = add_weights(w_D12, weight_t);

      // at the top level, the best weight here for each pid will be compared with the one produced from cases I
      if(w < best_weight)
      {
        // update the best weight
        best_weight = w;

        // get the best triangles
        Triangle t = {i, pid, k};
        best_triangles.swap(triangles_D1);
        best_triangles.insert(best_triangles.end(), triangles_D2.begin(), triangles_D2.end());
        best_triangles.push_back(t);
        best_pid=pid;

        // get the best bep
        best_bep.swap(bep_D1D2);

      } // w < best weight
    } // case 2 splits

    // store the best_weight
    if (weights_cache!=NULL)
    {
      weights_cache->put(e_D, best_pid, best_weight);
      triangulation_cache->put(e_D, best_triangles);
    }

    // delete W table
    if(table_created_here == true)
    {
      delete weights_cache;
      delete triangulation_cache;
      weights_cache = NULL;
      triangulation_cache = NULL;
    }

    // now copy the triangles from the best triangulation -- swap????
    triangles.insert(triangles.end(), best_triangles.begin(), best_triangles.end());

    // return the bep
    //boundary_edges_picked.insert(best_bep.begin(), best_bep.end());
    //boundary_edges_picked.add_from_triangles(best_triangles);     //use Pair_set
    boundary_edges_picked = best_bep;

    // useful when return for case II splits that have followed case I,
    // so as to compare different case I splits.
    return best_weight;
  }

  const Wpair calculate_weight(const int i, const int m, const int k)
  {
    const Weight& w_t = WC(points, Q, i, m, k, lambda);

    double angle = w_t.w.first;
    double area = w_t.w.second;

    // temp: handle degenerate edges - will be taken care with a new structure for the weight
    // which will produce the numeric limit instead of -1
    if(angle == -1)
      angle = std::numeric_limits<double>::max();
    if(area == -1)
      area = std::numeric_limits<double>::max();

    CGAL_assertion(angle >= 0);
    CGAL_assertion(area >= 0);
    return std::make_pair(angle, area);
  }

  // data
  const PointRange& points;
  const PointRange Q;

  WeightsTable* weights_cache;
  TrianglesTable* triangulation_cache;

  const Domain& domain;
  LambdaTable& lambda;
  const WeightCalculator& WC;

  const int n;

  // for DT3 filtering
  DT3 dt3;
  std::vector<typename DT3::Vertex_handle> dt3_vertices;

  // local dt3 faces and edges
  #ifdef LOCAL_DT_FACES
  std::set<Triangle, Compare_triangles> dt3_faces;
  #endif

  std::set<std::pair<int, int> > dt3_edges;

  bool correct_island_orientation;
  int count_DT_skips;
  int count_triangles;
  int count_table_skips;
  int count_avoiding_beps;

};


} // namespace internal
} // namsepace CGAL





#endif // CGAL_PMP_INTERNAL_HOLE_FILLING_ISLAND_TRIANGULATE_HOLE_POLYLINE_H
