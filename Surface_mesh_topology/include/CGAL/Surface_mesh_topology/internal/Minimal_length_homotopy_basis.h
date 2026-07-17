// Copyright (c) 2026 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Raphael Costil
//
#ifndef CGAL_MINIMAL_LENGTH_HOMOTOPY_BASIS_H
#define CGAL_MINIMAL_LENGTH_HOMOTOPY_BASIS_H

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Surface_mesh_topology/internal/Shortest_noncontractible_cycle.h>
#include <CGAL/Surface_mesh_topology/internal/Edge_weight_functor.h>

#include <vector>

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

// Computes a minimal length homotopy basis of a surface mesh, following the
// greedy tree-cotree algorithm of Eppstein (2003) / Erickson & Whittlesey
// (2005), as translated from cinolib's homotopy_basis.cpp.
//
// This class privately inherits Shortest_noncontractible_cycle, which
// already builds the closed local map and the root-based Dijkstra/BFS
// spanning tree that both algorithms need; it only adds the pieces that are
// specific to computing a whole basis (weighted tree-cotree selection of the
// generator edges, then one loop per generator) rather than a single
// shortest cycle. Private inheritance is used because a homotopy basis is
// not "a" shortest noncontractible cycle from the caller's point of view:
// Base's own public API (compute_cycle, compute_edge_width, ...) must not
// leak through.
template <class Mesh_, bool Copy=true>
class Minimal_length_homotopy_basis : private Shortest_noncontractible_cycle<Mesh_, Copy>
{
public:
  using Base = Shortest_noncontractible_cycle<Mesh_, Copy>;
  using Self = Minimal_length_homotopy_basis<Mesh_, Copy>;
  using Mesh = Mesh_;

  using Original_dart_const_descriptor = typename Base::Original_dart_const_descriptor;
  using Dart_descriptor       = typename Base::Dart_descriptor;
  using Dart_const_descriptor = typename Base::Dart_const_descriptor;
  using size_type             = typename Base::size_type;
  using Dart_container        = typename Base::Dart_container;
  using Path                  = typename Base::Path; // = Path_on_surface<Mesh>

  using Base::get_local_map;

  Minimal_length_homotopy_basis(const Mesh& amesh) :
    Base(amesh)
  {}

  // Step 1 of the algorithm: build the shortest-path spanning tree of the
  // local map, rooted at root_vertex (a descriptor of the *original* mesh
  // passed in by the caller). Uses Dijkstra if wf is weighted, plain BFS if
  // wf is Unit_weight_functor (dispatch happens inside compute_spanning_tree,
  // reused unchanged from Base).
  //
  // On return: this->get_local_map() darts have their vertex_info() set to
  // a 0-based visit index; Base::m_spanning_tree[i] and Base::m_trace_index[i]
  // describe, for the vertex with index i+1, which dart reaches it from its
  // tree-parent and what the parent's index is (-1 meaning the parent is the
  // root) -- this is exactly what is needed to walk any vertex back to the
  // root one step at a time. The returned vector holds, at index i, the tree
  // distance from the root to the vertex with that index.
  //
  // The weighted selection of which non-tree edges become basis generators,
  // and the assembly of the basis loops themselves, are implemented
  // separately (not yet done at this stage).
  template <class WeightFunctor>
  std::vector<typename WeightFunctor::Weight_t>
  compute_root_spanning_tree(Original_dart_const_descriptor root_vertex,
                             const WeightFunctor& wf)
  {
    Dart_descriptor root=this->m_origin_to_copy[root_vertex];

    std::vector<typename WeightFunctor::Weight_t> distance_from_root;
    this->m_spanning_tree.clear();
    this->m_trace_index.clear();
    this->initialize_vertex_info();
    this->compute_spanning_tree(root, this->m_spanning_tree, distance_from_root,
                                this->m_trace_index, wf);

    return distance_from_root;
  }

  template <class WeightFunctor=Unit_weight_functor>
  std::vector<typename WeightFunctor::Weight_t>
  compute_root_spanning_tree(Original_dart_const_descriptor root_vertex)
  { return compute_root_spanning_tree(root_vertex, WeightFunctor()); }

protected:
  // Disjoint-set (union-find) over face indices, used by the Kruskal-style
  // weighted cotree selection: two faces are in the same set iff they are
  // already connected through edges already assigned to the cotree.
  // No path compression / union by rank for now -- kept as simple as possible.
  struct Union_find
  {
    explicit Union_find(std::size_t n) : parent(n)
    { for (std::size_t i=0; i<n; ++i) parent[i]=static_cast<int>(i); }

    int find(int x)
    {
      while (parent[x]!=x) x=parent[x];
      return x;
    }

    // Returns true iff a and b were in different sets (and have now been merged).
    bool union_sets(int a, int b)
    {
      a=find(a); b=find(b);
      if (a==b) return false;
      parent[b]=a;
      return true;
    }

    std::vector<int> parent;
  };
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_MINIMAL_LENGTH_HOMOTOPY_BASIS_H
