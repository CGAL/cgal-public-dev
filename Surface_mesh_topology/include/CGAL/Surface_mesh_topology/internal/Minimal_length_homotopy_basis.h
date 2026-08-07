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
#include <CGAL/Union_find.h>

#include <vector>
#include <unordered_map>
#include <algorithm>

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
    Base(amesh),
    m_original_mesh(amesh)
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

  // Computes the minimal length homotopy basis rooted at root_vertex: runs
  // the spanning tree (step 1) and the weighted generator selection (steps
  // 2-5), then builds one closed loop per generator by walking the tree
  // from each of its two endpoints back to the root -- the same trace-back
  // logic already used by Base::compute_cycle for a single generator,
  // applied here once per generator instead.
  template <class WeightFunctor>
  std::vector<Path> compute_basis(Original_dart_const_descriptor root_vertex,
                                  const WeightFunctor& wf)
  {
    std::vector<typename WeightFunctor::Weight_t> distance_from_root=
      compute_root_spanning_tree(root_vertex, wf);
    Dart_container generators=compute_generators(distance_from_root, wf);

    std::vector<Path> basis;
    basis.reserve(generators.size());
    for (Dart_descriptor dh : generators)
    {
      Dart_descriptor a=dh, b=this->get_local_map().next(dh);
      int index_a=this->vertex_info(a), index_b=this->vertex_info(b);

      Path loop(m_original_mesh);
      // Trace back the path from `a` to root
      for (int ind=index_a-1; ind!=-1; ind=this->m_trace_index[ind])
      { this->add_to_cycle(this->m_spanning_tree[ind], loop, true); } 
      loop.reverse(); 
      // w root -> a
      this->add_to_cycle(dh, loop, false); // the generator edge itself
      // Trace back the path from `b` to root
      for (int ind=index_b-1; ind!=-1; ind=this->m_trace_index[ind])
      { this->add_to_cycle(this->m_spanning_tree[ind], loop, true); }

      loop.update_is_closed();
      basis.push_back(std::move(loop));
    }

    return basis;
  }

  template <class WeightFunctor=Unit_weight_functor>
  std::vector<Path> compute_basis(Original_dart_const_descriptor root_vertex)
  { return compute_basis(root_vertex, WeightFunctor()); }

protected:
  // Re-exposed for further subclasses (e.g. Minimal_length_homology_basis):
  // since this class privately inherits Shortest_noncontractible_cycle,
  // these are otherwise only visible to this class' own methods, not to
  // anything deriving from it -- private inheritance does not forward
  // access beyond the immediate derived class. Same idea as the
  // already-public `using Base::get_local_map;` above, just protected
  // instead of public since these are implementation details, not part of
  // the public API.
  using Base::m_spanning_tree;
  using Base::m_trace_index;
  using Base::m_copy_to_origin;
  using Base::vertex_info;

  // Assigns an integer id (0-based) to every face of the local map, and
  // fills face_id so that every dart of a given face maps to that face's
  // id. Assumes the input mesh is closed (no boundary, so close<2>() added
  // no artificial faces) -- see the class-level documented limitation.
  // Returns the number of faces.
  int compute_face_ids(std::unordered_map<Dart_descriptor, int>& face_id)
  {
    face_id.clear();
    int nb_faces=0;
    size_type face_seen=this->get_local_map().get_new_mark();
    for (auto it=this->get_local_map().darts().begin(),
              itend=this->get_local_map().darts().end(); it!=itend; ++it)
    {
      if (!this->get_local_map().is_marked(it, face_seen))
      {
        int id=nb_faces++;
        // NOTE: darts_of_cell_basic<2>'s internal marking of visited darts
        // was empirically found to not reliably persist for later
        // is_marked() checks (verified experimentally), so we mark each
        // dart ourselves here -- same pattern as Facewidth.h.
        for (auto itf=this->get_local_map().template darts_of_cell_basic<2>(it, face_seen).begin(),
                  itfend=this->get_local_map().template darts_of_cell_basic<2>(it, face_seen).end();
             itf!=itfend; ++itf)
        {
          this->get_local_map().mark(itf, face_seen);
          face_id[itf]=id;
        }
      }
    }
    this->get_local_map().free_mark(face_seen);
    return nb_faces;
  }

  // One non-tree edge, with the length of the loop it would form if kept
  // as a basis generator (see compute_candidate_edges).
  template <class WeightFunctor>
  struct Candidate
  {
    typename WeightFunctor::Weight_t weight;
    Dart_descriptor dart;

    // Descending order by weight: a plain std::sort(candidates.begin(),
    // candidates.end()) then already gives Kruskal the candidates from
    // longest to shortest candidate loop (see compute_generators).
    bool operator<(const Candidate& other) const
    { return weight>other.weight; }
  };

  // Lists every edge *not* in the spanning tree, together with the weight
  // Kruskal will sort on: length(edge) + distance_from_root[v0] +
  // distance_from_root[v1] -- i.e. the length of the loop this edge would
  // close if it ends up being kept as a generator instead of absorbed into
  // the cotree. distance_from_root must come from a prior call to
  // compute_root_spanning_tree() with the same WeightFunctor.
  template <class WeightFunctor>
  std::vector<Candidate<WeightFunctor>> compute_candidate_edges(
      const std::vector<typename WeightFunctor::Weight_t>& distance_from_root,
      const WeightFunctor& wf)
  {
    size_type in_tree=this->get_local_map().get_new_mark();
    for (Dart_descriptor dh : this->m_spanning_tree)
    {
      if (!this->get_local_map().is_marked(dh, in_tree))
      { this->get_local_map().template mark_cell<1>(dh, in_tree); }
    }

    std::vector<Candidate<WeightFunctor>> candidates;
    for (auto it=this->get_local_map().darts().begin(),
              itend=this->get_local_map().darts().end(); it!=itend; ++it)
    {
      Dart_descriptor opp=this->get_local_map().opposite2(it);
      if (it<opp && !this->get_local_map().is_marked(it, in_tree))
      {
        Dart_descriptor a=it, b=this->get_local_map().next(it);
        int ia=this->vertex_info(a), ib=this->vertex_info(b);
        typename WeightFunctor::Weight_t w=
          distance_from_root[ia]+distance_from_root[ib]
          +wf(Get_original_dart<Base, Copy>::run(this, this->nonhole_dart_of_same_edge(it)));
        candidates.push_back(Candidate<WeightFunctor>{w, it});
      }
    }

    this->get_local_map().free_mark(in_tree);
    return candidates;
  }

  // Kruskal's algorithm for a MAXIMUM spanning tree of the dual graph
  // (faces = nodes, candidate edges = dual edges): processing candidates
  // from longest to shortest candidate loop (Candidate::operator<) means
  // the longest ones get absorbed into the cotree first, leaving the
  // shortest ones as generators -- this is what gives the basis its
  // "minimal length" property.
  // Returns the darts that end up as generators (neither tree nor cotree).
  template <class WeightFunctor>
  Dart_container compute_generators(
      const std::vector<typename WeightFunctor::Weight_t>& distance_from_root,
      const WeightFunctor& wf)
  {
    std::vector<Candidate<WeightFunctor>> candidates=
      compute_candidate_edges(distance_from_root, wf);
    std::sort(candidates.begin(), candidates.end());

    std::unordered_map<Dart_descriptor, int> face_id;
    int nb_faces=compute_face_ids(face_id);

    CGAL::Union_find<int> uf;
    std::vector<typename CGAL::Union_find<int>::handle> handles;
    handles.reserve(nb_faces);
    for (int i=0; i<nb_faces; ++i)
    { handles.push_back(uf.make_set(i)); }

    Dart_container generators;
    for (const auto& c : candidates)
    {
      Dart_descriptor opp=this->get_local_map().opposite2(c.dart);
      int fa=face_id.at(c.dart), fb=face_id.at(opp);
      if (uf.same_set(handles[fa], handles[fb]))
      { generators.push_back(c.dart); }
      else
      { uf.unify_sets(handles[fa], handles[fb]); }
    }

    return generators;
  }

  const Mesh& m_original_mesh;
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_MINIMAL_LENGTH_HOMOTOPY_BASIS_H
