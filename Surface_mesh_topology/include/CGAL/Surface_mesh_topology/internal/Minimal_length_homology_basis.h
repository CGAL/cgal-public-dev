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
#ifndef CGAL_MINIMAL_LENGTH_HOMOLOGY_BASIS_H
#define CGAL_MINIMAL_LENGTH_HOMOLOGY_BASIS_H

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Surface_mesh_topology/internal/Minimal_length_homotopy_basis.h>

#include <unordered_map>
#include <queue>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

// Computes a minimum-weight basis of H_1(M; F_2), following the matroid
// greedy algorithm of Erickson & Whittlesey (2005), section 4. Unlike
// Minimal_length_homotopy_basis (a greedy tree-cotree approximation of a
// homotopy basis), this basis is provably minimum-weight, and its loops are
// guaranteed pairwise linearly independent in homology -- a strictly
// stronger property than the pairwise non-freely-homotopic guarantee the
// sibling class offers.
//
// Privately inherits Minimal_length_homotopy_basis (not
// Shortest_noncontractible_cycle directly), reusing its
// compute_root_spanning_tree (Dijkstra/BFS from an arbitrary root) and
// compute_face_ids (face numbering) unchanged.
//
// Same precondition as Minimal_length_homotopy_basis: the input mesh is
// assumed closed (no boundary).
template <class Mesh_, bool Copy=true>
class Minimal_length_homology_basis : private Minimal_length_homotopy_basis<Mesh_, Copy>
{
public:
  using Base = Minimal_length_homotopy_basis<Mesh_, Copy>;
  using Self = Minimal_length_homology_basis<Mesh_, Copy>;
  using Mesh = Mesh_;

  using Original_dart_const_descriptor = typename Base::Original_dart_const_descriptor;
  using Dart_descriptor       = typename Base::Dart_descriptor;
  using Dart_const_descriptor = typename Base::Dart_const_descriptor;
  using size_type             = typename Base::size_type;
  using Dart_container        = typename Base::Dart_container;
  using Path                  = typename Base::Path; // = Path_on_surface<Mesh>

  using Base::get_local_map;

  Minimal_length_homology_basis(const Mesh& amesh) : Base(amesh) {}

protected:
  // Face id of every dart, and the number of faces -- filled by
  // compute_genus(), reused as-is by Phase 1 (compute_dual_BFS_tree) to
  // number the faces of the dual spanning tree, instead of being
  // recomputed with a second full pass over the darts. Same idea as the
  // inherited m_spanning_tree/m_trace_index: state shared across the
  // different phases of the algorithm lives in a member, not threaded
  // through as function parameters.
  std::unordered_map<Dart_descriptor, int> m_face_id;
  int m_nb_faces;
  int m_genus;

  // Filled by compute_dual_BFS_tree() (Phase 1), reused directly by
  // compute_annotations() (Phase 2) instead of being threaded through as
  // parameters -- same idea as m_face_id/m_nb_faces. m_h_parent[f] is a
  // dart belonging to face f, sitting at the edge shared with f's parent
  // in the dual tree C (indexed like m_face_id). m_X holds the 2*genus
  // generator-candidate darts (neither in T nor C).
  Dart_container m_h_parent;
  Dart_container m_X;

  // Computes the genus of the surface from Euler's formula
  // 2-(V-E+F) = 2*genus (valid for closed orientable input, the class'
  // precondition), and fills m_face_id/m_nb_faces/m_genus. V and E come
  // from the generic BGL free functions num_vertices/num_edges on the
  // *original* mesh (found by ADL, so this works unchanged for every
  // Mesh_ the package supports -- Surface_mesh, Polyhedron_3,
  // LCC_for_CMap_2, LCC_for_GMap_2). F comes from the inherited
  // compute_face_ids instead of the BGL equivalent, since it numbers the
  // faces as a byproduct of counting them, and that numbering is needed
  // again in Phase 1.
  int compute_genus()
  {
    m_nb_faces=this->compute_face_ids(m_face_id);

    int nb_vertices=static_cast<int>(num_vertices(this->m_original_mesh));
    int nb_edges=static_cast<int>(num_edges(this->m_original_mesh));

    int euler_characteristic=nb_vertices-nb_edges+m_nb_faces;
    m_genus=(2-euler_characteristic)/2;

    return m_genus;
  }

  // Builds a spanning tree C of the dual graph (faces as nodes), restricted
  // to darts not already used by the primal tree T (this->m_spanning_tree,
  // from a prior compute_root_spanning_tree call). Faces are indexed via
  // m_face_id (filled once by compute_genus, reused here directly instead
  // of being recomputed). The root face is chosen arbitrarily as face id 0
  // -- the face containing the very first dart of the local map, since
  // compute_face_ids assigns id 0 to whichever face it discovers first
  // while scanning darts from the beginning. For each non-root face f,
  // h_parent[f] is a dart belonging to f itself, sitting at the edge
  // shared with f's parent -- both the entry point used to walk f's own
  // contour here, and (via opposite2) the parent face's own dart of that
  // same edge, needed in Phase 2. Darts that are neither in T nor end up
  // in C become X, the 2*genus generator-candidate edges (a byproduct of
  // the same traversal, same idea as
  // compute_candidate_edges/compute_generators in the sibling class).
  void compute_dual_BFS_tree()
  {
    size_type in_tree=this->get_local_map().get_new_mark();
    for (Dart_descriptor dh : this->m_spanning_tree)
    { if (!this->get_local_map().is_marked(dh, in_tree)) { this->get_local_map().template mark_cell<1>(dh, in_tree); } }

    size_type in_cotree=this->get_local_map().get_new_mark();
    size_type face_visited=this->get_local_map().get_new_mark();

    m_h_parent.assign(m_nb_faces, Dart_descriptor());
    m_X.clear();

    Dart_descriptor root_face_dart=this->get_local_map().darts().begin();
    int root_id=0;
    std::queue<int> q;
    q.push(root_id);
    this->get_local_map().template mark_cell<2>(root_face_dart, face_visited);

    while (!q.empty())
    {
      int f_index=q.front();
      q.pop();

      Dart_descriptor start;
      if (f_index==root_id) { start=root_face_dart; }
      else { start=m_h_parent[f_index]; }

      for (auto it=this->get_local_map().template darts_of_cell_basic<2>(start).begin(),
                itend=this->get_local_map().template darts_of_cell_basic<2>(start).end();
           it!=itend; ++it)
      {
        if (!this->get_local_map().is_marked(it, in_tree) &&
            !this->get_local_map().is_marked(it, in_cotree))
        {
          Dart_descriptor opp=this->get_local_map().opposite2(it);
          int g_index=m_face_id.at(opp);
          if (!this->get_local_map().is_marked(opp, face_visited))
          {
            m_h_parent[g_index]=opp;
            this->get_local_map().template mark_cell<1>(it, in_cotree);
            this->get_local_map().template mark_cell<2>(opp, face_visited);
            q.push(g_index);
          }
          else if (it<opp)
          { m_X.push_back(it); }
        }
      }
    }

    this->get_local_map().free_mark(face_visited);
    this->get_local_map().free_mark(in_cotree);
    this->get_local_map().free_mark(in_tree);
  }

  // Filled by compute_annotations() (Phase 2): the F_2^{2*genus}
  // annotation of every dart, the same value for both darts of an edge.
  std::unordered_map<Dart_descriptor, boost::dynamic_bitset<>> m_annotation;

  // Computes the F_2^{2*genus} annotation of every edge (Phase 2): a
  // bitvector of size 2*genus per dart, such that XOR-ing the annotations
  // along any cycle gives its class in this basis of H_1(M; F_2). T-edges
  // get the zero vector; X-edges get the 2*genus unit vectors of the
  // basis, one each, in the order they were collected in m_X. C-edges are
  // solved for by requiring the XOR around every face's contour to be
  // zero: propagated recursively over the dual tree C, children before
  // their parent (built from m_h_parent), so a face's parent edge is
  // always the only unknown annotation left on its contour when that face
  // is processed.
  void compute_annotations()
  {
    m_annotation.clear();

    boost::dynamic_bitset<> zero(2*m_genus);
    for (Dart_descriptor dh : this->m_spanning_tree)
    {
      m_annotation[dh]=zero;
      m_annotation[this->get_local_map().opposite2(dh)]=zero;
    }

    for (std::size_t i=0; i<m_X.size(); ++i)
    {
      boost::dynamic_bitset<> unit(2*m_genus);
      unit[i]=true;
      m_annotation[m_X[i]]=unit;
      m_annotation[this->get_local_map().opposite2(m_X[i])]=unit;
    }

    std::vector<std::vector<int>> children(m_nb_faces);
    for (int f=1; f<m_nb_faces; ++f)
    {
      int parent=m_face_id.at(this->get_local_map().opposite2(m_h_parent[f]));
      children[parent].push_back(f);
    }

    // First pass: discover every face starting from the root via an
    // explicit stack (not the call stack). This always discovers a face
    // strictly before its own children -- pushing a node's children onto
    // the stack only queues them for later, it does not wait for them to
    // finish first. discovery_order ends up in that same "ancestor before
    // descendants" order -- in particular, discovery_order[0] is always
    // the root (face 0), the very first thing pushed and popped.
    std::vector<int> discovery_order;
    std::vector<int> stack;
    stack.push_back(0);
    while (!stack.empty())
    {
      int f_index=stack.back();
      stack.pop_back();
      discovery_order.push_back(f_index);
      for (int child : children[f_index]) { stack.push_back(child); }
    }

    // Second pass: process discovery_order in reverse, which flips it to
    // "descendants before ancestors" -- exactly the order propagation
    // needs, since a face's own contour includes the edges to its
    // children, whose annotations must already be known. Stops at i>1
    // (not i>0) to skip discovery_order[0], the root, which has no parent
    // edge to compute -- no need for an if(f_index!=0) check inside.
    for (std::size_t i=discovery_order.size(); i>1; --i)
    {
      int f_index=discovery_order[i-1];
      Dart_descriptor entry=m_h_parent[f_index];
      boost::dynamic_bitset<> value(2*m_genus);
      for (Dart_descriptor it=this->get_local_map().next(entry); it!=entry;
           it=this->get_local_map().next(it))
      { value^=m_annotation.at(it); }
      m_annotation[entry]=value;
      m_annotation[this->get_local_map().opposite2(entry)]=value;
    }
  }

  // One non-tree (relative to T_y) edge, with the length and F_2^{2*genus}
  // vector of the candidate cycle it would close if kept as a basis
  // generator rooted at y (see compute_candidates, not yet written).
  // Sorted ascending by length, ready for Phase 4's greedy (unlike the
  // sibling's Candidate<WeightFunctor>, which sorts descending for its
  // Kruskal maximum-cotree).
  template <class WeightFunctor>
  struct Candidate
  {
    typename WeightFunctor::Weight_t weight;
    boost::dynamic_bitset<> vector;
    Dart_descriptor root; // y
    Dart_descriptor edge; // e

    bool operator<(const Candidate& other) const
    { return weight<other.weight; }
  };

  // The array index of v_index's parent in the current T_y. m_trace_index
  // already stores parent_vertex_index-1 (that's exactly why compute_basis
  // in the sibling class treats -1 as "the root"): adding 1 back always
  // recovers the correct 0-based vertex index, root included (-1+1=0), no
  // special case needed.
  int parent_index(int v_index)
  { return this->m_trace_index[v_index-1]+1; }

  // Walk-to-root-with-memoization: returns A_y[v] (XOR of edge annotations
  // along T_y from y to the vertex with vertex_info()==v_index). value[idx]
  // holds A_y[idx] for every vertex already computed so far, computed[idx]
  // says whether that entry is filled in yet; both are meant to be created
  // once per root y and passed back into this function for every candidate
  // edge of that y (see compute_candidates, not yet written), so a vertex
  // already computed for one candidate is an instant lookup for the next.
  // Computes v_index's own entry, if not already known, by following
  // parent_index up from v_index until an already-known vertex is found
  // (path ends up ordered from v_index towards that ancestor), then
  // filling in every visited vertex in the opposite order (ancestor-side
  // first), since a vertex's own entry needs its parent's to already be
  // known. i counts down from path.size() to 1 (not path.size()-1 to 0)
  // so it stays a valid size_t the whole time -- decrementing an unsigned
  // 0 would wrap around instead of going negative.
  boost::dynamic_bitset<> get_annotation_to_root(
      int v_index, std::vector<bool>& computed, std::vector<boost::dynamic_bitset<>>& value)
  {
    std::vector<int> path;
    int idx=v_index;
    while (!computed[idx])
    {
      path.push_back(idx);
      idx=parent_index(idx);
    }

    for (std::size_t i=path.size(); i>0; --i)
    {
      int cur=path[i-1];
      value[cur]=value[parent_index(cur)]^this->m_annotation.at(this->m_spanning_tree[cur-1]);
      computed[cur]=true;
    }

    return value[v_index];
  }

  // Phase 3: for every vertex y of the local map, reruns
  // compute_root_spanning_tree(y, wf) -- a fresh Dijkstra/BFS tree T_y,
  // unrelated to the primal tree T built once in Phase 1 -- and lists
  // every edge outside T_y whose candidate-loop homology class (via
  // get_annotation_to_root) is non-zero. Iterates over every dart of the
  // local map, not attributes<1>(), because this map type has no edge
  // attribute defined at all (checked: Items_for_shortest_noncontractible_
  // cycle's Dart_wrapper only declares a vertex attribute,
  // attributes<1>() does not even compile) -- same reason the sibling's
  // compute_candidate_edges does the same thing, canonicalizing each edge
  // via it<opp so it's only considered once.
  template <class WeightFunctor>
  std::vector<Candidate<WeightFunctor>> compute_candidates(const WeightFunctor& wf)
  {
    std::vector<Candidate<WeightFunctor>> candidates;

    for (auto ity=this->get_local_map().template attributes<0>().begin(),
              itYend=this->get_local_map().template attributes<0>().end();
         ity!=itYend; ++ity)
    {
      Dart_descriptor y=this->get_local_map().template dart_of_attribute<0>(ity);
      Original_dart_const_descriptor original_y=this->m_copy_to_origin.at(y);
      std::vector<typename WeightFunctor::Weight_t> distance_from_root=
        this->compute_root_spanning_tree(original_y, wf);

      std::vector<bool> computed(distance_from_root.size(), false);
      std::vector<boost::dynamic_bitset<>> annotation_from_y(
          distance_from_root.size(), boost::dynamic_bitset<>(2*m_genus));
      computed[0]=true;

      size_type in_tree_y=this->get_local_map().get_new_mark();
      for (Dart_descriptor dh : this->m_spanning_tree)
      { if (!this->get_local_map().is_marked(dh, in_tree_y)) { this->get_local_map().template mark_cell<1>(dh, in_tree_y); } }

      for (auto it=this->get_local_map().darts().begin(),
                itend=this->get_local_map().darts().end(); it!=itend; ++it)
      {
        Dart_descriptor opp=this->get_local_map().opposite2(it);
        if (it<opp && !this->get_local_map().is_marked(it, in_tree_y))
        {
          Dart_descriptor a=it, b=this->get_local_map().next(it);
          int ia=this->vertex_info(a), ib=this->vertex_info(b);

          boost::dynamic_bitset<> Aa=get_annotation_to_root(ia, computed, annotation_from_y);
          boost::dynamic_bitset<> Ab=get_annotation_to_root(ib, computed, annotation_from_y);
          boost::dynamic_bitset<> vec=Aa^this->m_annotation.at(it)^Ab;

          if (vec.any())
          {
            typename WeightFunctor::Weight_t len=
              distance_from_root[ia]+distance_from_root[ib]+wf(this->m_copy_to_origin.at(it));
            candidates.push_back(Candidate<WeightFunctor>{len, vec, y, it});
          }
        }
      }

      this->get_local_map().free_mark(in_tree_y);
    }

    return candidates;
  }

  // Phase 4: incremental Gauss elimination over GF(2) ("XOR basis"/"linear
  // basis" -- see the design discussion, no reusable implementation found
  // in CGAL or Boost). Sorts candidates by ascending length, then keeps a
  // persistent pivot[] array (row-echelon, not reduced -- see the design
  // discussion) reduced against every candidate in turn: pivot[i] starts
  // "empty" (default-constructed dynamic_bitset, size 0 -- distinct from
  // a same-size all-zero bitset, which is what an *assigned* pivot could
  // never be, since it's only ever assigned when its own bit i is set).
  // Full bit-scan for now (find_first()/find_next() optimization
  // deferred, see the design discussion). Stops once 2*genus candidates
  // are accepted.
  template <class WeightFunctor>
  std::vector<Candidate<WeightFunctor>> matroid_greedy(std::vector<Candidate<WeightFunctor>> candidates)
  {
    std::sort(candidates.begin(), candidates.end());

    std::vector<boost::dynamic_bitset<>> pivot(2*m_genus);
    std::vector<Candidate<WeightFunctor>> basis;

    for (std::size_t c=0; c<candidates.size() && basis.size()<pivot.size(); ++c)
    {
      boost::dynamic_bitset<> vec=candidates[c].vector;
      bool accepted=false;
      for (std::size_t i=0; i<vec.size(); ++i)
      {
        if (vec[i] && !accepted)
        {
          if (pivot[i].empty())
          {
            pivot[i]=vec;
            accepted=true;
          }
          else
          {
            vec^=pivot[i];
          }
        }
      }
      if (accepted) { basis.push_back(candidates[c]); }
    }

    return basis;
  }

  // Wires Phases 0-2 together: genus/face numbering (compute_genus), the
  // primal tree T from an arbitrary root (compute_root_spanning_tree,
  // unweighted -- Phases 1-2 are purely structural, so any spanning tree
  // works), the dual tree C and generator set X (compute_dual_BFS_tree),
  // then the F_2^{2*genus} edge annotations (compute_annotations). A
  // partial, early version of what the eventual public compute_basis(wf)
  // will do; used here so Phases 0-2 can be exercised end-to-end on an
  // actual instance.
  int compute_tree_cotree()
  {
    int genus=this->compute_genus();
    auto root=*(halfedges(this->m_original_mesh).begin());
    this->compute_root_spanning_tree(root);
    this->compute_dual_BFS_tree();
    this->compute_annotations();
    return genus;
  }
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_MINIMAL_LENGTH_HOMOLOGY_BASIS_H
