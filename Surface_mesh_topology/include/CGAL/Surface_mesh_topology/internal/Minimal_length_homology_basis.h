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
#ifndef CGAL_MINIMAL_HOMOLOGY_BASIS_H
#define CGAL_MINIMAL_HOMOLOGY_BASIS_H

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Surface_mesh_topology/internal/Minimal_length_homotopy_basis.h>

#include <unordered_map>
#include <queue>

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
class Minimal_homology_basis : private Minimal_length_homotopy_basis<Mesh_, Copy>
{
public:
  using Base = Minimal_length_homotopy_basis<Mesh_, Copy>;
  using Self = Minimal_homology_basis<Mesh_, Copy>;
  using Mesh = Mesh_;

  using Original_dart_const_descriptor = typename Base::Original_dart_const_descriptor;
  using Dart_descriptor       = typename Base::Dart_descriptor;
  using Dart_const_descriptor = typename Base::Dart_const_descriptor;
  using size_type             = typename Base::size_type;
  using Dart_container        = typename Base::Dart_container;
  using Path                  = typename Base::Path; // = Path_on_surface<Mesh>

  using Base::get_local_map;

  Minimal_homology_basis(const Mesh& amesh) : Base(amesh) {}

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

  // Computes the genus of the surface from Euler's formula
  // 2-(V-E+F) = 2*genus (valid for closed orientable input, the class'
  // precondition), and fills m_face_id/m_nb_faces. V and E come from the
  // generic BGL free functions num_vertices/num_edges on the *original*
  // mesh (found by ADL, so this works unchanged for every Mesh_ the
  // package supports -- Surface_mesh, Polyhedron_3, LCC_for_CMap_2,
  // LCC_for_GMap_2). F comes from the inherited compute_face_ids instead
  // of the BGL equivalent, since it numbers the faces as a byproduct of
  // counting them, and that numbering is needed again in Phase 1.
  int compute_genus()
  {
    m_nb_faces=this->compute_face_ids(m_face_id);

    int nb_vertices=static_cast<int>(num_vertices(this->m_original_mesh));
    int nb_edges=static_cast<int>(num_edges(this->m_original_mesh));

    int euler_characteristic=nb_vertices-nb_edges+m_nb_faces;

    return (2-euler_characteristic)/2;
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
  void compute_dual_BFS_tree(Dart_container& h_parent, Dart_container& X)
  {
    size_type in_tree=this->get_local_map().get_new_mark();
    for (Dart_descriptor dh : this->m_spanning_tree)
    { if (!this->get_local_map().is_marked(dh, in_tree)) { this->get_local_map().template mark_cell<1>(dh, in_tree); } }

    size_type in_cotree=this->get_local_map().get_new_mark();
    size_type face_visited=this->get_local_map().get_new_mark();

    h_parent.assign(m_nb_faces, Dart_descriptor());
    X.clear();

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
      else { start=h_parent[f_index]; }

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
            h_parent[g_index]=opp;
            this->get_local_map().template mark_cell<1>(it, in_cotree);
            this->get_local_map().template mark_cell<2>(opp, face_visited);
            q.push(g_index);
          }
          else if (it<opp)
          { X.push_back(it); }
        }
      }
    }

    this->get_local_map().free_mark(face_visited);
    this->get_local_map().free_mark(in_cotree);
    this->get_local_map().free_mark(in_tree);
  }

  // Wires Phase 0 and Phase 1 together: computes the genus (compute_genus,
  // also filling m_face_id/m_nb_faces), builds the primal tree T from an
  // arbitrary root (compute_root_spanning_tree, unweighted -- Phase 1 is
  // purely structural, so any spanning tree works), then the dual tree C
  // and the generator-candidate set X (compute_dual_BFS_tree). A partial,
  // early version of what the eventual public compute_basis(wf) will do;
  // used here so Phase 1 can be exercised end-to-end on an actual
  // instance.
  int compute_tree_cotree(Dart_container& h_parent, Dart_container& X)
  {
    int genus=this->compute_genus();
    auto root=*(halfedges(this->m_original_mesh).begin());
    this->compute_root_spanning_tree(root);
    this->compute_dual_BFS_tree(h_parent, X);
    return genus;
  }
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_MINIMAL_HOMOLOGY_BASIS_H
