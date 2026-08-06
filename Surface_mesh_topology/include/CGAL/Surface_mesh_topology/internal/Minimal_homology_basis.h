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
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_MINIMAL_HOMOLOGY_BASIS_H
