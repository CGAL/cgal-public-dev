// Copyright (c) 2019-2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H

//includes from cgal
#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/internal/Repair/helper.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
//#include <CGAL/Polygon_mesh_processing/repair_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/utility.h>

#include <CGAL/Surface_mesh/Surface_mesh.h>

//includes from in-meshing
#include <core/util.h>
#include <core/output.h>
#include <core/mesh_completion.h>
#include <core/adjust_geometry.h>
#include <core/dualcontouring/connectivity.h>


//includes from stl
#include <iterator>
#include <fstream>
#include <vector>

namespace CGAL {
namespace In_meshing_tools_for_repair_manifoldness {

  template <typename PolygonMesh>
  using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
  template <typename PolygonMesh>
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  template <typename PolygonMesh>
  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;
  template <typename PolygonMesh>
  using face_range = std::set<face_descriptor<PolygonMesh>>;
  template <typename PolygonMesh>
  using halfedge_pair = std::pair<halfedge_descriptor<PolygonMesh>, halfedge_descriptor<PolygonMesh>>;

  template <typename WorkZone, typename PolygonMesh>
  void erase_wz_faces(const WorkZone& wz_1,
                      const WorkZone& wz_2,
                      PolygonMesh& pmesh)
  {
    for(auto& f : wz_1.faces)
    {
      Euler::remove_face(halfedge(f, pmesh), pmesh);
    }
    for(auto& f : wz_2.faces)
    {
      Euler::remove_face(halfedge(f, pmesh), pmesh);
    }
  }

  template <typename WorkZone, typename PolygonMesh>
  face_range<PolygonMesh> make_rings(const WorkZone& wz_1,
                                     const WorkZone& wz_2,
                                     const PolygonMesh& pmesh,
                                     unsigned nb_expand = 3)
  {
    typedef CGAL::dynamic_face_property_t<bool>          Face_bool_tag;
    typedef CGAL::Face_filtered_graph<PolygonMesh>       Filtered_graph;

    face_range<PolygonMesh> rings;
    face_range<PolygonMesh> faces;
    std::merge(wz_1.faces.begin(), wz_1.faces.end(), wz_2.faces.begin(), wz_2.faces.end(), std::inserter(faces, faces.begin()));

    auto selected = get(Face_bool_tag(), pmesh);
    for(auto& f : faces)
    {
      put(selected, f, true);
    }
    CGAL::expand_face_selection(faces, pmesh, nb_expand, selected, std::inserter(rings, rings.end()));

    std::ofstream os;
    Filtered_graph ffg_second_rings(pmesh, rings);
    PolygonMesh second_rings_sm;
    CGAL::copy_face_graph(ffg_second_rings, second_rings_sm);
    os = std::ofstream("/home/felix/Bureau/Geo_Facto/PSR/tests-repair/dumps/make_first_ring.off");
    CGAL::write_off(os, second_rings_sm);
    os.close();

    return rings;
  }

  template <typename GeomTraits>
  void make_guides(point_mesh& pm,
                   std::vector<typename GeomTraits::Point_3>& guides,
                   std::vector<typename GeomTraits::Vector_3>& normals)
  {
    unsigned n = guides.size();
    for(unsigned i = 0; i < n; ++i)
    {
      Eigen::Vector3f p(guides[i].x(), guides[i].y(), guides[i].z());
      Eigen::Vector3f n(normals[i].x(), normals[i].y(), normals[i].z());
      pm.vertices.emplace_back(p, n);
    }
  }

  template <typename WorkZone, typename PolygonMesh, typename VPM, typename GeomTraits>
  point_mesh make_point_mesh_for_in_meshing(const WorkZone& wz_1,
                                            const WorkZone& wz_2,
                                            PolygonMesh& pmesh,
                                            const VPM vpm,
                                            const GeomTraits& gt,
                                            std::vector<unsigned>& hole_indices_1,
                                            std::vector<unsigned>& hole_indices_2,
                                            std::vector<typename GeomTraits::Point_3>& guides,
                                            std::vector<typename GeomTraits::Vector_3>& normals,
                                            unsigned nb_expand = 3)
  {
//    face_range<PolygonMesh> rings = make_rings(wz_1, wz_2, pmesh, nb_expand);

    face_range<PolygonMesh> rings;
    for(auto& f : wz_1.ring)
      rings.insert(f);
    for(auto& f : wz_2.ring)
      rings.insert(f);

    std::map<vertex_descriptor<PolygonMesh>, unsigned> pmesh_vertices_to_psoup_indices;

    unsigned i = 0;
    std::vector<oriented_point> points;
    std::set<vertex_descriptor<PolygonMesh>> treated;
    for(const auto& f : rings)
    {
      for(const auto& v : vertices_around_face(halfedge(f, pmesh), pmesh))
      {
        if(treated.count(v) > 0)
        {
          continue;
        }

        typename GeomTraits::Point_3 p = get(vpm, v);
        auto n = CGAL::Polygon_mesh_processing::compute_vertex_normal(v, pmesh);
        Eigen::Vector3f p_(p.x(), p.y(), p.z());
        Eigen::Vector3f n_(n.x(), n.y(), n.z());

        points.emplace_back(p_, n_);
        pmesh_vertices_to_psoup_indices[v] = i;
        ++i;
        treated.insert(v);
      }
    }

    std::vector<unsigned> faces;
    for(const auto& f : rings)
    {
      for(const auto& v : vertices_around_face(halfedge(f, pmesh), pmesh))
      {
        faces.push_back(pmesh_vertices_to_psoup_indices[v]);
      }
    }

    for(auto& h : wz_1.border)
    {
      vertex_descriptor<PolygonMesh> v = source(h, pmesh);
      hole_indices_1.push_back(pmesh_vertices_to_psoup_indices[v]);
    }
    for(auto& h : wz_2.border)
    {
      vertex_descriptor<PolygonMesh> v = source(h, pmesh);
      hole_indices_2.push_back(pmesh_vertices_to_psoup_indices[v]);
    }

    point_mesh pm = point_mesh(points, faces);
    make_guides<GeomTraits>(pm, guides, normals);

    return pm;
  }

  template <typename PolygonMesh, typename GeomTraits>
  PolygonMesh make_polygon_mesh_from_completed_mesh(const output& out, const completed_mesh& out_mesh, const GeomTraits& gt)
  {
    PolygonMesh poly_mesh;

    std::vector<typename GeomTraits::Point_3> points;
    for(auto& v : out_mesh.vertices)
    {
      double x = from_hpos<float>(out.back_transform * to_hpos(v.position))[0];
      double y = from_hpos<float>(out.back_transform * to_hpos(v.position))[1];
      double z = from_hpos<float>(out.back_transform * to_hpos(v.position))[2];
      points.emplace_back(x, y, z);
    }

    std::vector<std::vector<std::size_t>> polygons;
    for(unsigned i = 0; i < out_mesh.indices.size(); ++i)
    {
      if(i % 3 == 0)
      {
        polygons.emplace_back();
      }
      polygons.back().push_back(out_mesh.indices[i]);
    }

    CGAL::Polygon_mesh_processing::repair_polygon_soup(points, polygons);

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, poly_mesh);

    return poly_mesh;
  }

  template <typename GeomTraits>
  Surface_mesh<typename GeomTraits::Point_3> completed_mesh_to_surface_mesh(const completed_mesh& cm)
  {
    Surface_mesh<typename GeomTraits::Point_3> sm;

    std::vector<typename GeomTraits::Point_3> points;
    for(auto& v : cm.vertices)
    {
      double x = v.position[0];
      double y = v.position[1];
      double z = v.position[2];
      points.emplace_back(x, y, z);
    }

    std::vector<std::vector<std::size_t>> polygons;
    for(unsigned i = 0; i < cm.indices.size(); ++i)
    {
      if(i % 3 == 0)
      {
        polygons.emplace_back();
      }
      polygons.back().push_back(cm.indices[i]);
    }

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, sm);

    return sm;
  }

  template <typename GeomTraits>
  void surface_mesh_to_completed_mesh(const Surface_mesh<typename GeomTraits::Point_3>& sm, completed_mesh& mesh)
  {
    typedef typename GeomTraits::Point_3 Point_3;

    completed_mesh truc;

    std::vector<Point_3> points;
    std::vector<std::vector<unsigned>> polygons;

    CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(sm, points, polygons);

    mesh.vertices.clear();
    unsigned i = 0;
    for(auto& v : sm.vertices())
    {
      auto p = points[i];
      auto n = CGAL::Polygon_mesh_processing::compute_vertex_normal(v, sm);
      Eigen::Vector3f p_ (p.x(), p.y(), p.z());
      Eigen::Vector3f n_ (n.x(), n.y(), n.z());
      mesh.vertices.emplace_back(p_, n_);
      ++i;
    }

    mesh.indices.clear();
    for(auto& poly : polygons)
    {
      for(auto& j : poly)
      {
        mesh.indices.push_back(j);
      }
    }
  }

  template <typename GeomTraits>
  std::vector<typename Surface_mesh<typename GeomTraits::Point_3>::face_index>
      make_faces_for_remeshing(const Surface_mesh<typename GeomTraits::Point_3>& sm,
                               const completed_mesh& mesh)
  {
    typedef Surface_mesh<typename GeomTraits::Point_3> Surface_mesh;
    typedef std::vector<typename Surface_mesh::face_index> Face_range;

    Face_range faces;

    unsigned i = 0;
    for(auto& f : sm.faces())
    {
      if(i < mesh.num_base_triangles)
      {
        ++i;
        continue;
      }
      faces.push_back(f);
    }

    return faces;
  }

  template<typename GeomTraits>
  void remesh_with_CGAL(completed_mesh& mesh, float target_length)
  {
    typedef Surface_mesh<typename GeomTraits::Point_3>                         Surface_mesh;
    typedef CGAL::dynamic_edge_property_t<bool>                                Edge_bool_tag;
    typedef typename boost::property_map<Surface_mesh, Edge_bool_tag>::type    Mark_map;

    Surface_mesh sm = completed_mesh_to_surface_mesh<GeomTraits>(mesh);

    std::vector<typename Surface_mesh::face_index> faces = make_faces_for_remeshing<GeomTraits>(sm, mesh);

    Mark_map constrained_edges = get(Edge_bool_tag(), sm);
    unsigned i = 0;
    for(auto& f : sm.faces())
    {
      if(i < mesh.num_base_triangles + 10)
      {
        for(auto& h : halfedges_around_face(halfedge(f, sm), sm))
        {
          put(constrained_edges, edge(h, sm), true);
        }
        ++i;
      }
      else
      {
        break;
      }
    }

    namespace PMP = CGAL::Polygon_mesh_processing;
    CGAL::Polygon_mesh_processing::isotropic_remeshing(faces, target_length, sm,
                                                       PMP::parameters::edge_is_constrained_map(constrained_edges)
                                                       .protect_constraints(true));


    surface_mesh_to_completed_mesh<GeomTraits>(sm, mesh);
  }

  template <typename PolygonMesh, typename GeomTraits>
  PolygonMesh extract_reconstruction(const output& out, const completed_mesh& out_mesh,
                                     const GeomTraits& gt,
                                     const std::vector<unsigned>& hole_indices_1,
                                     const std::vector<unsigned>& hole_indices_2,
                                     std::vector<std::vector<halfedge_descriptor<PolygonMesh>>>& reconstruction_holes)
  {
    typedef CGAL::dynamic_face_property_t<bool>                                   Face_bool_tag;
    typedef typename boost::property_map<PolygonMesh, Face_bool_tag>::type        Mark_map;

    PolygonMesh reconstruction = make_polygon_mesh_from_completed_mesh<PolygonMesh, GeomTraits>(out, out_mesh, gt);
    std::ofstream os;

    std::string recon_file = "/home/felix/Bureau/Geo_Facto/PSR/tests-repair/dumps/inmsh-reconstruction.off";
    os = std::ofstream(recon_file);
    write_off(os, reconstruction);
    os.close();

    // mark map for the faces we will delete to keep only the reconstruction
    Mark_map to_delete = get(Face_bool_tag(), reconstruction);
    // vector used as a stack to treat the faces to delete and add their neighbours to expand the deletion zone
    std::vector<face_descriptor<PolygonMesh>> to_handle;

    // 1. mark all faces incident to the borders, and on the opposite side of the reconstruction
    for(auto& hole : {hole_indices_1, hole_indices_2})
    {
      reconstruction_holes.emplace_back();
      unsigned size = hole.size();
      for(unsigned i = 0; i < size; ++i)
      {
        vertex_descriptor<PolygonMesh> v1 = *(reconstruction.vertices().begin() + hole[i]);
        vertex_descriptor<PolygonMesh> v2 = *(reconstruction.vertices().begin() + hole[(i + 1) % size]);
        halfedge_descriptor<PolygonMesh> h = reconstruction.halfedge(v1, v2);
        reconstruction_holes.back().push_back(h);
        put(to_delete, face(opposite(h, reconstruction), reconstruction), true);
      }
    }

    // 2. for each border, add one face f next to the marked faces next to the border
    //    f will be used as initial point to start the delition zone expantion
    for(auto& hole : {hole_indices_1, hole_indices_2})
    {
      vertex_descriptor<PolygonMesh> v1 = *(reconstruction.vertices().begin() + hole[0]);
      vertex_descriptor<PolygonMesh> v2 = *(reconstruction.vertices().begin() + hole[1]);
      halfedge_descriptor<PolygonMesh> h = reconstruction.halfedge(v1, v2);
      halfedge_descriptor<PolygonMesh> h_ = opposite(next(opposite(h, reconstruction), reconstruction), reconstruction);
      if(!get(to_delete, face(h_, reconstruction)))
      {
        put(to_delete, face(h_, reconstruction), true);
        to_handle.push_back(face(h_, reconstruction));
      }
      else
      {
        put(to_delete, face(opposite(next(opposite(h_, reconstruction), reconstruction), reconstruction), reconstruction), true);
        to_handle.push_back(face(opposite(next(opposite(h_, reconstruction), reconstruction), reconstruction), reconstruction));
      }
    }

    // 3. expand the zone
    while(!to_handle.empty())
    {
      face_descriptor<PolygonMesh> f = to_handle.back();
      to_handle.pop_back();
      halfedge_descriptor<PolygonMesh> h = opposite(halfedge(f, reconstruction), reconstruction);
      for(unsigned i = 0; i < 3; ++i, h = opposite(next(opposite(h, reconstruction), reconstruction), reconstruction))
      {
        if(get(to_delete, face(h, reconstruction)))
        {
          continue;
        }
        to_handle.push_back(face(h, reconstruction));
        put(to_delete, face(h, reconstruction), true);
      }
    }

    // 4. delete the expanded zone
    for(auto& f : reconstruction.faces())
    {
      if(get(to_delete, f))
      {
        CGAL::Euler::remove_face(halfedge(f, reconstruction), reconstruction);
      }
    }

    return reconstruction;
  }

  template <typename WorkZone, typename PolygonMesh>
  std::vector<halfedge_pair<PolygonMesh>> make_associations(const WorkZone& wz_1,
                                                            const WorkZone& wz_2,
                                                            const PolygonMesh& pmesh,
                                                            const std::vector<std::vector<halfedge_descriptor<PolygonMesh>>>& reconstruction_holes,
                                                            const std::vector<halfedge_pair<PolygonMesh>>& h2h)
  {
    std::vector<halfedge_pair<PolygonMesh>> associations;

    for(unsigned i = 0; i < wz_1.border.size(); ++i)
    {
      halfedge_descriptor<PolygonMesh> h1 = wz_1.border[i];
      halfedge_descriptor<PolygonMesh> h2 = reconstruction_holes[0][i];
      for(auto& it : h2h)
      {
        if(it.first == h2)
        {
          h2 = pmesh.opposite(it.second);
          break;
        }
      }
      associations.emplace_back(h1, h2);
    }

    for(unsigned i = 0; i < wz_2.border.size(); ++i)
    {
      halfedge_descriptor<PolygonMesh> h1 = wz_2.border[i];
      halfedge_descriptor<PolygonMesh> h2 = reconstruction_holes[1][i];
      for(auto& it : h2h)
      {
        if(it.first == h2)
        {
          h2 = pmesh.opposite(it.second);
          break;
        }
      }
      associations.emplace_back(h1, h2);
    }

    return associations;
  }
}

namespace Polygon_mesh_processing {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Combinatorial treatment

namespace internal {

template <typename G>
struct Vertex_collector
{
  typedef typename boost::graph_traits<G>::vertex_descriptor      vertex_descriptor;

  bool has_old_vertex(const vertex_descriptor v) const { return collections.count(v) != 0; }
  void tag_old_vertex(const vertex_descriptor v)
  {
    CGAL_precondition(!has_old_vertex(v));
    collections[v];
  }

  void collect_vertices(vertex_descriptor v1, vertex_descriptor v2)
  {
    std::vector<vertex_descriptor>& verts = collections[v1];
    if(verts.empty())
      verts.push_back(v1);
    verts.push_back(v2);
  }

  template<typename OutputIterator>
  void dump(OutputIterator out)
  {
    typedef std::pair<const vertex_descriptor, std::vector<vertex_descriptor> > Pair_type;
    for(const Pair_type& p : collections)
      *out++ = p.second;
  }

  void dump(Emptyset_iterator) { }

  std::map<vertex_descriptor, std::vector<vertex_descriptor> > collections;
};

template <typename PolygonMesh, typename VPM, typename ConstraintMap>
typename boost::graph_traits<PolygonMesh>::vertex_descriptor
create_new_vertex_for_sector(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor sector_begin_h,
                             typename boost::graph_traits<PolygonMesh>::halfedge_descriptor sector_last_h,
                             PolygonMesh& pmesh,
                             const VPM& vpm,
                             const ConstraintMap& cmap)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  vertex_descriptor old_vd = target(sector_begin_h, pmesh);
  vertex_descriptor new_vd = add_vertex(pmesh);
  put(vpm, new_vd, get(vpm, old_vd));

  put(cmap, new_vd, true);

  set_halfedge(new_vd, sector_begin_h, pmesh);
  halfedge_descriptor h = sector_begin_h;
  do
  {
    set_target(h, new_vd, pmesh);

    if(h == sector_last_h)
      break;
    else
      h = prev(opposite(h, pmesh), pmesh);
  }
  while(h != sector_begin_h); // for safety
  CGAL_assertion(h != sector_begin_h);

  return new_vd;
}

template <typename PolygonMesh, typename NamedParameters>
std::size_t make_umbrella_manifold(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                   PolygonMesh& pmesh,
                                   internal::Vertex_collector<PolygonMesh>& dmap,
                                   const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, pmesh));

  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_is_constrained_t,
                                                       NamedParameters,
                                                       Constant_property_map<vertex_descriptor, bool> // default (no constraint pmap)
                                                       >::type                  VerticesMap;
  VerticesMap cmap = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                                      Constant_property_map<vertex_descriptor, bool>(false));

  std::size_t nb_new_vertices = 0;

  vertex_descriptor old_v = target(h, pmesh);
  put(cmap, old_v, true); // store the duplicates

  // count the number of borders
  int border_counter = 0;
  halfedge_descriptor ih = h, done = ih, border_h = h;
  do
  {
    if(is_border(ih, pmesh))
    {
      border_h = ih;
      ++border_counter;
    }

    ih = prev(opposite(ih, pmesh), pmesh);
  }
  while(ih != done);

  bool is_non_manifold_within_umbrella = (border_counter > 1);
  if(!is_non_manifold_within_umbrella)
  {
    const bool first_time_meeting_v = !dmap.has_old_vertex(old_v);
    if(first_time_meeting_v)
    {
      // The star is manifold, so if it is the first time we have met that vertex,
      // there is nothing to do, we just keep the same vertex.
      set_halfedge(old_v, h, pmesh); // to ensure halfedge(old_v, pmesh) stays valid
      dmap.tag_old_vertex(old_v); // so that we know we have met old_v already, next time, we'll have to duplicate
    }
    else
    {
      // This is not the canonical star associated to 'v'.
      // Create a new vertex, and move the whole star to that new vertex
      halfedge_descriptor last_h = opposite(next(h, pmesh), pmesh);
      vertex_descriptor new_v = create_new_vertex_for_sector(h, last_h, pmesh, vpm, cmap);
      dmap.collect_vertices(old_v, new_v);
      nb_new_vertices = 1;
    }
  }
  // if there is more than one sector, look at each sector and split them away from the main one
  else
  {
    // the first manifold sector, described by two halfedges
    halfedge_descriptor sector_start_h = border_h;
    CGAL_assertion(is_border(border_h, pmesh));

    bool should_stop = false;
    bool is_main_sector = true;
    do
    {
      CGAL_assertion(is_border(sector_start_h, pmesh));

      // collect the sector and split it away if it must be
      halfedge_descriptor sector_last_h = sector_start_h;
      do
      {
        halfedge_descriptor next_h = prev(opposite(sector_last_h, pmesh), pmesh);

        if(is_border(next_h, pmesh))
          break;

        sector_last_h = next_h;
      }
      while(sector_last_h != sector_start_h);
      CGAL_assertion(!is_border(sector_last_h, pmesh));
      CGAL_assertion(sector_last_h != sector_start_h);

      halfedge_descriptor next_start_h = prev(opposite(sector_last_h, pmesh), pmesh);

      // there are multiple CCs incident to this particular vertex, and we should create a new vertex
      // if it's not the first umbrella around 'old_v' or not the first sector, but only not if it's
      // both the first umbrella and first sector.
      bool must_create_new_vertex = (!is_main_sector || dmap.has_old_vertex(old_v));

      // In any case, we must set up the next pointer correctly
      set_next(sector_start_h, opposite(sector_last_h, pmesh), pmesh);

      if(must_create_new_vertex)
      {
        vertex_descriptor new_v = create_new_vertex_for_sector(sector_start_h, sector_last_h, pmesh, vpm, cmap);
        dmap.collect_vertices(old_v, new_v);
        ++nb_new_vertices;
      }
      else
      {
        // We are in the first sector and first star, ensure that halfedge(old_v, pmesh) stays valid
        set_halfedge(old_v, sector_start_h, pmesh);
      }

      is_main_sector = false;
      sector_start_h = next_start_h;
      should_stop = (sector_start_h == border_h);
    }
    while(!should_stop);
  }

  return nb_new_vertices;
}

} // end namespace internal

/// \ingroup PMP_repairing_grp
/// duplicates all the non-manifold vertices of the input mesh.
///
/// @tparam PolygonMesh a model of `HalfedgeListGraph` and `MutableHalfedgeGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param pm the surface mesh to be repaired
/// @param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
///       The type of this map is model of `ReadWritePropertyMap`.
///       If this parameter is omitted, an internal property map for
///       `CGAL::vertex_point_t` should be available in `PolygonMesh`
///    \cgalParamEnd
///   \cgalParamBegin{vertex_is_constrained_map} a writable property map with `vertex_descriptor`
///     as key and `bool` as `value_type`. `put(pmap, v, true)` will be called for each duplicated
///     vertices, as well as the original non-manifold vertex in the input mesh.
///  \cgalParamEnd
///   \cgalParamBegin{output_iterator} a model of `OutputIterator` with value type
///      `std::vector<vertex_descriptor>`. The first vertex of each vector is a non-manifold vertex
///       of the input mesh, followed by the new vertices that were created to fix this precise
///       non-manifold configuration.
///  \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \return the number of vertices created.
template <typename PolygonMesh, typename NamedParameters>
std::size_t duplicate_non_manifold_vertices(PolygonMesh& pmesh,
                                            const NamedParameters& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef boost::graph_traits<PolygonMesh>                            GT;
  typedef typename GT::halfedge_descriptor                            halfedge_descriptor;

  typedef typename internal_np::Lookup_named_param_def<internal_np::output_iterator_t,
                                                       NamedParameters,
                                                       Emptyset_iterator>::type Output_iterator;

  Output_iterator out = choose_parameter(get_parameter(np, internal_np::output_iterator),
                                         Emptyset_iterator());

  std::vector<halfedge_descriptor> non_manifold_umbrellas;
  non_manifold_vertices(pmesh, std::back_inserter(non_manifold_umbrellas));

  internal::Vertex_collector<PolygonMesh> dmap;
  std::size_t nb_new_vertices = 0;
  if(!non_manifold_umbrellas.empty())
  {
    for(halfedge_descriptor h : non_manifold_umbrellas)
      nb_new_vertices += internal::make_umbrella_manifold(h, pmesh, dmap, np);

    dmap.dump(out);
  }

  return nb_new_vertices;
}

template <class PolygonMesh>
std::size_t duplicate_non_manifold_vertices(PolygonMesh& pmesh)
{
  return duplicate_non_manifold_vertices(pmesh, parameters::all_default());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Geometrical treatment

// main:
// @todo handle overlapping geodesic disks
//       --> in the future: identify matching zones and fill hole with screened poisson?
// @todo VPM for heat method / heat method on a FFG
//
// optimization:
// @todo can make maps of Points lighter with a vertex as key, and a custom equal comparing actual points
// @todo small sets

enum NM_TREATMENT
{
  SEPARATE = 0,
  MERGE
};

namespace internal {

template <typename PolygonMesh>
struct Work_zone
{//Les faces à tej puis reconstruire
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;

  typedef CGAL::dynamic_vertex_property_t<double>                             Distance_tag;
  typedef typename boost::property_map<PolygonMesh, Distance_tag>::type       Ditance_map;

  std::set<face_descriptor> faces;
  std::set<face_descriptor> ring; // only used for in meshing
  std::vector<halfedge_descriptor> border;
  Ditance_map distances;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Merging treatment

template <typename PolygonMesh>
typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
split_edge_and_triangulate_incident_faces(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                          PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  halfedge_descriptor res = Euler::split_edge(h, pmesh);

  if(!is_border(res, pmesh))
    Euler::split_face(res, next(h, pmesh), pmesh);

  halfedge_descriptor opp_h = opposite(h, pmesh);
  if(!is_border(opp_h, pmesh))
  {
    Euler::split_face(opp_h, next(next(opp_h, pmesh), pmesh), pmesh);
  }

  return res;
}

template <typename PolygonMesh>
bool is_pinched_nm_vertex(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                          const PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  int number_of_border_halfedges = 0;

  for(const halfedge_descriptor ih : halfedges_around_target(h, pmesh))
    if(is_border(ih, pmesh))
      ++number_of_border_halfedges;

  return (number_of_border_halfedges > 1);
}

template <typename UmbrellaContainer, typename PolygonMesh>
bool treat_pinched_stars(UmbrellaContainer& umbrellas,
                         const NM_TREATMENT /*treatment*/, // @todo what should be done in merge treatment...?
                         PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  bool has_pinched_nm_vertices = false;
  std::set<halfedge_descriptor> treated_umbrellas;

  for(const halfedge_descriptor h : umbrellas)
  {
    if(!is_pinched_nm_vertex(h, pmesh))
      continue;

    has_pinched_nm_vertices = true;
    treated_umbrellas.insert(h);

    std::vector<halfedge_descriptor> faces_to_delete;
    for(const halfedge_descriptor ih : halfedges_around_target(h, pmesh))
      if(!is_border(ih, pmesh))
        faces_to_delete.push_back(ih);

    for(halfedge_descriptor ih : faces_to_delete)
       Euler::remove_face(ih, pmesh);
  }

  for(halfedge_descriptor h : treated_umbrellas)
    umbrellas.erase(h);

  return has_pinched_nm_vertices;
}

// Merging combinatorially changes the star, so we don't want to have two adjacent nm vertices
template <typename NMVM, typename PolygonMesh, typename VPM, typename GeomTraits>
void enforce_non_manifold_vertex_separation(NMVM nm_marks,
                                            PolygonMesh& pmesh,
                                            VPM vpm,
                                            const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor          edge_descriptor;

  typedef typename boost::property_traits<VPM>::reference                     Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                    Point;

  std::set<edge_descriptor> edges_to_split;
  for(const edge_descriptor e : edges(pmesh))
  {
    if(get(nm_marks, source(e, pmesh)) && get(nm_marks, target(e, pmesh)))
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "edge that needs to be combinatorially separated: " << e
                << " vertices: " << source(e, pmesh) << " " << target(e, pmesh) << std::endl;
#endif
      edges_to_split.insert(e);
    }
  }

  for(const edge_descriptor e : edges_to_split)
  {
    halfedge_descriptor h = halfedge(e, pmesh);
    if(is_border(h, pmesh))
      h = opposite(h, pmesh);

    const vertex_descriptor vs = source(h, pmesh);
    const vertex_descriptor vt = target(h, pmesh);
    const Point_ref sp = get(vpm, vs);
    const Point_ref tp = get(vpm, vt);
    const Point mp = gt.construct_midpoint_3_object()(sp, tp);

    halfedge_descriptor new_h = split_edge_and_triangulate_incident_faces(h, pmesh);
    put(vpm, target(new_h, pmesh), mp);
  }
}

template <typename FT>
bool need_to_split(FT small_radius, FT radius, FT large_radius, FT min_distance, FT max_distance)
{
  return ((min_distance - small_radius) * (max_distance - small_radius) < 0)
      || ((min_distance - radius) * (max_distance - radius) < 0)
      || ((min_distance - large_radius) * (max_distance - large_radius) < 0);
//  if(min_distance > max_distance)
//  {
//    return need_to_split(small_radius, radius, large_radius, max_distance, min_distance);
//  }
//  // here we know min_distance <= max_distance
//  if(min_distance < small_radius)
//  {
//    return max_distance > small_radius;
//  }
//  // here we know min_distance >= small_radius
//  if(min_distance < radius)
//  {
//    return max_distance > radius;
//  }
//  // here we know min_distance >= radius
//  if(min_distance < large_radius)
//  {
//    return max_distance > large_radius;
//  }
//  // here we know min_distance >= arge_radius, hence non need to split
//  return false;
}

template <typename PolygonMesh, typename VPM, typename GeomTraits>
void sample_faces(const PolygonMesh& pmesh,
                  std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor> to_sample,
                  std::vector<typename GeomTraits::Point_3>& points,
                  std::vector<typename GeomTraits::Vector_3>& normals,
                  VPM vpm)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor          edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;
  typedef typename GeomTraits::Point_3                                        Point_3;

  std::cout << std::endl << std::endl;

  for(const auto& f : to_sample)
  {
    typename GeomTraits::Point_3 p, q, r;

    unsigned i = 0;
    for(const auto& v : vertices_around_face(halfedge(f, pmesh), pmesh))
    {
      switch(i){
        case 0:
          p = get(vpm, v);
          break;
        case 1:
          q = get(vpm, v);
          break;
        case 2:
          r = get(vpm, v);
          break;
        default:
          break;
      }
      ++i;
    }

    Random_points_in_triangle_3<Point_3> random(p, q, r);

    unsigned nb_points_on_f = 100;

    std::copy_n(random, nb_points_on_f, std::back_inserter(points));
    for(unsigned j = 0; j < nb_points_on_f; ++j)
    {
      normals.push_back(compute_face_normal(f, pmesh));
    }
  }
}

// note that this also refines the mesh
template <typename PolygonMesh, typename VPM, typename GeomTraits>
void construct_work_zone(Work_zone<PolygonMesh>& wz,
                         std::vector<typename GeomTraits::Point_3>& points,
                         std::vector<typename GeomTraits::Vector_3>& normals,
                         const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor zh,
                         const typename GeomTraits::FT radius,
                         PolygonMesh& pmesh,
                         VPM vpm,
                         const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor          edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;

  typedef typename boost::property_traits<VPM>::reference                     Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                    Point;

  typedef typename GeomTraits::FT                                             FT;
  typedef typename GeomTraits::Vector_3                                       Vector;

  typedef CGAL::dynamic_edge_property_t<bool>                                 Considered_tag;
  typedef typename boost::property_map<PolygonMesh, Considered_tag>::type     Considered_edge_map;

  typedef CGAL::dynamic_vertex_property_t<FT>                                 Distance_tag;

  std::ofstream os;


  const typename GeomTraits::FT small_radius = radius / 5;
  const typename GeomTraits::FT large_radius = radius * 1.5;

  wz.distances = get(Distance_tag(), pmesh);

  vertex_descriptor source_v = target(zh, pmesh);

  CGAL::Heat_method_3::estimate_geodesic_distances(pmesh, wz.distances, source_v); // @fixme vpm?

  // Corefine the mesh with a piecewise(face)-linear approximation of the geodesic circle
  std::set<face_descriptor> faces_to_consider;
  std::vector<halfedge_descriptor> edges_to_split;

  std::stack<halfedge_descriptor> halfedges_to_consider;
  halfedges_to_consider.push(opposite(zh, pmesh));

  Considered_edge_map considered_edges = get(Considered_tag(), pmesh);

  while(!halfedges_to_consider.empty())
  {
    const halfedge_descriptor curr_h = halfedges_to_consider.top();
    halfedges_to_consider.pop();

    const edge_descriptor curr_e = edge(curr_h, pmesh);
    put(considered_edges, curr_e, true);

    const FT source_distance = get(wz.distances, source(curr_e, pmesh));
    const FT target_distance = get(wz.distances, target(curr_e, pmesh));
    if(need_to_split(small_radius, radius, large_radius, source_distance, target_distance))
    {
      // We want to keep zh pointing to the nm vertex after the edge split
      edges_to_split.push_back((curr_h == opposite(zh, pmesh)) ? zh : curr_h);
    }

    if(!is_border(curr_h, pmesh))
      faces_to_consider.insert(face(curr_h, pmesh));
    if(!is_border(opposite(curr_h, pmesh), pmesh))
      faces_to_consider.insert(face(opposite(curr_h, pmesh), pmesh));

    const vertex_descriptor curr_v = source(curr_h, pmesh);
    if(get(wz.distances, curr_v) > large_radius)
      continue;

    for(const halfedge_descriptor adj_h : CGAL::halfedges_around_source(curr_h, pmesh))
    {
      if(get(considered_edges, edge(adj_h, pmesh))) // already visited
        continue;

      // want the source of the new halfedges to be equal to 'source(curr_h, pmesh)'
      halfedges_to_consider.push(opposite(adj_h, pmesh));
    }
  }

  std::cout << edges_to_split.size() << " edges to split" << std::endl;

  std::vector<halfedge_descriptor> next_halfedges;
  std::vector<halfedge_descriptor> next_next_halfedges;

  unsigned i1 = 0;
  unsigned i2 = 0;
  unsigned i3 = 0;

  // Actual split
  for(halfedge_descriptor h : edges_to_split) {
    const vertex_descriptor vs = source(h, pmesh);
    const vertex_descriptor vt = target(h, pmesh);
    const Point_ref spt = get(vpm, vs);
    const Point_ref tpt = get(vpm, vt);
    Vector tsv = gt.construct_vector_3_object()(tpt, spt);

    const FT dist_at_vs = get(wz.distances, vs);
    const FT dist_at_vt = get(wz.distances, vt);

//    if(dist_at_vs == radius || dist_at_vt == radius) // nothing to do // will never occur anymore, see function need_to_split
//      continue;

//  FIRST SPLIT
    if ((dist_at_vs - large_radius) * (dist_at_vt - large_radius) < 0) // need to split for the big sphere
    {
      ++i1;
      Point new_p;
      if (dist_at_vs < dist_at_vt) {
        CGAL_assertion(dist_at_vs < large_radius && large_radius <= dist_at_vt);
        const FT lambda = (large_radius - dist_at_vs) / (dist_at_vt - dist_at_vs);
        new_p = gt.construct_translated_point_3_object()(
            spt, gt.construct_scaled_vector_3_object()(tsv, -lambda));
      } else {
        CGAL_assertion(dist_at_vt < large_radius && large_radius <= dist_at_vs);
        const FT lambda = (large_radius - dist_at_vt) / (dist_at_vs - dist_at_vt);
        new_p = gt.construct_translated_point_3_object()(
            tpt, gt.construct_scaled_vector_3_object()(tsv, lambda));
      }

      halfedge_descriptor new_h = split_edge_and_triangulate_incident_faces(h, pmesh);

      if (dist_at_vs < dist_at_vt)
        //    non mn vtx     radius
        //        *     --------|-------->  h
      {
        next_halfedges.push_back(new_h);
      } else
        //    non mn vtx     radius
        //        *     <-------|---------  h
      {
        next_halfedges.push_back(h);
      }

      put(vpm, target(new_h, pmesh), new_p);
      put(wz.distances, target(new_h, pmesh), large_radius);

      // @todo might be simplifiable?
      if (!is_border(h, pmesh)) {
        faces_to_consider.insert(face(h, pmesh));
        faces_to_consider.insert(face(new_h, pmesh));
      }

      if (!is_border(opposite(h, pmesh), pmesh)) {
        faces_to_consider.insert(face(opposite(h, pmesh), pmesh));
        faces_to_consider.insert(face(opposite(new_h, pmesh), pmesh));
      }
    }
    else
    {
      next_halfedges.push_back(h);
    }
  }


//  SECOND SPLIT
  for(halfedge_descriptor next_h : next_halfedges)
  {
    ++i2;
    const vertex_descriptor vs = source(next_h, pmesh);
    const vertex_descriptor vt = target(next_h, pmesh);
    const Point_ref spt = get(vpm, vs);
    const Point_ref tpt = get(vpm, vt);
    Vector tsv = gt.construct_vector_3_object()(tpt, spt);

    const FT dist_at_vs = get(wz.distances, vs);
    const FT dist_at_vt = get(wz.distances, vt);

    if((dist_at_vs - radius) * (dist_at_vt - radius) < 0) // need to split for the big sphere
    {
      Point new_p;
      if (dist_at_vs < dist_at_vt) {
        CGAL_assertion(dist_at_vs < radius && radius <= dist_at_vt);
        const FT lambda = (radius - dist_at_vs) / (dist_at_vt - dist_at_vs);
        new_p = gt.construct_translated_point_3_object()(
            spt, gt.construct_scaled_vector_3_object()(tsv, -lambda));
      } else {
        CGAL_assertion(dist_at_vt < radius && radius <= dist_at_vs);
        const FT lambda = (radius - dist_at_vt) / (dist_at_vs - dist_at_vt);
        new_p = gt.construct_translated_point_3_object()(
            tpt, gt.construct_scaled_vector_3_object()(tsv, lambda));
      }

      halfedge_descriptor new_h = split_edge_and_triangulate_incident_faces(next_h, pmesh);

      if (dist_at_vs < dist_at_vt)
        //    non mn vtx     radius
        //        *     --------|-------->  h
      {
        next_next_halfedges.push_back(new_h);
      } else
        //    non mn vtx     radius
        //        *     <-------|---------  h
      {
        next_next_halfedges.push_back(next_h);
      }

      put(vpm, target(new_h, pmesh), new_p);
      put(wz.distances, target(new_h, pmesh), radius);

      // @todo might be simplifiable?
      if (!is_border(next_h, pmesh)) {
        faces_to_consider.insert(face(next_h, pmesh));
        faces_to_consider.insert(face(new_h, pmesh));
      }

      if (!is_border(opposite(next_h, pmesh), pmesh)) {
        faces_to_consider.insert(face(opposite(next_h, pmesh), pmesh));
        faces_to_consider.insert(face(opposite(new_h, pmesh), pmesh));
      }
    }
    else
    {
      next_next_halfedges.push_back(next_h);
    }
  }

  os = std::ofstream("/home/felix/Bureau/Geo_Facto/PSR/tests-repair/dumps/split/second-split.off");
  CGAL::write_off(os, pmesh);
  os.close();

// THIRD SPLIT
  for(halfedge_descriptor next_h : next_next_halfedges)
  {
    ++i3;
    const vertex_descriptor vs = source(next_h, pmesh);
    const vertex_descriptor vt = target(next_h, pmesh);
    const Point_ref spt = get(vpm, vs);
    const Point_ref tpt = get(vpm, vt);
    Vector tsv = gt.construct_vector_3_object()(tpt, spt);

    const FT dist_at_vs = get(wz.distances, vs);
    const FT dist_at_vt = get(wz.distances, vt);

    if((dist_at_vs - small_radius) * (dist_at_vt - small_radius) < 0) // need to split for the big sphere
    {
      Point new_p;
      if (dist_at_vs < dist_at_vt) {
        CGAL_assertion(dist_at_vs < small_radius && small_radius <= dist_at_vt);
        const FT lambda = (small_radius - dist_at_vs) / (dist_at_vt - dist_at_vs);
        new_p = gt.construct_translated_point_3_object()(
            spt, gt.construct_scaled_vector_3_object()(tsv, -lambda));
      } else {
        CGAL_assertion(dist_at_vt < small_radius && small_radius <= dist_at_vs);
        const FT lambda = (small_radius - dist_at_vt) / (dist_at_vs - dist_at_vt);
        new_p = gt.construct_translated_point_3_object()(
            tpt, gt.construct_scaled_vector_3_object()(tsv, lambda));
      }

      halfedge_descriptor new_h = split_edge_and_triangulate_incident_faces(next_h, pmesh);

      put(vpm, target(new_h, pmesh), new_p);
      put(wz.distances, target(new_h, pmesh), small_radius);

      // @todo might be simplifiable?
      if (!is_border(next_h, pmesh)) {
        faces_to_consider.insert(face(next_h, pmesh));
        faces_to_consider.insert(face(new_h, pmesh));
      }

      if (!is_border(opposite(next_h, pmesh), pmesh)) {
        faces_to_consider.insert(face(opposite(next_h, pmesh), pmesh));
        faces_to_consider.insert(face(opposite(new_h, pmesh), pmesh));
      }
    }
  }


#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << faces_to_consider.size() << " faces to consider" << std::endl;

  static int fi = 0;
  std::stringstream oss;
  oss << "results/faces_to_consider_" << fi++ << ".off" << std::ends;
  dump_cc(faces_to_consider, pmesh, oss.str().c_str());
#endif

  std::set<face_descriptor> to_sample;

  for(const face_descriptor f : faces_to_consider)
  {
    bool is_face_in_ring = true;
    bool is_face_in = true;
    bool need_sample = true;
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, pmesh), pmesh))
    {
      if(get(wz.distances, target(h, pmesh)) > large_radius || get(wz.distances, target(h, pmesh)) < radius)
      {
        is_face_in_ring = false;
      }
      if(get(wz.distances, target(h, pmesh)) > radius)
      {
        is_face_in = false;
        need_sample = false;
      }
      if(get(wz.distances, target(h, pmesh)) < small_radius)
      {
        need_sample = false;
      }
    }

    if(is_face_in_ring)
      wz.ring.insert(f);
    if(is_face_in)
      wz.faces.insert(f);
    if(need_sample)
      to_sample.insert(f);
  }

//  sample_faces<PolygonMesh, VPM, GeomTraits>(pmesh, to_sample, points, normals, vpm);

}

template <typename PolygonMesh, typename VPM, typename GeomTraits>
Work_zone<PolygonMesh> construct_work_zone(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                           const typename GeomTraits::FT radius,
                                           PolygonMesh& pmesh,
                                           std::vector<typename GeomTraits::Point_3>& points,
                                           std::vector<typename GeomTraits::Vector_3>& normals,
                                           VPM vpm,
                                           const GeomTraits& gt)
{
  // 'set' complexity should be fine, it's not supposed to be a large number of faces
  Work_zone<PolygonMesh> wz;
  construct_work_zone(wz, points, normals, h, radius, pmesh, vpm, gt);
  CGAL_assertion(!wz.faces.empty());

  // If the selection is not a topological disk, just don't do anything
  // because we won't know how to fill it
  if(!is_selection_a_topological_disk(wz.faces, pmesh))
  {
    std::cerr << "Warning: selection is not a topological disk" << std::endl;
    return Work_zone<PolygonMesh>();
  }

  return wz;
}

template <typename PolygonMesh, typename VPM, typename GeomTraits>
void extract_border_of_work_zone(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                 Work_zone<PolygonMesh>& wz,
                                 const PolygonMesh& pmesh,
                                 const VPM vpm,
                                 const GeomTraits& /*gt*/) // @todo
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  typedef typename boost::property_traits<VPM>::reference                     Point_ref;

  wz.border.reserve(wz.faces.size());

  CGAL::Face_filtered_graph<PolygonMesh> ffg(pmesh, wz.faces);

  halfedge_descriptor bh = boost::graph_traits<PolygonMesh>::null_halfedge();
  for(halfedge_descriptor h : halfedges(ffg))
  {
    if(is_border(h, ffg))
    {
      bh = h;
      break;
    }
  }

  CGAL_assertion(bh != boost::graph_traits<PolygonMesh>::null_halfedge());

  // Below walks the border of the filtered graph, with a twist: if the star is open, we don't
  // care about the part of that border incident to the vertex:

  /*
      __________________               ______________
     |                  |             |              |
     |                  |             |              |
     |     nm vertex    |  should be  |              |
     |       /  \       |  -------->  |              |
     |      /    \      |             |              |
     |____b/      \c____|             |____b    c____|
  */

  // In other words, a sequence of border halfedges that are border halfedges for pmesh is only
  // acceptable if the nm vertex is not part of it.

  // start by trying to find a halfedge that is not on the border of pmesh to get full border sequences
  const Point_ref nmp = get(vpm, target(h, pmesh));
  halfedge_descriptor nmp_ih = boost::graph_traits<PolygonMesh>::null_halfedge();

  halfedge_descriptor done = bh;
  do
  {
    if(get(vpm, target(bh, pmesh)) == nmp)
      nmp_ih = bh;

    if(!is_border(bh, pmesh))
      break;
    bh = prev(bh, ffg);
  }
  while(bh != done);

  if(is_border(bh, pmesh) && // couldn't find a non-border halfedge: 'ffg' is a full CC of 'pmesh'
     nmp_ih != boost::graph_traits<PolygonMesh>::null_halfedge()) // nm vertex is on the border
  {
    // not really sure what to do here, so let's just remove the two halfedges
    // incident to the vertex and keep the rest

    halfedge_descriptor done = bh;
    {
      if(get(vpm, source(bh, pmesh)) != nmp && get(vpm, target(bh, pmesh)) != nmp)
        wz.border.push_back(opposite(bh, pmesh));

      bh = prev(bh, ffg); // walk backwards so that opposites are in the correct order
    }
    while(bh != done);
  }
  else // standard case now
  {
    halfedge_descriptor done = bh;
    do
    {
      if(is_border(bh, pmesh))
      {
        // start of a sequence of border (for pmesh) halfedges
        // we only want to add that to the complete border if it does not contain the nm vertex
        std::vector<halfedge_descriptor> seq;
        bool add_to_sequence = true;

        do
        {
          if(!is_border(bh, pmesh)) // end of the sequence
            break;

          std::cout << "seq hd: " << get(vpm, source(bh, pmesh)) << " || " << get(vpm, target(bh, pmesh)) << std::endl;
          std::cout << "nmp: " << nmp << std::endl;
          std::cout << "get(vpm, source(bh, pmesh)) == nmp : " << (get(vpm, source(bh, pmesh)) == nmp) << std::endl;

          if(get(vpm, source(bh, pmesh)) == nmp)
            add_to_sequence = false;
          else if(add_to_sequence)
            seq.push_back(opposite(bh, pmesh));

          bh = prev(bh, ffg); // walk backwards so that opposites are in the correct order
        }
        while(bh != done);

        CGAL_assertion(!is_border(bh, pmesh) ||
                       (bh == done && is_border(prev(bh, ffg), pmesh)));

        std::cout << "add_to_sequence = " << add_to_sequence << std::endl;
        if(add_to_sequence)
          wz.border.insert(wz.border.end(), std::begin(seq), std::end(seq)); // not worth move-ing pointers
      }
      else
      {
        wz.border.push_back(opposite(bh, ffg));
        bh = prev(bh, ffg); // walk backwards so that opposites are in the correct order
      }
    }
    while(bh != done);
  }
}

// compute the faces to remove
template <typename UmbrellaContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
std::vector<Work_zone<PolygonMesh> >
construct_work_zones(UmbrellaContainer& umbrellas,
                     const typename GeomTraits::FT radius,
                     PolygonMesh& pmesh,
                     std::vector<typename GeomTraits::Point_3>& points,
                     std::vector<typename GeomTraits::Vector_3>& normals,
                     VPM vpm,
                     const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  std::vector<Work_zone<PolygonMesh> > wzs;
  wzs.reserve(umbrellas.size());

  for(const halfedge_descriptor h : umbrellas)
  {
    std::cout << "### treating zone incident to: " << get(vpm, source(h, pmesh)) << " " << get(vpm, target(h, pmesh)) << std::endl;

    Work_zone<PolygonMesh> wz = construct_work_zone(h, radius, pmesh, points, normals, vpm, gt);
    if(!wz.faces.empty())
      wzs.push_back(wz);

    std::cout << "### extract border of zone incident to: " << get(vpm, source(h, pmesh)) << " " << get(vpm, target(h, pmesh)) << std::endl;

    extract_border_of_work_zone(h, wzs.back(), pmesh, vpm, gt);

    std::cout << "border of zone " << wzs.size() - 1 << std::endl;
    for(const halfedge_descriptor bh : wzs.back().border)
      std::cout << get(vpm, source(bh, pmesh)) << " " << get(vpm, target(bh, pmesh)) << std::endl;
  }

  return wzs;
}

template <typename HolePointContainer, typename ThirdPointContainer, typename FaceIndexContainer>
bool fill_polyline_hole(const HolePointContainer& hole_points,
                        const ThirdPointContainer& third_points,
                        FaceIndexContainer& patch)
{
  triangulate_hole_polyline(hole_points, third_points, std::back_inserter(patch));

  if(patch.empty())
  {
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Failed to fill a hole using Delaunay search space.\n";
#endif

    triangulate_hole_polyline(hole_points, third_points, std::back_inserter(patch),
                              parameters::use_delaunay_triangulation(false));
#endif // CGAL_HOLE_FILLING_DO_NOT_USE_DT3
    if(patch.empty())
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "Failed to fill a hole using the whole search space.\n";
#endif
      return false;
    }
  }

  return true;
}

// Currently performed with a hack: transform the borders into a single border
// by manually adding a patch between the closest edges
template <typename HalfedgeContainer_A, typename HalfedgeContainer_B,
          typename VPM_A, typename VPM_B,
          typename PolygonMesh,
          typename OutputIterator,
          typename GeomTraits>
bool two_borders_hole_fill(const HalfedgeContainer_A& bhv_A,
                           const PolygonMesh& pmesh_A,
                           const VPM_A vpm_A,
                           const HalfedgeContainer_B& bhv_B,
                           const PolygonMesh& pmesh_B,
                           const VPM_B vpm_B,
                           OutputIterator out,
                           const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  typedef typename boost::property_traits<VPM_A>::value_type                  Point;
  typedef typename boost::property_traits<VPM_A>::reference                   Point_ref_A;
  typedef typename boost::property_traits<VPM_B>::reference                   Point_ref_B;

  typedef typename HalfedgeContainer_A::const_iterator                        HCit_A;
  typedef typename HalfedgeContainer_B::const_iterator                        HCit_B;

  typedef CGAL::Triple<int, int, int>                                         Face_indices;

  CGAL_precondition(bhv_A.size() >= 3);
  CGAL_precondition(bhv_B.size() >= 3);

  HCit_A canon_A;
  HCit_B canon_B;

  // @todo avoid the O(n^2) complexity (but the borders are likely small, so...)
  // @todo this type of matching might not be the best idea

  double best_score = std::numeric_limits<double>::max();
  for(HCit_A it_A=std::begin(bhv_A), A_end=std::end(bhv_A); it_A!=A_end; ++it_A)
  {
    const Point_ref_A sap = get(vpm_A, source(*it_A, pmesh_A));
    const Point_ref_A tap = get(vpm_A, target(*it_A, pmesh_A));

    for(HCit_B it_B=std::begin(bhv_B), B_end=std::end(bhv_B); it_B!=B_end; ++it_B)
    {
      const Point_ref_B sbp = get(vpm_B, source(*it_B, pmesh_B));
      const Point_ref_B tbp = get(vpm_B, target(*it_B, pmesh_B));

      const double score = gt.compute_squared_distance_3_object()(sap, tbp)
                         + gt.compute_squared_distance_3_object()(tap, sbp);
      if(score < best_score)
      {
        best_score = score;
        canon_A = it_A;
        canon_B = it_B;
      }
    }
  }

  // polyline
  std::vector<Point> hole_points, third_points;

  // _____   _________   ________
  //      | | <------ | |
  //      | | canon_A | |
  //      | |         | |
  //      | | canon_B | |
  //      | | ------> | |
  // -----   --------   -------

  const Point_ref_A sap = get(vpm_A, source(*canon_A, pmesh_A));
  const Point_ref_A tap = get(vpm_A, target(*canon_A, pmesh_A));
  const Point_ref_B sbp = get(vpm_B, source(*canon_B, pmesh_B));
  const Point_ref_B tbp = get(vpm_B, target(*canon_B, pmesh_B));

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "best score:" << std::endl;
  std::cout << "sap: " << sap << std::endl;
  std::cout << "tap: " << tap << std::endl;
  std::cout << "sbp: " << sbp << std::endl;
  std::cout << "tbp: " << tbp << std::endl;
#endif

  std::ofstream border_out("results/hole_border.polylines.txt");
  border_out.precision(17);

  // Only a single ID is needed, the rest is 0, hole.size()-1, and quad_id + 1
  std::size_t quad_id = static_cast<std::size_t>(-1);

  // Walk A's border
  HCit_A last_A = std::prev(bhv_A.end());
  HCit_A it_A = (canon_A == last_A) ? std::begin(bhv_A) : std::next(canon_A);
  do
  {
    const halfedge_descriptor h = *it_A;

    if(!hole_points.empty())
      border_out << "2 " << hole_points.back() << " " << get(vpm_A, source(h, pmesh_A)) << std::endl;

    hole_points.push_back(get(vpm_A, source(h, pmesh_A)));
    std::cout << "point on hole: " << hole_points.back() << std::endl;

    const halfedge_descriptor oh = opposite(h, pmesh_A);
    if(is_border(oh, pmesh_A))
      third_points.push_back(construct_artificial_third_point(h, pmesh_A, vpm_A, gt));
    else
      third_points.push_back(get(vpm_A, target(next(oh, pmesh_A), pmesh_A)));
    it_A = (it_A == last_A) ? std::begin(bhv_A) : std::next(it_A);
  }
  while(it_A != canon_A);

  // vertical
  border_out << "2 " << hole_points.back() << " " << sap << std::endl;

  quad_id = hole_points.size();
  hole_points.push_back(sap);
  third_points.push_back(CGAL::midpoint(sbp, tap));

  // Walk B's border
  HCit_B last_B = std::prev(bhv_B.end());
  HCit_B it_B = (canon_B == last_B) ? std::begin(bhv_B) : std::next(canon_B);
  do
  {
    const halfedge_descriptor h = *it_B;

    border_out << "2 " << hole_points.back() << " " << get(vpm_B, source(h, pmesh_B)) << std::endl;
    hole_points.push_back(get(vpm_B, source(h, pmesh_B)));

    const halfedge_descriptor oh = opposite(h, pmesh_B);
    if(is_border(oh, pmesh_A))
      third_points.push_back(construct_artificial_third_point(h, pmesh_B, vpm_B, gt));
    else
      third_points.push_back(get(vpm_B, target(next(opposite(h, pmesh_B), pmesh_B), pmesh_B)));
    it_B = (it_B == last_B) ? std::begin(bhv_B) : std::next(it_B);
  }
  while(it_B != canon_B);

  border_out << "2 " << hole_points.back() << " " << sbp << std::endl;
  hole_points.push_back(sbp);
  third_points.push_back(CGAL::midpoint(tbp, sap));

  CGAL_assertion(hole_points.size() == third_points.size());

  std::vector<Face_indices> patch;
  fill_polyline_hole(hole_points, third_points, patch);

  // add the missing quad
  patch.emplace_back(quad_id, 0, quad_id+1);
  patch.emplace_back(quad_id+1, 0, hole_points.size() - 1);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::ofstream pout("results/patch.off");
  pout << std::setprecision(17);
  pout << 3 * patch.size() << " " << patch.size() << " 0\n";

  for(const Face_indices& f : patch)
  {
    pout << hole_points[f.first] << "\n";
    pout << hole_points[f.second] << "\n";
    pout << hole_points[f.third] << "\n";
  }

  for(std::size_t i=0, ps=patch.size(); i<ps; ++i)
    pout << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
#endif

  for(const Face_indices& face : patch)
  {
    *out++ = std::initializer_list<Point>{ hole_points[face.first],
                                           hole_points[face.second],
                                           hole_points[face.third] };
  }

  return true;
}

template <typename Patch,
          typename PolygonMesh,
          typename VPM,
          typename GeomTraits>
bool fix_patch_orientation(Patch& point_patch,
                           const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h1,
                           const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h2,
                           const PolygonMesh& pmesh,
                           const VPM vpm,
                           const GeomTraits& gt)
{
  typedef typename boost::property_traits<VPM>::reference                       Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  bool ok_orientation_1 = false, ok_orientation_2 = false;

  const Point_ref h1sp = get(vpm, source(h1, pmesh));
  const Point_ref h1tp = get(vpm, target(h1, pmesh));
  const Point_ref h2sp = get(vpm, source(h2, pmesh));
  const Point_ref h2tp = get(vpm, target(h2, pmesh));

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "h1sp: " << h1sp << std::endl;
  std::cout << "h1tp: " << h1tp << std::endl;
  std::cout << "h2sp: " << h2sp << std::endl;
  std::cout << "h2tp: " << h2tp << std::endl;
#endif

  for(const auto& face : point_patch)
  {
    for(int i=0; i<3; ++i)
    {
      const Point& p1 = face[i];
      const Point& p2 = face[(i+1)%3];

      if(gt.equal_3_object()(p1, h1sp) && gt.equal_3_object()(p2, h1tp))
        ok_orientation_1 = true;
      if(gt.equal_3_object()(p1, h2sp) && gt.equal_3_object()(p2, h2tp))
        ok_orientation_2 = true;
    }
  }

  std::cout << "orientations: " << ok_orientation_1 << " " << ok_orientation_2 << std::endl;

  if(ok_orientation_1 != ok_orientation_2)
    return false;

  if(!ok_orientation_1)
  {
    CGAL_assertion(false); // @tmp
    for(auto& face : point_patch)
      std::swap(face[0], face[1]);
  }

  return true;
}

template <typename WorkZone, typename PolygonMesh, typename VPM, typename GeomTraits>
bool merge_zones(const WorkZone& wz_1,
                 const WorkZone& wz_2,
                 PolygonMesh& pmesh,
                 std::vector<typename GeomTraits::Point_3>& points,
                 std::vector<typename GeomTraits::Vector_3>& normals,
                 const VPM vpm,
                 const GeomTraits& gt)
{
  namespace IMT = CGAL::In_meshing_tools_for_repair_manifoldness;

  std::ofstream os("/home/felix/Bureau/Geo_Facto/PSR/tests-repair/dumps/avant.off");
  write_off(os, pmesh);
  os.close();

  output out;
  out.enable_export = true;
  out.enable_log_data = true;
  out.enable_log_quality = true;
  out.enable_log_timing = true;
  bool allow_recompute_normals = false;
  unsigned max_depth = 8;
  bool export_raw_poisson_surface = true;
  bool allow_same_boundary_chains = false;
  float salient_angle_deg = 84;

  out.start_timing();

  std::vector<unsigned> hole_indices_1, hole_indices_2;
  point_mesh pm = IMT::make_point_mesh_for_in_meshing(wz_1, wz_2, pmesh, vpm, gt, hole_indices_1, hole_indices_2, points, normals);
//  pm.add_guide("/home/felix/Bureau/Geo_Facto/PSR/tests-repair/jeux-de-test/guide_spherique_1p2.ply");

  IMT::erase_wz_faces(wz_1, wz_2, pmesh);

  out.back_transform = pm.transform_to_unit(1.25f);
  out.stop_timing("Loading & transforming points");

  if(allow_recompute_normals)
  {
    size_t invalid_normals = 0;
    for(size_t v = 0; v < pm.vertices.size(); v++)
      if(pm.vertices[v].normal.squaredNorm() <= 0.1)
        invalid_normals++;
    if(invalid_normals > 0)
    {
      printf("%i invalid normals => recompute them\n", (int)invalid_normals);
      const std::vector<Eigen::Vector3f> normals = compute_normals(pm.vertices, pm.indices);
      for(size_t v = 0; v < pm.vertices.size(); v++)
        pm.vertices[v].normal = normals[v];
    }
  }

  /// 1) Mesh completion
  completed_mesh out_mesh = std::move(pm);
  mesh_in_filling(out_mesh, out, max_depth, export_raw_poisson_surface);

  /// 2) Large scale geometric adjustments
  out.log_metric("Total vertices", out_mesh.vertices.size());
  out.log_metric("Total triangles", out_mesh.indices.size() / 3);
  out.log_metric("Base vertices", out_mesh.num_base_vertices);
  out.log_metric("Base triangles", out_mesh.num_base_triangles);
//  geometry_new(out, out_mesh, allow_same_boundary_chains, salient_angle_deg, hole_indices_1, hole_indices_2);
  float target_length = geometry_new1(out, out_mesh, allow_same_boundary_chains, salient_angle_deg, hole_indices_1, hole_indices_2);
  IMT::remesh_with_CGAL<GeomTraits>(out_mesh, target_length * 3);
  geometry_new2(out, out_mesh, allow_same_boundary_chains, salient_angle_deg, hole_indices_1, hole_indices_2);

  /// 3) Final output
  out.enable_export = true;
  out.log_metric("Output/Vertices", out_mesh.vertices.size());
  out.log_metric("Output/Triangles", out_mesh.indices.size() / 3);
//  out.save_object("Final object", cmd.output_filename(), out_mesh);

  typedef IMT::halfedge_descriptor<PolygonMesh> halfedge_descriptor;
  std::vector<std::vector<halfedge_descriptor>> reconstruction_holes;
  PolygonMesh reconstruction = IMT::extract_reconstruction<PolygonMesh, GeomTraits>
      (out, out_mesh, gt, hole_indices_1, hole_indices_2, reconstruction_holes);

  typedef std::pair<halfedge_descriptor, halfedge_descriptor> halfedge_pair;
  std::vector<halfedge_pair> h2h;
  CGAL::copy_face_graph(reconstruction, pmesh, CGAL::parameters::halfedge_to_halfedge_output_iterator(std::back_inserter(h2h)));

  std::vector<halfedge_pair> associations = IMT::make_associations(wz_1, wz_2, pmesh, reconstruction_holes, h2h);
  stitch_borders(pmesh, associations);

  os = std::ofstream("/home/felix/Bureau/Geo_Facto/PSR/tests-repair/dumps/apres.off");
  write_off(os, pmesh);
  os.close();

  return true;
}

template <typename WorkZone, typename PolygonMesh, typename VPM, typename GeomTraits>
bool merge_zones_(const WorkZone& wz_1,
                  const WorkZone& wz_2,
                  PolygonMesh& pmesh,
                  const VPM vpm,
                  const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor          face_descriptor;

  typedef typename GeomTraits::Point_3                                        Point;

  const std::vector<halfedge_descriptor>& bhv_1 = wz_1.border;
  const std::vector<halfedge_descriptor>& bhv_2 = wz_2.border;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << "Merge holes (" << wz_1.faces.size() << " and " << wz_2.faces.size() << ")" << std::endl;
  std::cout << "holes of size: " << bhv_1.size() << " && " << bhv_2.size() << std::endl;
#endif

  // make sure that the holes are topological disks
  std::vector<std::vector<Point> > point_patch;
  point_patch.reserve(2 * bhv_1.size());

  bool success = two_borders_hole_fill(bhv_1, pmesh, vpm, bhv_2, pmesh, vpm,
                                       std::back_inserter(point_patch), gt);
  if(!success)
    return false;

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << point_patch.size() << " new faces" << std::endl;
  static int merge_id = 0;
  std::stringstream oss;
  oss << "results/tentative_merge_patch_" << merge_id << ".off" << std::ends;
  dump_tentative_hole(point_patch, oss.str().c_str());
#endif

  if(!check_patch_sanity<PolygonMesh>(point_patch))
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "this patch is INSANE!" << std::endl;
#endif
    return false;
  }

  // @todo shouldn't the orientation always be correct by construction?
  // Do this with combinatorics if it turns out it can be not correct
  success = fix_patch_orientation(point_patch, bhv_1.front(), bhv_2.front(), pmesh, vpm, gt);
  if(!success)
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Incompatible orientations" << std::endl;
#endif
    return false; // can't find an orientation compatible with both borders
  }

  std::set<face_descriptor> fs(std::begin(wz_1.faces), std::end(wz_1.faces));
  fs.insert(std::begin(wz_2.faces), std::end(wz_2.faces));

  return replace_faces_with_patch(fs, point_patch, pmesh, vpm);
}

/*
 The idea is to walk the polylines incident to the nm vertex that are on the border of the mesh,
 and modify the boundary so that it doesn't pass through that vertex anymore.
 for example:

*/
template <typename PolygonMesh, typename VPM, typename GeomTraits>
bool fill_zone_with_open_border(const Work_zone<PolygonMesh>& wz,
                                PolygonMesh& pmesh,
                                const VPM vpm,
                                const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  typedef typename boost::property_traits<VPM>::value_type                    Point;
  typedef typename boost::property_traits<VPM>::reference                     Point_ref;

  typedef CGAL::Triple<int, int, int>                                         Face_indices;

  std::vector<Point> hole_points, third_points;
  for(const halfedge_descriptor bh : wz.border)
  {
    const halfedge_descriptor obh = opposite(bh, pmesh);

    hole_points.push_back(get(vpm, source(bh, pmesh)));

    if(is_border(obh, pmesh))
      third_points.push_back(construct_artificial_third_point(bh, pmesh, vpm, gt));
    else
      third_points.push_back(get(vpm, target(next(obh, pmesh), pmesh)));
  }

  // last one
  const halfedge_descriptor lh = wz.border.back();
  hole_points.push_back(get(vpm, target(lh, pmesh)));

  const halfedge_descriptor olh = opposite(lh, pmesh);
  if(is_border(olh, pmesh))
    third_points.push_back(construct_artificial_third_point(lh, pmesh, vpm, gt));
  else
    third_points.push_back(get(vpm, target(next(olh, pmesh), pmesh)));

  CGAL_assertion(hole_points.size() == third_points.size());

  std::vector<Face_indices> patch;
  fill_polyline_hole(hole_points, third_points, patch);

  std::vector<std::vector<Point> > point_patch;
  point_patch.reserve(patch.size());
  for(const Face_indices& face : patch)
  {
    point_patch.emplace_back(std::initializer_list<Point>{ hole_points[face.first],
                                                           hole_points[face.second],
                                                           hole_points[face.third] });
  }

  return replace_faces_with_patch(wz.faces, point_patch, pmesh, vpm);
}

template <typename WorkZoneContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
bool fill_zones(const WorkZoneContainer& wzs,
                PolygonMesh& pmesh,
                VPM vpm,
                const GeomTraits& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor        halfedge_descriptor;

  for(std::size_t i=0, n=wzs.size(); i<n; ++i)
  {
    CGAL_assertion(!wzs[i].faces.empty());
    CGAL_assertion(!wzs[i].border.empty());

    if(wzs[i].faces.size() == 1)
    {
      Euler::remove_face(halfedge(*(std::begin(wzs[i].faces)), pmesh), pmesh);
      continue;
    }

    const halfedge_descriptor f = wzs[i].border.front();
    const halfedge_descriptor b = wzs[i].border.back();

    const bool is_loop = source(f, pmesh) == target(b, pmesh);
    if(!is_loop)
    {
      fill_zone_with_open_border(wzs[i], pmesh, vpm, gt);
    }
    else if(!fill_hole(wzs[i].faces, pmesh, vpm, gt))
    {
      std::cerr << "Failed to fill hole of work zone #" << i << std::endl;
      return false;
    }
  }

  return true;
}

template <typename UmbrellaContainer, typename PolygonMesh, typename VPM, typename GeomTraits>
bool treat_umbrellas(UmbrellaContainer& umbrellas,
                     const NM_TREATMENT treatment,
                     const typename GeomTraits::FT radius,
                     PolygonMesh& pmesh,
                     VPM vpm,
                     const GeomTraits& gt)
{
  // can't merge if there were any pinched stars
  const bool can_employ_merge_strategy = !treat_pinched_stars(umbrellas, treatment, pmesh);

  std::vector<typename GeomTraits::Point_3> points;
  std::vector<typename GeomTraits::Vector_3> normals;

  std::vector<Work_zone<PolygonMesh> > wzs = construct_work_zones(umbrellas, radius, pmesh, points, normals, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << wzs.size() << " work zones" << std::endl;
  static int i = 0;
  for(const Work_zone<PolygonMesh>& wz : wzs)
  {
    std::stringstream oss;
    oss << "results/zone_" << i++ << ".off" << std::ends;
    dump_cc(wz.faces, pmesh, oss.str().c_str());
  }

  std::ofstream("results/post_construction.off") << std::setprecision(17) << pmesh;
#endif

  if(treatment == SEPARATE)
    return fill_zones(wzs, pmesh, vpm, gt);

  if(!can_employ_merge_strategy || wzs.size() != 2)
  {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "Merging treatment requested, but configuration makes it impossible" << std::endl;
#endif
    return fill_zones(wzs, pmesh, vpm, gt);
  }
  else
  {
    const bool success = merge_zones(wzs.front(), wzs.back(), pmesh, points, normals, vpm, gt);
    if(!success)
    {
#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
      std::cout << "merging failed, falling back to separating strategy" << std::endl;
#endif
      return fill_zones(wzs, pmesh, vpm, gt);
    }

    return true;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Collect

// pretty much the same as 'PMP::non_manifold_vertices()', but consider the geometry instead of the combinatorics
// ({combinatorial non-manifold vertices} <= {geometrical non-manifold vertices}) so that creates a few changes
template <typename NMPContainer, typename NMVM,
          typename PolygonMesh, typename VPM, typename GeomTraits>
void geometrically_non_manifold_vertices(NMPContainer& nm_points, // m[vertex] = {halfedges}
                                         NMVM nm_marks,
                                         const PolygonMesh& pmesh,
                                         const VPM vpm,
                                         const GeomTraits& /*gt*/)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor          vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor        halfedge_descriptor;

  typedef typename boost::property_traits<VPM>::reference                       Point_ref;
  typedef typename boost::property_traits<VPM>::value_type                      Point;

  typedef CGAL::dynamic_halfedge_property_t<bool>                                       Halfedge_property_tag;
  typedef typename boost::property_map<PolygonMesh, Halfedge_property_tag>::const_type  Visited_halfedge_map;

  CGAL_precondition(nm_points.empty());

  Visited_halfedge_map visited_halfedges = get(Halfedge_property_tag(), pmesh);
  std::unordered_map<Point, halfedge_descriptor> visited_points;

  for(halfedge_descriptor h : halfedges(pmesh))
  {
    // If 'h' is not visited yet, we walk around the target of 'h' and mark these
    // halfedges as visited. Thus, if we are here and the target is already marked as visited,
    // it means that the vertex is non manifold.
    if(get(visited_halfedges, h))
      continue;

    put(visited_halfedges, h, true);

    bool is_non_manifold_due_to_multiple_umbrellas = false;
    const vertex_descriptor v = target(h, pmesh);
    const Point_ref p = get(vpm, v);

    const auto visited_itb = visited_points.emplace(p, h);
    if(!visited_itb.second) // already seen this point, but not from this star
    {
      is_non_manifold_due_to_multiple_umbrellas = true;
      put(nm_marks, v, true);

      // if this is the second time we visit that vertex and the first star was manifold, we have
      // not marked the vertex as non-manifold from the first star
      const auto nm_itb = nm_points.emplace(p, std::set<halfedge_descriptor>{h});
      if(nm_itb.second) // successful insertion
      {
        const halfedge_descriptor h_from_another_star = visited_itb.first->second;
        nm_itb.first->second.insert(h_from_another_star);
        put(nm_marks, target(h_from_another_star, pmesh), true);
      }
      else
      {
        nm_itb.first->second.insert(h);
      }
    }

    // While walking the star of this halfedge, if we meet a border halfedge more than once,
    // it means the mesh is pinched and we are also in the case of a non-manifold situation
    halfedge_descriptor ih = h, done = ih;
    int border_counter = 0;
    do
    {
      put(visited_halfedges, ih, true);
      if(is_border(ih, pmesh))
        ++border_counter;

      ih = prev(opposite(ih, pmesh), pmesh);
    }
    while(ih != done);

    if(border_counter > 1 && !is_non_manifold_due_to_multiple_umbrellas)
    {
      put(nm_marks, v, true);
      nm_points[p].insert(h); // might or might not have been an empty vector before
    }
  }
}

} // namespace internal

template <typename PolygonMesh, typename NamedParameters>
void treat_non_manifold_vertices(PolygonMesh& pmesh,
                                 const NM_TREATMENT treatment,
                                 const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor      halfedge_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type      VertexPointMap;
  VertexPointMap vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                        get_property_map(vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type          Geom_traits;
  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  typedef typename Geom_traits::FT                                            FT;
  typedef typename boost::property_traits<VertexPointMap>::value_type         Point;

  const FT radius = 0.06; // @todo automatic or np

  // Collect the non-manifold vertices
  typedef CGAL::dynamic_vertex_property_t<bool>                               Mark;
  typedef typename boost::property_map<PolygonMesh, Mark>::type               Marked_vertices;
  Marked_vertices nm_marks = get(Mark(), pmesh);
  for(vertex_descriptor v : vertices(pmesh))
    put(nm_marks, v, false);

  std::unordered_map<Point, std::set<halfedge_descriptor> > nm_points;
  internal::geometrically_non_manifold_vertices(nm_points, nm_marks, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
  std::cout << nm_points.size() << " cases to treat" << std::endl;
  for(auto& e : nm_points)
    std::cout << e.first << std::endl;
#endif

  internal::enforce_non_manifold_vertex_separation(nm_marks, pmesh, vpm, gt);
  std::ofstream("results/enforced_separation.off") << std::setprecision(17) << pmesh;

  for(auto& e : nm_points)
  {
    std::set<halfedge_descriptor>& umbrellas = e.second;
    CGAL_assertion(!umbrellas.empty());

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "NM vertex at pos " << get(vpm, target(*(std::begin(umbrellas)), pmesh))
              << " with " << umbrellas.size() << " incident umbrellas:" << std::endl;
    for(const halfedge_descriptor h : umbrellas)
      std::cout << h << " (" << get(vpm, source(h, pmesh)) << " " << get(vpm, target(h, pmesh)) << ")" << std::endl;
#endif

    internal::treat_umbrellas(umbrellas, treatment, radius, pmesh, vpm, gt);

#ifdef CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
    std::cout << "done with that nm vertex" << std::endl << std::endl;
    std::ofstream("results/intermediary.off") << std::setprecision(17) << pmesh;
#endif
  }

  CGAL_postcondition(is_valid_polygon_mesh(pmesh, true));

  std::cout << "done with all" << std::endl;
  CGAL_postcondition_code(nm_points.clear();)
  CGAL_postcondition_code(internal::geometrically_non_manifold_vertices(nm_points, nm_marks, pmesh, vpm, gt);)
  CGAL_postcondition(nm_points.empty());
}

template <typename PolygonMesh>
void treat_non_manifold_vertices(PolygonMesh& pmesh,
                                 const NM_TREATMENT treatment)
{
  return treat_non_manifold_vertices(pmesh, treatment, parameters::all_default());
}

template <typename PolygonMesh>
void treat_non_manifold_vertices(PolygonMesh& pmesh)
{
  return treat_non_manifold_vertices(pmesh, MERGE, parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_MANIFOLDNESS_H
