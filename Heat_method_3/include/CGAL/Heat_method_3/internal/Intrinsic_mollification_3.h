// Copyright (c) 2018  Carnegie Mellon University (USA), GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Christina Vaz, Keenan Crane, Andreas Fabri


#ifndef CGAL_INTRINSIC_MOLLIFICATION_3_H
#define CGAL_INTRINSIC_MOLLIFICATION_3_H

#include <CGAL/license/Heat_method_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/property_map.h>
#include <CGAL/double.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Heat_method_3/internal/V2V.h>


#include <boost/iterator/transform_iterator.hpp>

#include <unordered_map>
#include <set>
#include <stack>
#include <cmath>
#include <list>

#ifndef DOXYGEN_RUNNING

namespace CGAL {
namespace Heat_method_3 {

// forward declaration
template <typename IM>
struct IM_vertex_point_property_map;

template <typename TriangleMesh,
    typename Traits = typename Kernel_traits<typename boost::property_traits<typename boost::
            property_map<TriangleMesh, vertex_point_t>::const_type>::value_type>::Kernel>
struct Constant_mollification
{
  typedef Constant_mollification<TriangleMesh, Traits> Self;

  typedef boost::graph_traits<TriangleMesh>                        graph_traits;
  typedef typename graph_traits::vertex_descriptor            vertex_descriptor;
  typedef typename graph_traits::edge_descriptor                edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor        halfedge_descriptor;
  typedef typename graph_traits::face_descriptor                face_descriptor;
  typedef typename std::set<vertex_descriptor>::iterator        vertex_iterator;
  typedef typename std::set<edge_descriptor>::iterator            edge_iterator;
  /// Geometric typedefs
  typedef typename Traits::Point_3                                      Point_3;
  typedef typename Traits::FT                                                FT;
  typedef typename Traits::Vector_3                                    Vector_3;

  typedef int Index;

public:
  Constant_mollification(
      const TriangleMesh &input_tm, std::vector<double> &edge_lengths_, const double delta)
      : m_input_tm(input_tm), edge_lengths(edge_lengths_)
  {
    mollify(delta);
  }

private:
  // Mollification strategy to avoid degeneracies
  void
  mollify(const double delta)
  {
    // compute smallest length epsilon we can add to
    // all edges to ensure that the strict triangle
    // inequality holds with a tolerance of delta
    double epsilon = 0;
    for(halfedge_descriptor hd : halfedges(m_input_tm)) {
      halfedge_descriptor hd2 = next(hd, m_input_tm);
      halfedge_descriptor hd3 = next(hd2,m_input_tm);
      Index i = m_input_tm.edge_index(edge(hd, m_input_tm));
      Index j = m_input_tm.edge_index(edge(hd2, m_input_tm));
      Index k = m_input_tm.edge_index(edge(hd3, m_input_tm));
      double ineq = edge_lengths[j] + edge_lengths[k] - edge_lengths[i];
      epsilon = (std::max)(epsilon, (std::max)(0., delta-ineq));
    }
    // update edge lengths
    for(edge_descriptor ed : edges(m_input_tm)) {
        Index i = m_input_tm.edge_index(ed);
        edge_lengths[i] += epsilon;
    }
  }

  const TriangleMesh& m_input_tm; // this is the reference to the original
  TriangleMesh &edge_lengths;     // this is the reference to the original
};

/**
 * \ingroup PkgHeatMethodRef
 *
 * Class `Intrinsic_Delaunay_triangulation_3` is a remeshing algorithm to improve the approximation of the `Surface_mesh_geodesic_distances_3`.
 * It internally makes a copy of the triangle mesh, performs edge flips, and computes 2D vertex coordinates per face
 * which are stored in the halfedge with the vertex as target.
 *
 * The BGL API of this class .....
 *
 *
 * \tparam TriangleMesh a triangulated surface mesh, model of `FaceListGraph` and `HalfedgeListGraph`
 * \tparam Traits a model of `HeatMethodTraits_3`
 *
 * \cgalModels{FaceListGraph}
 */

template <typename TriangleMesh,
    typename Traits = typename Kernel_traits<typename boost::property_traits<typename boost::
            property_map<TriangleMesh, vertex_point_t>::const_type>::value_type>::Kernel,
          typename Scheme = Constant_mollification<TriangleMesh>>
class Intrinsic_mollification_3
{
  typedef Intrinsic_mollification_3<TriangleMesh, Traits, Scheme> Self;
  typedef Scheme Mollification_scheme;

  typedef boost::graph_traits<TriangleMesh>                        graph_traits;
  typedef typename graph_traits::vertex_descriptor            vertex_descriptor;
  typedef typename graph_traits::edge_descriptor                edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor        halfedge_descriptor;
  typedef typename graph_traits::face_descriptor                face_descriptor;
  typedef typename std::set<vertex_descriptor>::iterator        vertex_iterator;
  typedef typename std::set<edge_descriptor>::iterator            edge_iterator;
  /// Geometric typedefs
  typedef typename Traits::Point_3                                      Point_3;
  typedef typename Traits::FT                                                FT;
  typedef typename Traits::Vector_3                                    Vector_3;

  typedef std::pair<double,double>                                      Point_2;

  typedef int Index;

  typedef CGAL::dynamic_halfedge_property_t<Point_2> Halfedge_coordinate_tag;
  typedef typename boost::property_map<TriangleMesh, Halfedge_coordinate_tag >::type HalfedgeCoordinateMap;

  typedef typename boost::graph_traits<TriangleMesh>::vertices_size_type vertices_size_type;
  typedef typename boost::graph_traits<TriangleMesh>::edges_size_type edges_size_type;
  typedef typename boost::graph_traits<TriangleMesh>::faces_size_type faces_size_type;



private:
  friend struct IM_vertex_point_property_map<Self>;
  template <class IM, class VDM> friend struct IM_vertex_distance_property_map;

public: // for the BGL functions below. They should maybe become friend?
  typedef CGAL::Heat_method_3::IM_vertex_point_property_map<Self> Vertex_point_map;


  typedef typename TriangleMesh::Vertex_descriptor Vertex_descriptor;
  typedef typename TriangleMesh::Vertex_iterator_functor<TriangleMesh> Vertex_iterator_functor;

public:
  /// Constructor
  /// \param input_tm the triangle mesh
  Intrinsic_mollification_3(const TriangleMesh& input_tm)
    :m_input_tm(input_tm), hcm(get(Halfedge_coordinate_tag(), m_input_tm))
  {
    build();
  }

  template <class VertexPointMap>
  Intrinsic_mollification_3(const TriangleMesh& input_tm, VertexPointMap vpm)
    :m_input_tm(input_tm), hcm(get(Halfedge_coordinate_tag(), m_input_tm))
  {
    build(vpm);
  }

  typedef TriangleMesh Triangle_mesh;

  const Triangle_mesh&
  triangle_mesh() const
  {
    return m_input_tm;
  }


  Triangle_mesh&
  triangle_mesh()
  {
    return m_input_tm;
  }

  const HalfedgeCoordinateMap&
  hcmap() const
  {
    return hcm;
  }

private:
  template <class VertexPointMap>
  void
  build(VertexPointMap vpm)
  {
    CGAL_precondition(is_triangle_mesh(m_input_tm));

    typename Traits::Compute_squared_distance_3 squared_distance = Traits().compute_squared_distance_3_object();

    std::size_t number_of_edges = num_edges(m_input_tm);
    edge_lengths.resize(number_of_edges);
    Index edge_i = 0;

    double min_length = (std::numeric_limits<double>::max)();
    for(edge_descriptor ed : edges(m_input_tm)) {
      edge_lengths[edge_i] = CGAL::sqrt(to_double(
          squared_distance(get(vpm, source(ed, m_input_tm)), get(vpm, target(ed, m_input_tm)))));
      if (edge_lengths[edge_i] != 0 && edge_lengths[edge_i] < min_length) min_length = edge_lengths[edge_i++];
    }

    // mollify
    Mollification_scheme mollification_scheme(
        m_input_tm, edge_lengths, min_length * 1e-4);  // TODO: expose this parameter

    //now that edges are calculated, go through and for each face, calculate the vertex positions around it

    for(face_descriptor f : faces(m_input_tm)) {
      CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend, vmiddle;

      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f,m_input_tm),m_input_tm);
      halfedge_descriptor hd = halfedge(f,m_input_tm);
      if(face(hd,m_input_tm) != f) {
        hd = opposite(hd,m_input_tm);
      }
      hd = next(hd,m_input_tm);
      //each 'local' set of coordinates will have 0,0 at the first vertex/halfedge
      Point_2 p11(0,0);
      put(hcm, prev(hd,m_input_tm),p11);
      edge_descriptor ed1 = edge(hd, m_input_tm);
      hd = next(hd,m_input_tm);
      //the second local coordinate will be edge_length(first edge),0
      Point_2 p21(edge_lengths[m_input_tm.edge_index(ed1)], 0);
      put(hcm,prev(hd,m_input_tm),p21);

      //use basic trigonometry to compute third coordinate
      edge_descriptor ed2 = edge(hd, m_input_tm);
      hd = next(hd,m_input_tm);
      edge_descriptor ed3 = edge(hd, m_input_tm);
      Index e1 = m_input_tm.edge_index(ed1);
      Index e2 = m_input_tm.edge_index(ed2);
      Index e3 = m_input_tm.edge_index(ed3);
      double e1_len = edge_lengths[e1];
      double e2_len = edge_lengths[e2];
      double e3_len = edge_lengths[e3];
      double angle_a = -(e2_len*e2_len) + e3_len*e3_len + e1_len*e1_len;
      angle_a = acos(angle_a/(2*e3_len*e1_len));
      Point_2 p31(e3_len*std::cos(angle_a), e3_len*std::sin(angle_a));
      put(hcm,prev(hd,m_input_tm),p31);

    }
  }

  void
  build()
  {
    build( get(boost::vertex_point, m_input_tm) );
  }


  //todo:: determine which can be const
  const TriangleMesh& m_input_tm; // this is the reference to the original
  HalfedgeCoordinateMap hcm;

  std::vector<double> edge_lengths;
};

} // namespace Heat_method_3

} // namespace CGAL
namespace CGAL {
namespace Heat_method_3 {

template <typename IM>
struct IM_vertex_point_property_map {
  const IM& im;
  typedef typename IM::Triangle_mesh TM;
  typedef typename boost::graph_traits<IM>::vertex_descriptor key_type;
  typedef typename IM::Point_3 value_type;
  typedef typename IM::Point_2 Point_2;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  /**
   * Default constructor for vertex/point property map
   */
  IM_vertex_point_property_map(const IM& im)
    : im(im)
    {}


  /**
   * friend function for Heat method to get vertex descriptor's coordinates in iM's local coordinate system
   */
  friend value_type get(const IM_vertex_point_property_map<IM>& pm,
                        key_type vd)
  {
    const Point_2& p = get(pm.im.hcmap(), vd.hd);
    return value_type(p.first, p.second, 0);
  }
};

} // namespace Heat_method_3
} // namespace CGAL

#endif // DOXYGEN_RUNNING

#include <CGAL/enable_warnings.h>
#endif // CGAL_INTRINSIC_MOLLIFICATION_3_H
