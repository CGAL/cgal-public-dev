// Copyright (c) 2017, 2018 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H

#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph_node_dictionary.h>
#include <CGAL/Polyline_tracing/Motorcycle_priority_queue.h>
#include <CGAL/Polyline_tracing/Track.h>
#include <CGAL/Polyline_tracing/internal/Collision_information.h>
#include <CGAL/Polyline_tracing/internal/trace_validity_check.h>
#include <CGAL/Polyline_tracing/internal/Motorcycle_graph_builder.h>
#include <CGAL/Polyline_tracing/internal/robust_collinear.h>
#include <CGAL/Polyline_tracing/internal/robust_intersections.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/enum.h>
#include <CGAL/iterator.h>
#include <CGAL/number_utils.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/result_of.h>
#include <CGAL/use.h>

#include <boost/container/slist.hpp>
#include <boost/foreach.hpp> // @todo CGAL_foreach everywhere
#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/variant.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <limits>
#include <list>
#include <sstream>
#include <utility>
#include <vector>

namespace CGAL {

namespace Polyline_tracing {

// @todo handle degenerate faces in input
// @todo crash a motorcycle if new dest = current_pos ? (currently being done)
template<typename MotorcycleGraphTraits,
         typename MotorcycleType = Motorcycle<MotorcycleGraphTraits> >
class Motorcycle_graph
{
  typedef Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>           Self;

  typedef internal::Collision_information<Self>                             Collision_information;
  typedef typename MotorcycleGraphTraits::Kernel                            K; // needed for AABB_traits

public:
  typedef MotorcycleGraphTraits                                             Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                               Triangle_mesh;
  typedef typename Geom_traits::Face_graph                                  Face_graph;

  // Geometric types
  typedef typename Geom_traits::FT                                          FT;

  typedef typename Geom_traits::Point_2                                     Point_2;
  typedef typename Geom_traits::Segment_2                                   Segment_2;
  typedef typename Geom_traits::Vector_2                                    Vector_2;

  typedef typename Geom_traits::Point_d                                     Point;
  typedef typename Geom_traits::Segment_d                                   Segment;
  typedef typename Geom_traits::Vector_d                                    Vector;
  typedef typename Geom_traits::Ray_d                                       Ray;

  typedef typename Geom_traits::Bbox_d                                      Bbox;

  // Point types
  typedef Motorcycle_graph_node_dictionary<Geom_traits>                     Nodes;
  typedef typename Nodes::Node                                              Node;
  typedef typename Nodes::Node_ptr                                          Node_ptr;

  typedef typename Geom_traits::Barycentric_coordinates                     Barycentric_coordinates;
  typedef typename Geom_traits::Face_location                               Face_location;

  // Motorcycles
  typedef MotorcycleType                                                    Motorcycle;
  typedef std::deque<Motorcycle>                                            Motorcycle_container;
  typedef typename Motorcycle_container::iterator                           MCC_it;
  typedef typename Motorcycle_container::const_iterator                     MCC_cit;

  typedef Motorcycle_priority_queue<Self>                                   Motorcycle_PQ;
  typedef Motorcycle_priority_queue_entry<Self>                             Motorcycle_PQE;

  // Location-related types
  typedef boost::variant<Point, Face_location>                              Point_or_location;
  typedef boost::optional<Point_or_location>                                Optional_point_or_location;
  typedef boost::variant<Node_ptr, Face_location>                           Node_ptr_or_Face_location;

  // tuple of:
    // - #1: whether we have found a destination or not
    // - #2: the origin of the next path (might be different from mc.current_position() if on the border)
    // - #3: the destination of the next path
    // - #4: the time at the destination
    // - #5: whether the motorcycle should continue once at the destination
  // @todo change that to face locations? or variant of both?
  typedef boost::tuple<bool, Node_ptr, Node_ptr, FT, bool>                  Tracer_result;

  typedef Polygon_mesh_processing::internal::Point_to_Point_3<Triangle_mesh, Point>  Point_to_Point_3;
  typedef Polygon_mesh_processing::internal::Point_to_Point_3_VPM<Triangle_mesh>     AABB_tree_VPM;

  typedef CGAL::AABB_face_graph_triangle_primitive<Triangle_mesh, AABB_tree_VPM>     AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>                            AABB_face_graph_traits;
  typedef CGAL::AABB_tree<AABB_face_graph_traits>                                    AABB_tree;

  // BGL
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::edge_descriptor      edge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;
  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                   descriptor_variant;

  // Tracks and targets
  typedef typename Motorcycle::TPC_iterator                                 TPC_iterator;

  typedef Motorcycle_track<Geom_traits>                                     Track;
  typedef typename Track::Track_segment                                     Track_segment;
  typedef typename Track::iterator                                          Track_segment_ptr;
  typedef boost::container::slist<Track_segment_ptr>                        Track_segment_ptr_container;

  // To check collisions with adjacent faces, we need to know for a face the given track
  // track segments in this face. Pointers because tracks live in their respective
  // motorcycle and there is no need to duplicate information.
  typedef boost::unordered_map<face_descriptor,
                               Track_segment_ptr_container>                 Track_face_map;
  typedef typename Track_face_map::iterator                                 TFM_iterator;
  typedef typename Track_face_map::const_iterator                           TFM_const_iterator;

  typedef typename Collision_information::Foreign_collision_information     Foreign_collision_information;
  typedef typename Collision_information::Foreign_collisions_container      Foreign_collisions_container;
  typedef typename Foreign_collisions_container::const_iterator             FCC_cit;

  // Constructor mesh/graph building functions
  template<typename Face_graph_ptr>
  void initialize_graph(Face_graph_ptr graph) { graph_ = graph; }
  void initialize_graph(const boost::param_not_found);

  template<typename Triangle_mesh_ptr, typename NamedParameters>
  void initialize_mesh_and_graph(Triangle_mesh_ptr _mesh, const NamedParameters& np);
  template<typename NamedParameters>
  void initialize_mesh_and_graph(const boost::param_not_found, const NamedParameters& np);

  // Constructor
  template<typename NamedParameters>
  explicit Motorcycle_graph(const NamedParameters& np);

  // Destructor (needed because of pointers to mesh/graph)
  ~Motorcycle_graph();

  // Main function
  template<typename VertexNodeMap, typename EdgeTrackMap>
  void construct_motorcycle_graph(VertexNodeMap& vnmap,
                                  EdgeTrackMap& etmap,
                                  const bool construct_faces = false);

  void construct_motorcycle_graph(const bool construct_faces = false)
  {
    boost::dummy_property_map dummy;
    return construct_motorcycle_graph(dummy, dummy, construct_faces);
  }

  // Function to add motorcycles, forwards the arguments to the constructor of 'Motorcycle'
  template<typename ... Args>
  int add_motorcycle(const Args& ... args);

  void trace_graph();

  template<typename VertexNodeMap, typename EdgeTrackMap>
  void assemble_graph(VertexNodeMap& vnmap, EdgeTrackMap& etmap, const bool construct_faces = false);

  // Access
  const Geom_traits& geom_traits() const { return gt_; }
  Triangle_mesh& mesh() { return *mesh_; }
  const Triangle_mesh& mesh() const { return *mesh_; }
  AABB_tree& aabb_tree() { return aabb_tree_; }
  const AABB_tree& aabb_tree() const { return aabb_tree_; }
  Face_graph& graph() { return *graph_; }
  const Face_graph& graph() const { return *graph_; }
  Nodes& nodes() { return nodes_; }
  const Nodes& nodes() const { return nodes_; }
  Motorcycle_container& motorcycles() { return motorcycles_; }
  const Motorcycle_container& motorcycles() const { return motorcycles_; }

  FT latest_event_time() const { return latest_event_time_; }

  Motorcycle& motorcycle(const int id) {
    CGAL_precondition(id >= 0 && id < static_cast<int>(motorcycles_.size()));
    return motorcycles_[id];
  }
  const Motorcycle& motorcycle(const int id) const {
    CGAL_precondition(id >= 0 && id < static_cast<int>(motorcycles_.size()));
    return motorcycles_[id];
  }

  Motorcycle& motorcycle(MCC_it it) {
    CGAL_precondition(it >= motorcycles().begin() && it < motorcycles().end());
    return *it;
  }
  const Motorcycle& motorcycle(MCC_cit it) const {
    CGAL_precondition(it >= motorcycles().begin() && it < motorcycles().end());
    return *it;
  }

  std::size_t number_of_motorcycles() const { return motorcycles_.size(); }

  // Output
  void output_all_points() const;
  void output_tracks() const;
  void print_motorcycle_graph(std::ofstream& out) const;

private:
  /// \param fd face in which the segment belongs
  /// \param id the id of the motorcycle
  /// \param s, t the source and target of the oriented segment
  ///
  /// \return iterator in the tracking map
  TFM_iterator add_track_segment_to_map(face_descriptor fd, const Track_segment_ptr ts);

  void add_origin_node(Motorcycle& mc, const Point_or_location& input_origin);
  void add_destination_node(Motorcycle& mc, const Optional_point_or_location& input_destination);

  /// \param p, q first and second points
  /// \param p_time, q_time times at the first and second points
  ///
  /// \return new point and time at the new point
  std::pair<Node_ptr, FT> compute_halving_point(const Motorcycle& mc,
                                                Node_ptr p, const FT p_time,
                                                Node_ptr q, const FT q_time);

  /// \param p, q first and second points
  /// \param p_time, q_time times at the first and second points
  ///
  /// \return new point and time at the new point
  std::pair<Node_ptr, FT> compute_middle_point(Node_ptr p, const FT p_time,
                                               Node_ptr q, const FT q_time);

  bool compute_and_set_next_destination(Motorcycle& mc);
  bool initialize_next_path(Motorcycle& mc);
  void crash_motorcycle(Motorcycle& mc, const bool kill_regardless_of_stop_predicate = false); // @todo add visitors/callbacks ?
  void drive_to_closest_target(Motorcycle& mc);

  // Snapping functions
  std::pair<Node_ptr, bool> find_close_existing_point(const Face_location& location,
                                                      const Point& p,
                                                      const bool allow_same_location = true) const;
  void snap_to_existing_point(Node_ptr& e, const Node_ptr& t);

  typedef typename Collision_information::Collision_return         Collision_return;

  // Collisions between two motorcycles in different faces
  Collision_return find_collision_with_foreign_motorcycles(Motorcycle& mc, Collision_information& tc);

  // Below, only the target of the tentative track is on a border
  // ---------------------------------------------------------------------------------
  // collect the different faces in which we seek a collision depending on the location 'dv'
  Collision_return find_collision_with_tentative_track_extremity_on_border(const Motorcycle& mc,
                                                                           const Node_ptr extremity,
                                                                           const FT time_at_extremity,
                                                                           Collision_information& tc) const;
  // collect the motorcycles and tracks that we need to seek collisions with in the face 'ffd'
  Collision_return find_collision_with_tentative_track_extremity_on_border(const Motorcycle& mc,
                                                                           const Node_ptr extremity,
                                                                           const FT time_at_extremity,
                                                                           const face_descriptor ffd,
                                                                           Collision_information& tc) const;
  // triage based on the validity of 'fmc' and build the foreign track
  Collision_return find_collision_with_tentative_track_extremity_on_border_with_live_motorcycle_on_foreign_face(const Motorcycle& mc,
                                                                                                                const Node_ptr extremity,
                                                                                                                const FT time_at_extremity,
                                                                                                                const face_descriptor ffd,
                                                                                                                const Motorcycle& fmc,
                                                                                                                Collision_information& tc) const;
  // try to find a collision between the tentative tracks's target and the foreign track
  Collision_return find_collision_with_tentative_track_extremity_on_border_with_track_on_foreign_face(const Motorcycle& mc,
                                                                                                      const Node_ptr extremity,
                                                                                                      const FT time_at_extremity,
                                                                                                      const Track_segment& fmc_track,
                                                                                                      Collision_information& tc) const;
  // ---------------------------------------------------------------------------------

  // Below, both the source and the target of the tentative track are on the same halfedge
  // ---------------------------------------------------------------------------------
  // collect the motorcycles and tracks that we need to seek collisions with in the face 'opposite(hd, mesh)'
  Collision_return find_foreign_collision_with_tentative_track_on_border(const Motorcycle& mc,
                                                                         const halfedge_descriptor hd,
                                                                         Collision_information& tc);
  // discard 'fmc' if it's improper and build the foreign track
  Collision_return find_collision_with_live_motorcycle_on_foreign_face(const Motorcycle& mc,
                                                                       const halfedge_descriptor hd,
                                                                       const Motorcycle& fmc,
                                                                       Collision_information& tc) const;
  // distinguish between collinear tracks and foreign tracks with a single extremity on the halfedge
  Collision_return find_collision_with_track_on_foreign_face(const Motorcycle& mc,
                                                             const halfedge_descriptor hd,
                                                             const Track_segment& fmc_track,
                                                             const bool is_fmc_moving_on_track,
                                                             Collision_information& tc) const;
  // case of a foreign track and mc's tentative tracks being collinear
  Collision_return find_collision_with_collinear_tracks_on_different_faces(const Motorcycle& mc,
                                                                           const halfedge_descriptor hd,
                                                                           const Track_segment& fmc_track,
                                                                           const bool is_fmc_moving_on_track,
                                                                           Collision_information& tc) const;
  // case of a foreign track only having a single extremity on the halfedge 'opposite(hd, mesh)'
  Collision_return find_collision_with_foreign_track_extremity(const Motorcycle& mc,
                                                               const halfedge_descriptor hd,
                                                               const Motorcycle& fmc,
                                                               const Node_ptr foreign_extremity,
                                                               const FT foreign_time_at_collision,
                                                               Collision_information& tc) const;
  // ---------------------------------------------------------------------------------

  // Below, find collisions in a common face
  // ---------------------------------------------------------------------------------
  // collisions between two motorcycles in the same face
  Collision_return find_collision_at_tentative_track_destination(const Motorcycle& mc,
                                                                 const Motorcycle& fmc,
                                                                 const FT fmc_visiting_time,
                                                                 Collision_information& tc) const;
  Collision_return find_collision_between_collinear_tracks(const Motorcycle& mc,
                                                           const Segment_2& mcs,
                                                           const Motorcycle& fmc,
                                                           const Track_segment& fmc_track,
                                                           const Segment_2& fmcs,
                                                           const bool is_fmc_moving_on_track,
                                                           Collision_information& tc) const;

  Collision_return find_collision_between_tracks(const Motorcycle& mc,
                                                 const Track_segment& fmc_track,
                                                 const bool is_fmc_moving_on_track,
                                                 Collision_information& tc) const;
  Collision_return find_collision_with_complete_track(const Motorcycle& mc,
                                                      const Track_segment& fmc_track,
                                                      Collision_information& tc);
  Collision_return find_collision_with_live_motorcycle(Motorcycle& mc,
                                                       const Motorcycle& fmc,
                                                       Collision_information& tc);

  // \return collision point (if any), time at the collision for `mc`, id of
  //         the foreign intersecting motorcycle, time at the collision for
  //         the foreign motorcycle.
  Collision_return find_collision(Motorcycle& mc, Collision_information& tc);

  // \brief If in 2D and no mesh is passed, generate a large-enough triangle that encloses
  //        all the points and all the collisions.
  void generate_enclosing_face();

  bool has_simultaneous_collision(const Node_ptr node) const
  {
    CGAL_precondition(!node->visiting_motorcycles().empty());

    if(node->visiting_motorcycles().size() <= 1)
      return false;

    // Traverse the list of visiting motorcycles in order of visiting time
    // and check if the two earliest are visiting at (roughly) the same time.
    //
    // BUT: ignore the visiting time of a motorcycle if it is its starting time,
    //      because a motorcycle starting from a point at t0 should not block
    //      another motorcycle from going through that point at t0 (as long as
    //      they have different directions)

    FT earliest_visiting_time = std::numeric_limits<FT>::max();
    int number_of_motorcycles_at_earliest_time = 0; // ignoring motorcycles starting at node

    // Accessing 'right' of the bimap gives ordered time values.
    typename Node::VMC_right_cit mc_it = node->visiting_motorcycles().right.begin(),
                                 mc_end = node->visiting_motorcycles().right.end(),
                                 first = mc_it;
    for(; mc_it!=mc_end; ++mc_it)
    {
      const FT time = mc_it->first;

      const Motorcycle& mc = motorcycle(mc_it->second);
      const bool is_earliest_motorcycle = (mc_it == first);

      if(is_earliest_motorcycle)
        earliest_visiting_time = time;

      const bool motorcycle_can_create_simultaneous_collision = (!node->is_sibling(mc.origin()));
      if(!motorcycle_can_create_simultaneous_collision)
        continue;

      const bool is_at_earliest_time =
        (CGAL::abs(time - earliest_visiting_time) <= (tolerance_ * CGAL::abs(time + earliest_visiting_time)));

      if(is_at_earliest_time)
        ++number_of_motorcycles_at_earliest_time;
      else
        return false;

      if(number_of_motorcycles_at_earliest_time == 2)
        return true;
    }

    return false;
  }

  bool has_motorcycle_reached_crashing_point(const Motorcycle& mc) const;
  void initialize_motorcycle(Motorcycle& mc);

  void initialize_tracing();
  Face_location locate(const Point& p) const;

  void treat_collision(Motorcycle& mc, const Collision_information& collision);
  void treat_collision(Motorcycle& mc, Node_ptr collision_point, const FT time_at_collision);
  void treat_foreign_collision(Motorcycle& fmc, Node_ptr foreign_collision_point,
                               const FT foreign_time_at_collision, Track_segment_ptr foreign_track);

private:
  Geom_traits gt_;

  Nodes nodes_; // points that will be used throughout the algorithm
  Motorcycle_container motorcycles_;

  Motorcycle_PQ motorcycle_pq_; // motorcycle priority queue
  FT latest_event_time_; // useful to create multiple waves of motorcycles

  Triangle_mesh* mesh_; // input mesh
  Face_graph* graph_; // output graph
  bool is_mesh_provided, is_graph_provided;

  AABB_tree aabb_tree_;
  AABB_tree_VPM aabb_tree_vpm_;

  // map to store the completed tracks of the motorcycles for each face of the mesh
  Track_face_map track_face_map_;

  // snapping tolerance
  const FT tolerance_ = 1e-13;

private:
  // disable copy
  Motorcycle_graph& operator=(const Motorcycle_graph& other);
  Motorcycle_graph(const Motorcycle_graph& other);

  // disable move operators
  Motorcycle_graph(Motorcycle_graph&& other);
  Motorcycle_graph& operator=(Motorcycle_graph&& other);
};

// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
initialize_graph(const boost::param_not_found)
{
  graph_ = new Face_graph();
  is_graph_provided = false;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
template<typename Triangle_mesh_ptr, typename NamedParameters>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
initialize_mesh_and_graph(Triangle_mesh_ptr _mesh, const NamedParameters& np)
{
  mesh_ = _mesh;
  CGAL_precondition(num_vertices(mesh()) > 0);
  CGAL_precondition(CGAL::is_triangle_mesh(mesh()));

  initialize_graph(boost::get_param(np, CGAL::internal_np::output_graph));
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
template<typename NamedParameters>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
initialize_mesh_and_graph(const boost::param_not_found, const NamedParameters& np)
{
  // If no mesh is provided, we must be in 2D (technically, all points could be within a plane in 3D, but...)
  CGAL_precondition(Geom_traits::dimension() == 2);

  mesh_ = new Triangle_mesh();
  is_mesh_provided = false;

  //@tmp disabled while I find out what to do with the "no mesh provided option"
  // The issue is that the points are identified by a location described with barycentric
  // coordinates. I guess, I could generate a bbox, then a triangle that includes
  // the box ? Pretty ugly, though...
  CGAL_assertion(false);

  initialize_graph(boost::get_param(np, CGAL::internal_np::output_graph));
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
template <typename NamedParameters>
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
Motorcycle_graph(const NamedParameters& np)
  :
    gt_(),
    nodes_(),
    motorcycles_(),
    motorcycle_pq_(),
    latest_event_time_(),
    mesh_(),
    graph_(),
    is_mesh_provided(true),
    is_graph_provided(true),
    aabb_tree_(),
    aabb_tree_vpm_(),
    track_face_map_()
{
  using boost::choose_param;
  using boost::get_param;

  gt_ = choose_param(get_param(np, internal_np::geom_traits), Geom_traits());

  initialize_mesh_and_graph(get_param(np, CGAL::internal_np::input_mesh), np);

  CGAL_precondition(num_vertices(graph()) == 0);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
~Motorcycle_graph()
{
  if(!is_mesh_provided)
    delete mesh_;
  if(!is_graph_provided)
    delete graph_;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
template<typename ... Args>
int
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
add_motorcycle(const Args& ... args)
{
  int new_id = number_of_motorcycles();
  motorcycles_.emplace_back(args...);
  motorcycles_.back().set_id(new_id);

  return new_id;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::TFM_iterator
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
add_track_segment_to_map(face_descriptor fd, const Track_segment_ptr ts)
{
  CGAL_precondition(ts->source()->face() == fd);
  CGAL_precondition(ts->target()->face() == fd);
  CGAL_precondition(ts->motorcycle_id() >= 0 &&
                    ts->motorcycle_id() < static_cast<int>(number_of_motorcycles()));

  Track_segment_ptr_container l;
  l.push_front(ts);

  std::pair<TFM_iterator, bool> is_insert_success = track_face_map_.insert(std::make_pair(fd, l));

  if(!is_insert_success.second)
    is_insert_success.first->second.push_front(ts);

  return is_insert_success.first;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
add_origin_node(Motorcycle& mc, const Point_or_location& input_origin)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  Node_ptr origin;
  Face_location origin_location;
  Point origin_point;

  if(const Point* origin_point_ptr = boost::get<Point>(&input_origin))
  {
    origin_point = *origin_point_ptr;

    if(aabb_tree().empty())
    {
      aabb_tree_vpm_ = AABB_tree_VPM(mesh());
      PMP::build_AABB_tree(mesh(), aabb_tree(), parameters::vertex_point_map(aabb_tree_vpm_));
    }

    origin_location = locate(origin_point);
  }
  else
  {
    origin_location = boost::get<Face_location>(input_origin);
    origin_point = PMP::location_to_point(origin_location, mesh());
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  // handle nasty origin input
  bool snapped = PMP::internal::snap_location_to_border<Triangle_mesh>(origin_location, tolerance_);
  if(snapped)
    origin_point = PMP::location_to_point(origin_location, mesh());
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Origin point: " << origin_point << std::endl;
  std::cout << "Origin location: fd (" << origin_location.first
            << ") bc: [" << origin_location.second[0] << " "
                         << origin_location.second[1] << " "
                         << origin_location.second[2] << "]" << std::endl;
#endif

#define CGAL_MOTORCYCLE_GRAPH_SNAPPING_CODE
#ifdef CGAL_MOTORCYCLE_GRAPH_SNAPPING_CODE
  // Try to find an existing point close to that location
  std::pair<Node_ptr, bool> is_snappable = find_close_existing_point(origin_location, origin_point);
  if(is_snappable.second)
  {
    origin = is_snappable.first;
  }
  else
#endif
  {
    std::pair<Node_ptr, bool> is_insert_successful = nodes().insert(origin_location, origin_point, mesh());
    origin = is_insert_successful.first;
  }

  CGAL_postcondition(origin != Node_ptr());
  mc.origin() = origin;
  mc.current_position() = mc.origin();

  mc.origin()->add_motorcycle(mc.id(), mc.current_time());
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
add_destination_node(Motorcycle& mc,
                     const Optional_point_or_location& input_destination)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  // At the start of this function, mc.origin() must already be initialized
  CGAL_precondition(mc.origin() != Node_ptr());

  if(input_destination == boost::none) // A destination was not provided
  {
    compute_and_set_next_destination(mc);
    return;
  }

  // From now on, the destination is known, the time of arrival must be computed
  Node_ptr destination;
  FT time_at_destination;
  Face_location destination_location;
  Point destination_point;

  Face_location origin_location = mc.origin()->location();

  if(const Point* p = boost::get<Point>(&(*input_destination)))
  {
    destination_point = *p;

    // If the origin is on the border of the mesh, we must find a common face
    if(PMP::is_on_face_border(origin_location, mesh()))
    {
      PMP::locate_in_common_face(origin_location, destination_point, destination_location, mesh());
    }
    else // The origin is located strictly within a face
    {
      // Must ensure that the origin and destination are on the same face
      destination_location = PMP::locate_in_face(destination_point, origin_location.first, mesh());
    }
  }
  else // A 'Face_location' was provided in input
  {
    destination_location = boost::get<Face_location>(*input_destination);

    // The origin and destination must live in the same face
    if(origin_location.first != destination_location.first)
      PMP::locate_in_common_face(origin_location, destination_location, mesh());

    destination_point = PMP::location_to_point(destination_location, mesh());
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  // handle nasty origin input
  bool snapped = PMP::internal::snap_location_to_border<Triangle_mesh>(destination_location, tolerance_);
  if(snapped)
    destination_point = PMP::location_to_point(destination_location, mesh());
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Destination point: " << destination_point << std::endl;
  std::cout << "Destination location fd: (" << destination_location.first
            << ") bc: [" << destination_location.second[0] << " "
            << destination_location.second[1] << " "
            << destination_location.second[2] << "]" << std::endl;
#endif

  // 'origin_location' might have changed to find a common face
  if(origin_location != mc.origin()->location())
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cerr << "Warning: origin has changed!" << std::endl;
#endif
    mc.origin() = mc.origin()->sibling(origin_location.first);
    mc.current_position() = mc.origin();
    CGAL_assertion(mc.origin()->has_motorcycle(mc.id(), mc.time_at_origin()));
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_SNAPPING_CODE
  // Try to find an existing point close to that location
  std::pair<Node_ptr, bool> is_snappable = find_close_existing_point(destination_location, destination_point);

  if(is_snappable.second)
  {
    destination = is_snappable.first;
  }
  else
#endif
  {
    std::pair<Node_ptr, bool> is_insert_successful = nodes().insert(destination_location, destination_point, mesh());
    CGAL_assertion(is_insert_successful.second);
    destination = is_insert_successful.first;
  }

  CGAL_postcondition(destination != Node_ptr());
  mc.destination() = destination;

  const FT speed = mc.speed();
  const Point& origin_point = mc.origin()->point();
  time_at_destination = mc.time_at_origin() +
                        CGAL::sqrt(CGAL::squared_distance(origin_point, destination_point)) / speed;
  mc.time_at_destination() = time_at_destination;

  if(mc.origin() != mc.destination())
    mc.destination()->add_motorcycle(mc.id(), mc.time_at_destination());
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
std::pair<typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Node_ptr,
          typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::FT>
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
compute_halving_point(const Motorcycle& m, Node_ptr p, const FT p_time,
                                           Node_ptr q, const FT q_time)
{
  CGAL_USE(m);
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "***/***" << std::endl;
  std::cout << " Computing halving point on motorcycle #" << m.id() << "'s track."
            << " Points are:" << std::endl << "  " << *p << std::endl
                                           << "  " << *q << std::endl;
#endif
  CGAL_precondition(p != q);
  CGAL_precondition(p->face() == q->face());

#ifdef CGAL_MOTORCYCLE_GRAPH_USE_ADVANCED_HALVING_STRUCTURE
  // interface with the halving data structure @todo
#else
  return compute_middle_point(p, p_time, q, q_time);
#endif
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
std::pair<typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Node_ptr,
          typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::FT>
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
compute_middle_point(Node_ptr p, const FT p_time, Node_ptr q, const FT q_time)
{
  CGAL_assertion(p->face() == q->face());

  const Barycentric_coordinates& p_coords = p->location().second;
  const Barycentric_coordinates& q_coords = q->location().second;

  Barycentric_coordinates middle_coords = CGAL::make_array(0.5*(p_coords[0] + q_coords[0]),
                                                           0.5*(p_coords[1] + q_coords[1]),
                                                           0.5*(p_coords[2] + q_coords[2]));
  Face_location middle_loc = std::make_pair(p->face(), middle_coords);
  const Point middle_p = Polygon_mesh_processing::location_to_point(middle_loc, mesh());
  const FT time_at_r = 0.5 * (p_time + q_time);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  New middle point: (" << middle_p  << ") at time: " << time_at_r << std::endl;
  std::cout << "Location: " << p->face()
            << " bc: " << middle_coords[0] << " "
                       << middle_coords[1] << " "
                       << middle_coords[2] << std::endl;
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_SNAPPING_CODE
  std::pair<Node_ptr, bool> is_snappable = find_close_existing_point(middle_loc, middle_p);
  if(is_snappable.second)
  {
    // @todo should probably change the time here...
    return std::make_pair(is_snappable.first, time_at_r);
  }
  else
#endif
  {
    std::pair<Node_ptr, bool> entry = nodes().insert(middle_loc, mesh());
    return std::make_pair(entry.first, time_at_r);
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
compute_and_set_next_destination(Motorcycle& mc)
{
  Tracer_result res = mc.compute_next_destination(nodes(), mesh());

  if(!res.template get<0>()) // couldn't find a next path
  {
    mc.destination() = mc.current_position();
    mc.time_at_destination() = mc.current_time();

    return false;
  }

  const Node_ptr& next_origin = res.template get<1>();
  Node_ptr next_destination = res.template get<2>();
  const FT time_at_next_destination = res.template get<3>();
  const bool is_destination_final = res.template get<4>();

#ifdef CGAL_MOTORCYCLE_GRAPH_SNAPPING_CODE
  // 'false' because we want to ignore the location of 'next_destination'
  std::pair<Node_ptr, bool> is_snappable =
    find_close_existing_point(next_destination->location(), next_destination->point(), false);

  if(is_snappable.second && is_snappable.first != mc.current_position())
    snap_to_existing_point(next_destination, is_snappable.first);
#endif

  mc.destination() = next_destination;
  mc.time_at_destination() = time_at_next_destination;
  mc.is_destination_final() = is_destination_final;

  if(mc.destination() == mc.current_position())
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    if(!is_destination_final)
      std::cerr << "Warning: new destination is the current position but 'is_final' is set to 'no'!" << std::endl;
#endif

    // @todo allow this anyway? (for point set tracer or flow walking with degenerate faces)
    return false;
  }
  else
  {
    mc.destination()->add_motorcycle(mc.id(), mc.time_at_destination());
  }

  // The tracer might change the origin to ensure that the origin and the destination are on the same face
  if(mc.current_position() != next_origin)
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cerr << "Origin has changed!" << std::endl
              << "Previously: " << std::endl << *(mc.origin()) << std::endl
              << "Now: " << std::endl << *(res.template get<1>()) << std::endl;
#endif

    mc.current_position() = next_origin;
    mc.origin() = mc.current_position();
    CGAL_assertion(mc.current_position()->is_sibling(next_origin->location()));
  }

  return true;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
initialize_next_path(Motorcycle& mc)
{
  if(!compute_and_set_next_destination(mc))
     return false;

  mc.origin() = mc.current_position();
  mc.time_at_origin() = mc.current_time();

  // Add the next origin as target, even if it is equal to the current position.
  // This allows the new path to be treated with highest priority to compute
  // intersections on the new tentative track [origin-destination].
  mc.add_target(mc.origin(), mc.current_time());

  if(mc.origin() != mc.destination())
  {
    // No need to add the same information twice
    mc.add_target(mc.destination(), mc.time_at_destination());
  }

  CGAL_postcondition(mc.origin() != Node_ptr());
  CGAL_postcondition(mc.destination() != Node_ptr());
  CGAL_postcondition(mc.current_position() == mc.origin());
  CGAL_postcondition(mc.current_time() == mc.time_at_origin());
  CGAL_postcondition(mc.origin()->face() == mc.destination()->face());
  CGAL_postcondition(mc.has_target(mc.origin()).second);
  CGAL_postcondition(mc.has_target(mc.destination()).second);
  CGAL_postcondition(mc.origin()->has_motorcycle(mc.id(), mc.time_at_origin()));
  CGAL_postcondition(mc.destination()->has_motorcycle(mc.id(), mc.time_at_destination()));

  return true;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
crash_motorcycle(Motorcycle& mc, const bool crash_regardless_of_stop_predicate)
{
  if(mc.status() == Motorcycle::CRASHED)
  {
    std::cerr << "Warning: trying to crash an already crashed moto" << std::endl;
    return;
  }

  if(crash_regardless_of_stop_predicate || mc.stop_predicate()(mc.id(), *this))
  {
    mc.crash();
    motorcycle_pq_.erase(mc);
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
drive_to_closest_target(Motorcycle& mc)
{
  bool created_new_track_segment = mc.drive_to_closest_target();

  if(created_new_track_segment)
  {
    Track_segment_ptr ts_ptr = mc.newest_track_segment();
    CGAL_postcondition(mc.current_position() == ts_ptr->target());

    add_track_segment_to_map(mc.current_face(), ts_ptr);

    // @todo If we have just added a second track segment, we could remove
    // the first (degenerate) track segment (but a few code segments rely on
    // track.size() == 1 <=> we did not leave the origin, so gotta fix that too then).
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_foreign_motorcycles(Motorcycle& mc, Collision_information& tc)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "~~~~~~~~~X ?" << std::endl;
  std::cout << "Checking for collisions on motorcycle #" << mc.id() << "'s track with foreign faces" << std::endl;
#endif

  Collision_return res = Collision_information::NO_COLLISION;

  const bool is_source_on_face_border = PMP::is_on_face_border(mc.current_location(), mesh());
  const bool is_target_on_face_border = PMP::is_on_face_border(mc.closest_target()->location(), mesh());
  const bool check_for_collision_at_source = (is_source_on_face_border && !mc.has_left_starting_position());

  if(!is_target_on_face_border && !check_for_collision_at_source)
    return res;

  // If 'mc' is just leaving its starting position, we need to look for collisions at the source.
  if(check_for_collision_at_source)
  {
    res = res | find_collision_with_tentative_track_extremity_on_border(mc, mc.current_position(), mc.current_time(), tc);

    // Any collision found at the source is going to be the best collision
    if(res == Collision_information::COLLISION)
      return res;
  }

  if(is_target_on_face_border)
  {
    if(is_source_on_face_border)
    {
      // @todo make a function of out of below (same code exists in another function)
      // check if the source and target lie on the same halfedge of the border
      halfedge_descriptor hd = halfedge(mc.current_face(), mesh()), done(hd);
      bool are_on_same_halfedge = false;

      do
      {
        if(PMP::is_on_halfedge(mc.current_location(), hd, mesh()) &&
           PMP::is_on_halfedge(mc.closest_target()->location(), hd, mesh()))
        {
          are_on_same_halfedge = true;
          break;
        }

        hd = next(hd, mesh());
      }
      while(hd != done);
      // @todo till here (see above)

  #ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
      std::cout << "Tentative track on the same halfedge: " << are_on_same_halfedge << std::endl;
  #endif

      if(are_on_same_halfedge)
      {
        // We now look for collisions on the full tentative track
        // since "source" and "target" are both on the (same) border of the face.
        res = find_foreign_collision_with_tentative_track_on_border(mc, hd, tc);

        descriptor_variant target_dv = PMP::get_descriptor_from_location(mc.closest_target()->location(), mesh());
        if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&target_dv))
        {
          // if the target is on a vertex, we also need to check the faces incident to 'vd'
          res = res | find_collision_with_tentative_track_extremity_on_border(mc, mc.closest_target(), mc.time_at_closest_target(), tc);
        }

        return res;
      }
      else // both on the border, but not on the same halfedge, --> only look for collisions at the destination
      {
        return find_collision_with_tentative_track_extremity_on_border(mc, mc.closest_target(), mc.time_at_closest_target(), tc);
      }
    }
    else // is_target_on_face_border && !is_source_on_face_border
    {
      // Tentative track's source is not on a border --> only the target is on the border of the face

      // @todo generalize that, somewhere else
      // Small skip: if we have already found an intersection strictly within the face,
      // there's no point to check adjacent faces, since the intersection will be
      // at a later time.
      if(tc.time_at_closest_collision < mc.time_at_closest_target() -  tolerance_)
        return Collision_information::NO_COLLISION;

      return find_collision_with_tentative_track_extremity_on_border(mc, mc.closest_target(), mc.time_at_closest_target(), tc);
    }
  }

  return Collision_information::NO_COLLISION;
}

// Below, only the target of the tentative track is on a border
// ---------------------------------------------------------------------------------
template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_tentative_track_extremity_on_border(const Motorcycle& mc,
                                                        const Node_ptr extremity,
                                                        const FT time_at_extremity,
                                                        Collision_information& tc) const
{
  namespace PMP = CGAL::Polygon_mesh_processing;

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "¤ Find collision with tentative track extremity of motorcycle #" << mc.id() << " on border" << std::endl;
#endif

  descriptor_variant dv = PMP::get_descriptor_from_location(extremity->location(), mesh());
  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
  {
    const vertex_descriptor vd = *vd_ptr;

    // Check all faces incident to 'vd' and intersections at vd
    const halfedge_descriptor hd = halfedge(vd, mesh());
    BOOST_FOREACH(face_descriptor ffd, CGAL::faces_around_target(hd, mesh()))
    {
      if(ffd == mc.current_face() || ffd == boost::graph_traits<Triangle_mesh>::null_face())
        continue;

      return find_collision_with_tentative_track_extremity_on_border(mc, extremity, time_at_extremity, ffd, tc);
    }
  }
  else // the extremity is on a halfedge
  {
    const halfedge_descriptor hd = boost::get<halfedge_descriptor>(dv);

    if(is_border(edge(hd, mesh()), mesh()))
      return Collision_information::NO_COLLISION;

    // Check opposite face for a possible intersection at the extremity
    const face_descriptor ffd = face(opposite(hd, mesh()), mesh());
    return find_collision_with_tentative_track_extremity_on_border(mc, extremity, time_at_extremity, ffd, tc);
  }

  return Collision_information::NO_COLLISION;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_tentative_track_extremity_on_border(const Motorcycle& mc,
                                                        const Node_ptr extremity,
                                                        const FT time_at_extremity,
                                                        const face_descriptor ffd,
                                                        Collision_information& tc) const
{
  CGAL_precondition(ffd != boost::graph_traits<Triangle_mesh>::null_face());
  CGAL_precondition(mc.current_face() != ffd);

  Collision_return res = Collision_information::NO_COLLISION;

  // Step 1: check complete tracks
  TFM_const_iterator it = track_face_map_.find(ffd);
  if(it != track_face_map_.end())
  {
    const Track_segment_ptr_container& face_tracks = it->second;

    typename Track_segment_ptr_container::const_iterator ts_it = face_tracks.begin();
    typename Track_segment_ptr_container::const_iterator tl_end = face_tracks.end();
    for(; ts_it!=tl_end; ++ts_it)
    {
      const Track_segment& ts = *(*ts_it);
      Collision_return r = find_collision_with_tentative_track_extremity_on_border_with_track_on_foreign_face(mc, extremity, time_at_extremity, ts, tc);

      if(r == Collision_information::COLLISION) // Need to keep the foreign track in memory to add a new point on the confirmed track...
        tc.foreign_collisions.front().foreign_track_ptr = *ts_it;

      res = res | r;
    }
  }

  // Step 2: check incomplete tracks (path of a motorcycle currently moving in the same face)
  MCC_cit fmc_it = motorcycles().begin(), fmc_end = motorcycles().end();
  for(; fmc_it!=fmc_end; ++fmc_it)
  {
    const Motorcycle& fmc = motorcycle(fmc_it);
    res = res | find_collision_with_tentative_track_extremity_on_border_with_live_motorcycle_on_foreign_face(mc, extremity, time_at_extremity, ffd, fmc, tc);
  }

  return res;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_tentative_track_extremity_on_border_with_live_motorcycle_on_foreign_face(const Motorcycle& mc,
                                                                                             const Node_ptr extremity,
                                                                                             const FT time_at_extremity,
                                                                                             const face_descriptor ffd,
                                                                                             const Motorcycle& fmc,
                                                                                             Collision_information& tc) const
{
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "¤ Checking for foreign intersection at tentative track extremity "
            << "with live motorcycle #" << fmc.id() << std::endl;
#endif

  CGAL_precondition(ffd != boost::graph_traits<Triangle_mesh>::null_halfedge());
  CGAL_precondition(mc.current_face() != ffd);

  // uninitialized motorcycles have a degenerate track living in the future)
  if(!fmc.is_initialized())
    return Collision_information::NO_COLLISION;

  if(// the foreign motorcycle must be in the foreign face 'ffd'
     fmc.current_face() != ffd ||
     // the foreign motorcycle must be in motion
     fmc.status() == Motorcycle::CRASHED)
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << " ignoring 'fmc' in foreign face..." << std::endl;
    std::cout << "  > motorcycles #" << mc.id() << " and #" << fmc.id() << std::endl;
    std::cout << "  > faces: " << fmc.current_face() << " and " << fmc.current_face() << std::endl;
    std::cout << "  > crashed status: " << fmc.status() << std::endl;
#endif
    return Collision_information::NO_COLLISION;
  }

  CGAL_assertion(fmc.id() != mc.id());
  CGAL_assertion(!fmc.targets().empty());

  Track_segment fmc_track(fmc.id(),
                          fmc.current_position(), fmc.current_time(),
                          fmc.closest_target(), fmc.time_at_closest_target());

  return find_collision_with_tentative_track_extremity_on_border_with_track_on_foreign_face(mc, extremity, time_at_extremity, fmc_track, tc);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_tentative_track_extremity_on_border_with_track_on_foreign_face(const Motorcycle& mc,
                                                                                   const Node_ptr extremity,
                                                                                   const FT time_at_extremity,
                                                                                   const Track_segment& fmc_track,
                                                                                   Collision_information& tc) const
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  const int fmc_id = fmc_track.motorcycle_id();

  const Motorcycle& fmc = motorcycle(fmc_id);
  const Node_ptr fmc_track_source = fmc_track.source();
  const Node_ptr fmc_track_target = fmc_track.target();

  const bool is_fmcs_degenerate = (fmc_track_source == fmc_track_target);

  const face_descriptor ffd = fmc_track_source->face();
  CGAL_assertion(ffd == fmc_track_target->face());

  const Node_ptr extremity_in_ffd = extremity->sibling(ffd);
  CGAL_postcondition(ffd == extremity_in_ffd->face());

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "-----------------------------" << std::endl;
  std::cout << "¤¤ Checking collision with single point on border "
            << "with foreign motorcycle #" << fmc_id << std::endl;
  std::cout << " + extremity target: " << &*extremity << std::endl << *extremity << std::endl;
  std::cout << " + location in foreign face: " << " " << extremity_in_ffd->face() << " bc: "
                                               << extremity_in_ffd->barycentric_coordinate(0) << " "
                                               << extremity_in_ffd->barycentric_coordinate(1) << " "
                                               << extremity_in_ffd->barycentric_coordinate(2) << std::endl;
  std::cout << " + source: " << &*fmc_track_source << std::endl << *fmc_track_source << std::endl;
  std::cout << " + target: " << &*fmc_track_target << std::endl << *fmc_track_target << std::endl;
#endif

  const FT time_at_collision = time_at_extremity;
  const FT time_at_fmc_track_source = fmc_track.time_at_source();
  const FT time_at_fmc_track_target = fmc_track.time_at_target();

  FT foreign_visiting_time;
  if(extremity->has_motorcycle(fmc.id(), time_at_fmc_track_source, time_at_fmc_track_target, foreign_visiting_time))
  {
    if(extremity == mc.current_position())
      return Collision_information::NO_COLLISION; // nothing to do
    else // extremity == mc.closest_target()
      return tc.treat_potential_collision(extremity, time_at_collision, fmc.id(), foreign_visiting_time);
  }

  if(is_fmcs_degenerate)
  {
    // the only possible intersection is the extremity is the same point
    // as the degenerate foreign track
    if(extremity == fmc_track_source)
      return tc.treat_potential_collision(extremity, time_at_extremity, fmc.id(), time_at_fmc_track_source);
    else
      return Collision_information::NO_COLLISION;
  }

  descriptor_variant dv = PMP::get_descriptor_from_location(extremity->location(), mesh());
  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
  {
    // If the extremity is at a vertex, then the only possible intersection
    // is at 'fmc_track_source' or 'fmc_track_target' and it will (should)
    // have been detected with the checks above (if it exists).
    return Collision_information::NO_COLLISION;
  }
  else if(const halfedge_descriptor* hd_ptr = boost::get<halfedge_descriptor>(&dv))
  {
    halfedge_descriptor hd = *hd_ptr;

    // If the extremities of the foreign track are not on a border halfedge,
    // then there can't be an intersection with a point on the border (except
    // for the foreign track's source and target, which have been checked above)

    // Check if source and targets lie on the same halfedge
    halfedge_descriptor cfhd = halfedge(ffd, mesh()), done(cfhd);
    bool are_on_same_halfedge = false;

    do
    {
      if(PMP::is_on_halfedge(fmc_track_source->location(), cfhd, mesh()) &&
         PMP::is_on_halfedge(fmc_track_target->location(), cfhd, mesh()))
      {
        are_on_same_halfedge = true;
        break;
      }

      cfhd = next(cfhd, mesh());
    }
    while(cfhd != done);

    if(!are_on_same_halfedge)
      return Collision_information::NO_COLLISION;

    // 'hd' is in the non-foreign face, and we want the halfedge in the foreign face
    halfedge_descriptor opp_hd = opposite(hd, mesh());

    if(cfhd != opp_hd)
      return Collision_information::NO_COLLISION;

    // We are now in the configuration of 'mc' having an extremity on a halfedge,
    // and the foreign track is on the opposite halfedge
    const Point_2 s = geom_traits().construct_point_2_object()(fmc_track_source->location().second[0],
                                                               fmc_track_source->location().second[1]);
    const Point_2 t = geom_traits().construct_point_2_object()(fmc_track_target->location().second[0],
                                                               fmc_track_target->location().second[1]);
    const Point_2 ct2 = geom_traits().construct_point_2_object()(extremity_in_ffd->barycentric_coordinate(0),
                                                                 extremity_in_ffd->barycentric_coordinate(1));

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "s-ct2-t: " << s << " || " << ct2 << " || " << t << std::endl;
#endif

    CGAL_assertion(s != ct2 && t != ct2);

    // Below might fail due to numerical errors, but it is supposed to be 'true'
#ifdef CGAL_POLYLINE_TRACING_ENABLE_RIGOROUS_PRECONDITIONS
    CGAL_assertion(geom_traits().collinear_2_object()(s, ct2, t));
#endif

    // Check if the closest target is in between the source and the target
    if(!geom_traits().collinear_are_strictly_ordered_along_line_2_object()(s, ct2, t))
      return Collision_information::NO_COLLISION;

    // From here on, 'ct2' is strictly in between 's' and 't'

    // No choice but to compute the foreign time
    const FT time_at_fmc_track_source = fmc_track.time_at_source();
    const FT foreign_time_at_collision = time_at_fmc_track_source +
      CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(),
                                         extremity->point())) / fmc.speed();

    return tc.treat_potential_collision(extremity, time_at_extremity, fmc_id, foreign_time_at_collision);
  }
  else
  {
    // Motorcycle is strictly within a face and we shouldn't be here
    CGAL_assertion(false);
  }

  return Collision_information::NO_COLLISION;
}
// ---------------------------------------------------------------------------------

// Below, both the source and the target of the tentative track are on the same halfedge
// ---------------------------------------------------------------------------------
template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_foreign_collision_with_tentative_track_on_border(const Motorcycle& mc,
                                                      const halfedge_descriptor hd,
                                                      Collision_information& tc)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "¤ Checking collision with tentative track on border" << std::endl;
#endif

  const halfedge_descriptor opp_hd = opposite(hd, mesh());
  if(is_border(opp_hd, mesh()))
    return Collision_information::NO_COLLISION;

  Collision_return res = Collision_information::NO_COLLISION;
  const face_descriptor ffd = face(opp_hd, mesh());

  // Step 1: check complete tracks
  TFM_const_iterator it = track_face_map_.find(ffd);
  if(it != track_face_map_.end())
  {
    const Track_segment_ptr_container& face_tracks = it->second;

    typename Track_segment_ptr_container::const_iterator ts_it = face_tracks.begin();
    typename Track_segment_ptr_container::const_iterator ts_end = face_tracks.end();
    for(; ts_it!=ts_end; ++ts_it)
    {
      const Track_segment& ts = *(*ts_it);

      Collision_return r = find_collision_with_track_on_foreign_face(mc, hd, ts, false /*is_fmc_moving_on_track*/, tc);

      // Need to keep the foreign track in memory to add a new point on the confirmed track...
      if(r == Collision_information::COLLISION)
        tc.foreign_collisions.front().foreign_track_ptr = *ts_it;

      res = res | r;
    }
  }

  // Step 2: check incomplete tracks (path of a motorcycle currently moving in the same face)
  MCC_it fmc_it = motorcycles().begin(), fmc_end = motorcycles().end();
  for(; fmc_it!=fmc_end; ++fmc_it)
  {
    Motorcycle& fmc = motorcycle(fmc_it);
    res = res | find_collision_with_live_motorcycle_on_foreign_face(mc, hd, fmc, tc);
  }

  return res;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_live_motorcycle_on_foreign_face(const Motorcycle& mc,
                                                    const halfedge_descriptor hd,
                                                    const Motorcycle& fmc,
                                                    Collision_information& tc) const
{
  const face_descriptor ffd = face(opposite(hd, mesh()), mesh());

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "¤ Checking for foreign intersection with live motorcycle #" << fmc.id()
            << " in foreign face: " << ffd << std::endl;
#endif

  CGAL_precondition(!is_border(edge(hd, mesh()), mesh()));
  CGAL_precondition(mc.current_face() != ffd);

  // uninitialized motorcycles have a degenerate track living in the future)
  if(!fmc.is_initialized())
    return Collision_information::NO_COLLISION;

  if(// the foreign motorcycle must be in the foreign face 'ffd'
     fmc.current_face() != ffd ||
     // the foreign motorcycle must be in motion
     fmc.status() == Motorcycle::CRASHED)
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << " ignoring 'fmc' in foreign face..." << std::endl;
    std::cout << "  > motorcycles #" << mc.id() << " and #" << fmc.id() << std::endl;
    std::cout << "  > faces: " << mc.current_face() << " and " << fmc.current_face() << std::endl;
    std::cout << "  > crashed status: " << fmc.status() << std::endl;
#endif
    return Collision_information::NO_COLLISION;
  }

  CGAL_assertion(fmc.id() != mc.id());
  CGAL_assertion(!fmc.targets().empty());

  Track_segment fmc_track(fmc.id(), fmc.current_position(), fmc.current_time(),
                          fmc.closest_target(), fmc.time_at_closest_target());

  return find_collision_with_track_on_foreign_face(mc, hd, fmc_track, true /*is_fmc_moving_on_track*/, tc);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_track_on_foreign_face(const Motorcycle& mc,
                                          const halfedge_descriptor hd,
                                          const Track_segment& fmc_track,
                                          const bool is_fmc_moving_on_track,
                                          Collision_information& tc) const
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  const int fmc_id = fmc_track.motorcycle_id();

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "¤¤ Checking collision with tentative track on border "
            << "and foreign motorcycle #" << fmc_id << std::endl;
#endif

  CGAL_precondition(!is_border(edge(hd, mesh()), mesh()));

  const Motorcycle& fmc = motorcycle(fmc_id);
  const Node_ptr fmc_track_source = fmc_track.source();
  const FT time_at_fmc_track_source = fmc_track.time_at_source();
  const Node_ptr fmc_track_target = fmc_track.target();
  const FT time_at_fmc_track_target = fmc_track.time_at_target();

  const halfedge_descriptor opp_hd = opposite(hd, mesh());

  bool is_fts_on_opp_hd = PMP::is_on_halfedge(fmc_track_source->location(), opp_hd, mesh());
  bool is_ftd_on_opp_hd = PMP::is_on_halfedge(fmc_track_target->location(), opp_hd, mesh());

  if(is_fts_on_opp_hd)
  {
    if(is_ftd_on_opp_hd)
    {
      // foreign track is a subset (or the whole) of 'opp_hd'
      return find_collision_with_collinear_tracks_on_different_faces(mc, hd, fmc_track,
                                                                     is_fmc_moving_on_track, tc);
    }
    else // is_fts_on_opp_hd && !is_ftd_on_opp_hd
    {
      // the only possible intersection is at the source
      return find_collision_with_foreign_track_extremity(mc, hd, fmc, fmc_track_source,
                                                         time_at_fmc_track_source, tc);
    }
  }
  else if(is_ftd_on_opp_hd) // !is_fts_on_opp_hd && is_ftd_on_opp_hd
  {
    // only possible intersection is at the destination
    return find_collision_with_foreign_track_extremity(mc, hd, fmc, fmc_track_target,
                                                       time_at_fmc_track_target, tc);
  }

  return Collision_information::NO_COLLISION;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_collinear_tracks_on_different_faces(const Motorcycle& mc,
                                                        const halfedge_descriptor hd,
                                                        const Track_segment& fmc_track,
                                                        const bool is_fmc_moving_on_track,
                                                        Collision_information& tc) const
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  const int fmc_id = fmc_track.motorcycle_id();
  const Motorcycle& fmc = motorcycle(fmc_id);
  const Node_ptr fmc_track_source = fmc_track.source();
  const Node_ptr fmc_track_target = fmc_track.target();

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "¤¤¤ Find collision between collinear tracks on different faces of motorcycles #"
            << mc.id() << " and #" << fmc.id() << std::endl;
  std::cout << "   foreign track:" << std::endl << *fmc_track_source << std::endl
                                                << *fmc_track_target << std::endl;
#endif

  CGAL_precondition(PMP::is_on_halfedge(mc.current_location(), hd, mesh()));
  CGAL_precondition(PMP::is_on_halfedge(mc.closest_target()->location(), hd, mesh()));

  const halfedge_descriptor opp_hd = opposite(hd, mesh());
  CGAL_precondition(!is_border(opp_hd, mesh()));
  const face_descriptor ffd = face(opp_hd, mesh());

  // Convert the motorcycle track to locations on the foreign face
  const Face_location& cp_in_ffd = mc.current_position()->sibling(ffd)->location();
  const Face_location& ct_in_ffd = mc.closest_target()->sibling(ffd)->location();

  const Point_2 s = geom_traits().construct_point_2_object()(cp_in_ffd.second[0], cp_in_ffd.second[1]);
  const Point_2 t = geom_traits().construct_point_2_object()(ct_in_ffd.second[0], ct_in_ffd.second[1]);
  const Segment_2 mcs = geom_traits().construct_segment_2_object()(s, t);

  const Point_2 fs = geom_traits().construct_point_2_object()(fmc_track_source->location().second[0],
                                                              fmc_track_source->location().second[1]);
  const Point_2 ft = geom_traits().construct_point_2_object()(fmc_track_target->location().second[0],
                                                              fmc_track_target->location().second[1]);
  const Segment_2 fmcs = geom_traits().construct_segment_2_object()(fs, ft);

  return find_collision_between_collinear_tracks(mc, mcs, fmc, fmc_track, fmcs,
                                                 is_fmc_moving_on_track, tc);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_foreign_track_extremity(const Motorcycle& mc,
                                            const halfedge_descriptor hd,
                                            const Motorcycle& fmc,
                                            const Node_ptr foreign_extremity,
                                            const FT foreign_time_at_collision,
                                            Collision_information& tc) const
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  // this is the case of 'mc' tentative track being on a border, and a foreign
  // track with a single point on this same border

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "¤¤¤ Checking collision with tentative track on border"
            << " with foreign motorcycle " << fmc.id()
            << " and single foreign point on border: " << std::endl;
#endif

  // mc's track is non-degenerate
  CGAL_precondition(mc.current_position() != mc.closest_target());
  // mc's track in on the halfedge
  CGAL_precondition(PMP::is_on_halfedge(mc.current_location(), hd, mesh()));
  CGAL_precondition(PMP::is_on_halfedge(mc.closest_target()->location(), hd, mesh()));
  // the foreign extremity is on a halfedge
  CGAL_precondition(PMP::is_on_face_border(foreign_extremity->location(), mesh()));

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "foreign extremity: " << &*foreign_extremity
                                     << " (" << foreign_extremity->point() << ")" << std::endl;
#endif

  const halfedge_descriptor opp_hd = opposite(hd, mesh());
  CGAL_precondition(!is_border(opp_hd, mesh()));
  const face_descriptor ffd = face(opp_hd, mesh());
  CGAL_precondition(foreign_extremity->face() == ffd);

  const Face_location& cp_in_ffd = mc.current_position()->sibling(ffd)->location();
  const Face_location& ct_in_ffd = mc.closest_target()->sibling(ffd)->location();

  const Point_2 s = geom_traits().construct_point_2_object()(cp_in_ffd.second[0],
                                                             cp_in_ffd.second[1]);
  const Point_2 t = geom_traits().construct_point_2_object()(ct_in_ffd.second[0],
                                                             ct_in_ffd.second[1]);
  const Point_2 e = geom_traits().construct_point_2_object()(foreign_extremity->location().second[0],
                                                             foreign_extremity->location().second[1]);

  if(s == e) // intersection at mc's current_position
  {
    // ignore it, 'mc' would have been stopped before if that intersection was meaningful
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "s == e" << std::endl;
#endif
    return Collision_information::NO_COLLISION;
  }
  else if(t == e) // intersection at mc's closest target
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "t == e" << std::endl;
#endif

    return tc.treat_potential_collision(mc.closest_target(), mc.time_at_closest_target(), fmc.id(), foreign_time_at_collision);
  }
  else // general case
  {
    // The assertion below might fail due to numerical errors, but it is, logically,
    // a correct statement (case of three points on the same halfedge)
#ifdef CGAL_POLYLINE_TRACING_ENABLE_RIGOROUS_PRECONDITIONS
      CGAL_assertion(geom_traits().collinear_2_object()(s, e, t));
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
      std::cout << "Intersection not at source or target" << std::endl;
#endif

    if(!geom_traits().collinear_are_strictly_ordered_along_line_2_object()(s, e, t))
      return Collision_information::NO_COLLISION;

    // From here on, e is on ]s;t[
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "e is on ]s;t[" << std::endl;
#endif

    const Point collision_point = foreign_extremity->point();
    const FT time_at_collision = mc.current_time() +
      CGAL::sqrt(CGAL::squared_distance(mc.current_position()->point(),
                                        collision_point)) / mc.speed();

    return tc.treat_potential_collision(foreign_extremity, time_at_collision, fmc.id(), foreign_time_at_collision);
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_at_tentative_track_destination(const Motorcycle& mc,
                                              const Motorcycle& fmc,
                                              const FT fmc_visiting_time,
                                              Collision_information& tc) const
{
  FT time_at_collision = mc.time_at_closest_target();

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "  /!\\ Tentative path collides with track of motorcycle #" << fmc.id()
            << " at the closest target. Time: " << time_at_collision << std::endl;
#endif

  return tc.treat_potential_collision(mc.closest_target(), time_at_collision, fmc.id(), fmc_visiting_time);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_between_collinear_tracks(const Motorcycle& mc,
                                        const Segment_2& mcs,
                                        const Motorcycle& fmc,
                                        const Track_segment& fmc_track,
                                        const Segment_2& fmcs,
                                        const bool is_fmc_moving_on_track,
                                        Collision_information& tc) const
{
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "¤¤¤¤ Find collision between collinear tracks of motorcycles #"
            << mc.id() << " and #" << fmc.id() << std::endl;
#endif

#ifdef CGAL_POLYLINE_TRACING_ENABLE_RIGOROUS_PRECONDITIONS
  // Below might fail due to numerical errors, but we are treating here the case
  // of two collinear tracks, possibly on two different faces incident to the same edge.
  CGAL_precondition(geom_traits().collinear_2_object()(mcs.source(), fmcs.source(), mcs.target()));
  CGAL_precondition(geom_traits().collinear_2_object()(mcs.source(), fmcs.target(), mcs.target()));
#endif

  const Node_ptr fmc_track_source = fmc_track.source();
  const FT time_at_fmc_track_source = fmc_track.time_at_source();
  const Node_ptr fmc_track_target = fmc_track.target();
  const FT time_at_fmc_track_target = fmc_track.time_at_target();

  const face_descriptor ffd = fmc_track_source->face();
  const bool are_motorcycles_on_the_same_face = (mc.current_face() == ffd);

  // Some sanity checks
  if(!are_motorcycles_on_the_same_face)
  {
    // Check that all the track points are indeed on the same halfedge
    CGAL_precondition_code
    (
      boost::optional<halfedge_descriptor> hd =
        CGAL::Polygon_mesh_processing::common_halfedge(fmc_track_source->face(),
                                                       mc.current_face(),
                                                       mesh());
    )
    CGAL_precondition(bool(hd));
    CGAL_precondition_code(halfedge_descriptor opp_hd = opposite(*hd, mesh());)
    CGAL_precondition(CGAL::Polygon_mesh_processing::is_on_halfedge(fmc_track_source->location(), *hd, mesh()));
    CGAL_precondition(CGAL::Polygon_mesh_processing::is_on_halfedge(fmc_track_target->location(), *hd, mesh()));
    CGAL_precondition(CGAL::Polygon_mesh_processing::is_on_halfedge(mc.current_location(), opp_hd, mesh()));
    CGAL_precondition(CGAL::Polygon_mesh_processing::is_on_halfedge(mc.closest_target()->location(), opp_hd, mesh()));
  }
  // end of sanity checks -----

  // What happens below:
  // #1: Handle the very particular case of both segments having the same sources
  // #2: Handle the very particular case of 'mc' starting on a confirmed foreign track
  // #3: General case, split in two subcases whether the motorcycles are moving in the same direction or not

  const bool is_fmcs_degenerate = fmc_track.is_degenerate();

  // The respective direction of the two motorcycles
  bool are_motorcycles_moving_in_the_same_direction =
         (geom_traits().angle_2_object()(mcs.source(), mcs.target(),
                                         fmcs.source(), fmcs.target()) == CGAL::ACUTE);

  // #1: Handle the particular case of identical source points
  if(mcs.source() == fmcs.source())
  {
    if(is_fmcs_degenerate)
    {
      // Ignore the case of a degenerate fmc track starting at the same source as mc's
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
      std::cout << "degenerate fmc and mcs.source() == fmcs.source()" << std::endl;
#endif
      return Collision_information::NO_COLLISION;
    }

    if(!are_motorcycles_moving_in_the_same_direction)
      return Collision_information::NO_COLLISION;

    // The very specific case of two motorcycles starting from the same point,
    // in the same direction, and with the same time... (both motorcycles crash)
    if(is_fmc_moving_on_track)
    {
      std::cout << "Same origins, times, and directions!" << std::endl;

      CGAL_assertion(mc.current_time() == time_at_fmc_track_source);
      CGAL_assertion(mc.track().size() == 1);

      return tc.treat_potential_collision(mc.current_position(), mc.current_time(), fmc.id(), mc.current_time(),
                                          true /*crash mc*/, true /*crash fmc*/);
    }
  }

  // #2: Handle the case of a 'mc' starting within a (non-degenerate) confirmed foreign track.
  // This can happen when 'mc' is created on a pre-exisiting track
  if(!is_fmcs_degenerate && !is_fmc_moving_on_track)
  {
    if(mcs.source() == fmcs.source())
    {
      CGAL_assertion(are_motorcycles_moving_in_the_same_direction); // other case has been treated above

      return tc.treat_potential_collision(mc.current_position(), mc.current_time(),
                                          fmc.id(), time_at_fmc_track_source, true /*crash mc*/);
    }
    else if(mcs.source() == fmcs.target())
    {
      if(are_motorcycles_moving_in_the_same_direction)
      {
        return Collision_information::NO_COLLISION;
      }
      else
      {
        return tc.treat_potential_collision(mc.current_position(), mc.current_time(),
                                            fmc.id(), time_at_fmc_track_target, true /*crash mc*/);
      }
    }
    else if(geom_traits().collinear_are_strictly_ordered_along_line_2_object()(fmcs.source(),
                                                                               mcs.source(),
                                                                               fmcs.target()))
    {
      const FT foreign_time_at_collision = time_at_fmc_track_source +
        CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(),
                                          mc.current_position()->point())) / fmc.speed();

      return tc.treat_potential_collision(mc.current_position(), mc.current_time(),
                                          fmc.id(), foreign_time_at_collision, true /*crash mc*/);
    }
    else
    {
      return Collision_information::NO_COLLISION;
    }
  }

  // #3:
  // Many different configurations exist, e.g. (_S is for source, _T for target):
  //  MC_S  ---- FMC_S ---- FMC_T ---- MC_T
  //  FMC_T ---- MC_S  ---- FMC_S ---- MC_T
  // etc.
  // If, on the ray MC_S --> MC_T,
  // - FMC_S is "before" MC_S, then it doesn't matter for MC whichever respective
  //   direction the motorcycles are moving in.
  // - FMC_S is MC_S, then it only matters if they are moving in the same direction
  // - FMC_S is "after" MC_S, then it depends on the motorcycles' directions.

  // A degenerate track and a foreign motorcycle moving in the same direction as 'mc' is the same case
  are_motorcycles_moving_in_the_same_direction = (are_motorcycles_moving_in_the_same_direction || is_fmcs_degenerate);

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "  is degen: " << is_fmcs_degenerate << std::endl;
  std::cout << "  angle: " << geom_traits().angle_2_object()(mcs.source(), mcs.target(),
                                                             fmcs.source(), fmcs.target()) << std::endl;
  std::cout << "  are motorcycles moving in the same direction: "
            << are_motorcycles_moving_in_the_same_direction << std::endl;
#endif

  FT time_at_collision = 0.;

  // The motorcycles move in the same direction; if there's an intersection,
  // 'mc' will impact fmcs' source (except mc's source is on a confirmed track
  // and then it crashes).
  if(are_motorcycles_moving_in_the_same_direction)
  {
    if(mcs.target() == fmcs.source())
    {
      time_at_collision = mc.time_at_closest_target();
    }
    // Note that here, we know that fmcs.source() != mcs.source() and != mcs.target()
    else if(geom_traits().collinear_are_strictly_ordered_along_line_2_object()(mcs.source(),
                                                                               fmcs.source(),
                                                                               mcs.target()))
    {
      time_at_collision = mc.current_time() +
        CGAL::sqrt(CGAL::squared_distance(mc.current_position()->point(),
                                          fmc_track_source->point())) / mc.speed();
    }
    else
    {
      // fmcs.source() is either 'before' mcs.source() or 'after' mcs.target().
      // Either way, we don't care about any potential intersection.
      return Collision_information::NO_COLLISION;
    }

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "  Motorcycles #" << mc.id() << " crashes into the source of Motorcycle #"
                                   << fmc.id() << " at time: " << time_at_collision << std::endl;
#endif

    return tc.treat_potential_collision(fmc_track_source, time_at_collision, fmc.id(), time_at_fmc_track_source);
  }
  else // Motorcycles are moving in opposite directions
  {
    // Note that here we know that:
    // - fmcs is not degenerate
    // - mcs.source() != fmcs.source()

    // If the foreign source is 'before' mc's source, then there is no intersection
    if(mcs.target() != fmcs.source() && // to be able to use strictly (and there is an intersection if 'true')
       geom_traits().collinear_are_strictly_ordered_along_line_2_object()(fmcs.source(),
                                                                          mcs.source(),
                                                                          mcs.target()))
      return Collision_information::NO_COLLISION;

    // If mc's target is in [mcs, fmcs], then there is no intersection
    if(mcs.target() != fmcs.source() && // to be able to use strictly (and there is an intersection if 'true')
       mcs.target() != fmcs.target() &&
       geom_traits().collinear_are_strictly_ordered_along_line_2_object()(mcs.target(),
                                                                          fmcs.target(),
                                                                          fmcs.source()))
      return Collision_information::NO_COLLISION;

    // Now, we know that on the imaginary axis on which 'mc' is driving:
    // - fmcs is in ]mcs; infinity[
    // - fmct is in ]-infinity; mct]
    // - fmct is 'before' fmcs
    // Thus there is an intersection (except if fmcs = mcs, but we have already
    // discarded that case).
    // There are two cases to distinguish: moving 'fmc' and stationary 'fmc'.

    if(!is_fmc_moving_on_track) // stationary 'fmc'
    {
      // The foreign motorcycle is not moving on its track, thus 'mc' crashes
      // into the final position of the foreign track.

      // Check some known cases to avoid having to compute the collision time
      if(mcs.source() == fmcs.target())
      {
        time_at_collision = mc.current_time();
      }
      else if(mcs.target() == fmcs.target())
      {
        time_at_collision = mc.time_at_closest_target();
      }
      // Note that we know that fmcs.target() != mcs.source() and mcs.target()
      else if(geom_traits().collinear_are_strictly_ordered_along_line_2_object()(mcs.source(),
                                                                                 fmcs.target(),
                                                                                 mcs.target()))
      {
        // No choice but to compute the collision time
        time_at_collision = mc.current_time() +
          CGAL::sqrt(CGAL::squared_distance(mc.current_position()->point(),
                                            fmc_track_target->point())) / mc.speed();

        // @todo time snapping here ?

        CGAL_assertion(!mc.has_target_at_time(time_at_collision).second);
      }
      else
      {
        // fmcs.target() can't be 'before' mcs.source() because 'not_moving' means
        // that we are on a confirmed track and if fmcs.target() is 'after' mcs.target(),
        // then there is no intersection.
        return Collision_information::NO_COLLISION;
      }

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
      std::cout << "  Motorcycles #" << mc.id() << " crashes into the final position of Motorcycle #"
                << fmc.id() << " at time: " << time_at_collision << std::endl;
#endif

      return tc.treat_potential_collision(fmc_track_target, time_at_collision, fmc.id(), time_at_fmc_track_target);
    }
    else // The foreign motorcycle is (also) moving
    {
      // The collision is at the middle point and both motorcycles reach it at the same time.
      // Note that this point might not actually be reached by either motorcycle,
      // e.g. if a motorcycle crashes before reaching it.

      const FT sqd = CGAL::squared_distance(mc.current_position()->point(),
                                            fmc_track_source->point());
      time_at_collision = mc.current_time() +
        (CGAL::sqrt(sqd) - fmc.speed() * (mc.current_time() - time_at_fmc_track_source)) / (mc.speed() + fmc.speed());

      // @todo complete time snapping here
#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
      // By construction, the time and foreign_time should be greater
      // than the times at the sources of the tracks. Some numerical errors
      // can sneak it and, if so, fix the times.
      //
      // It should only be a numerical error, that is a very small error
      if(time_at_collision < mc.current_time())
      {
        CGAL_assertion(time_at_collision + tolerance_ >= mc.current_time()); // can't be too far off
        time_at_collision = mc.current_time();
      }

      if(time_at_collision < time_at_fmc_track_source)
      {
        CGAL_assertion(time_at_collision + tolerance_ >= time_at_fmc_track_source);
        time_at_collision = time_at_fmc_track_source;
      }
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
      std::cout << "  sqd: " << sqd << std::endl;
      std::cout << "  speeds: " << mc.speed() << " " << fmc.speed() << std::endl;
      std::cout << "  current times: " << mc.current_time() << " " << time_at_fmc_track_source << std::endl;
      std::cout << "  final time: " << time_at_collision << std::endl;
      std::cout << "  § mc and fmc would meet at time: " << time_at_collision << std::endl;
#endif

      CGAL_postcondition(time_at_collision >= time_at_fmc_track_source);
      CGAL_postcondition(time_at_collision >= mc.current_time());

      if(tc.compare_collision_time_to_closest(time_at_collision) != Collision_information::LATER_THAN_CURRENT_CLOSEST_TIME)
      {
        // both values are used later when we attempt to snap times/points
        const FT time_at_closest_collision_memory = tc.time_at_closest_collision;

        // Temporal snapping ---------------------------------------------------
        // Try to find the collision point by checking if any of the motorcycles
        // has a point at that time.
        std::pair<TPC_iterator, bool> mc_res = mc.has_target_at_time(time_at_collision);
        if(mc_res.second) // there is already a target at that time
        {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
          std::cout << "Motorcycle #" << mc.id() << " already has a target at time: " << time_at_collision << std::endl;
#endif

          TPC_iterator target_point = mc_res.first;
          CGAL_assertion(target_point->second == time_at_collision);
          Node_ptr alternate_collision = target_point->first;

          return tc.treat_potential_collision(alternate_collision, time_at_collision, fmc.id(), time_at_collision);
        }

        // Same check, but with the foreign time at collision
        std::pair<TPC_iterator, bool> fmc_res = fmc.has_target_at_time(time_at_collision);
        if(fmc_res.second) // there is already a target at that time
        {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
          std::cout << "Foreign motorcycle #" << fmc.id() << " already has a target at time: " << time_at_collision << std::endl;
#endif

          TPC_iterator target_point = fmc_res.first;
          Node_ptr alternate_foreign_collision = target_point->first;
          CGAL_assertion(alternate_foreign_collision->face() == fmc.current_face());
          CGAL_assertion(target_point->second == time_at_collision);

          return tc.treat_potential_collision(alternate_foreign_collision, time_at_collision, fmc.id(), time_at_collision);
        }

        // No choice but to construct the collision location
        const Vector_2 mcv(mcs);
        const FT ratio = (time_at_collision - mc.current_time()) /
                         (mc.time_at_closest_target() - mc.current_time());
        const Point_2 collision = mcs.source() + ratio * mcv;

        Face_location collision_location = std::make_pair(fmc_track_source->face(),
                                                          CGAL::make_array(collision[0], collision[1],
                                                                           1. - collision[0] - collision[1]));
#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
        // 1-x-y can result in some nasty "1e-17" imprecisions...
        CGAL::Polygon_mesh_processing::internal::snap_location_to_border<Triangle_mesh>(collision_location, tolerance_);
#endif

        // Couldn't find it through visiting times, but check if the new location
        // is already visited by 'mc' or 'fmc' (can happen due to numerical imprecisions)
        std::pair<Node_ptr, bool> collision_entry = nodes().find(collision_location);
        if(collision_entry.second) // the point already existed
        {
          CGAL_assertion(collision_entry.first != Node_ptr());

          // We previously searched by time but couldn't find anything but the
          // point existed. Check if that point is visited by either 'mc' or 'fmc';
          // if it's the case, we need to repare the time to be that of the existing
          // point.

          // Add a small tolerance on the time since we previously didn't find any target at the exact time
          FT visiting_time; // will be filled by the call to 'has_motorcycle'
          if(collision_entry.first->has_motorcycle(mc.id(), time_at_collision - tolerance_,
                                                   time_at_collision + tolerance_, visiting_time))
          {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
            std::cout << "Motorcycle #" << mc.id() << " already has a target at time: " << visiting_time << std::endl;
#endif

            // Assert that we are still the closest collision (not sure what to do otherwise)
            CGAL_assertion(visiting_time < time_at_closest_collision_memory);

            return tc.treat_potential_collision(collision_entry.first, visiting_time, fmc.id(), visiting_time);
          }

          // Try with 'fmc'
          FT foreign_visiting_time;
          if(collision_entry.first->has_motorcycle(fmc.id(), time_at_collision - tolerance_,
                                                   time_at_collision + tolerance_, foreign_visiting_time))
          {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
            std::cout << "Foreign motorcycle #" << fmc.id() << " already has a target at time: " << foreign_visiting_time << std::endl;
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
            std::cout << "found: fmc.id(): " << fmc.id() << " in pt: " << std::endl << *(collision_entry.first) << std::endl;
            std::cout << "foreign_visiting_time: " << foreign_visiting_time << std::endl;
#endif
            return tc.treat_potential_collision(collision_entry.first, foreign_visiting_time, fmc.id(), foreign_visiting_time);
          }
        }
        else
        {
          // At this point, we have a new location at an unknown time...
#ifdef CGAL_MOTORCYCLE_GRAPH_SNAPPING_CODE
          // But maybe there exists another point that is very close! Check for it,
          // and if needed, snap the new location (and the time) to it.

          Point collision_point = CGAL::Polygon_mesh_processing::location_to_point(collision_location, mesh());

          std::pair<Node_ptr, bool> is_snappable = find_close_existing_point(collision_location, collision_point);
          if(is_snappable.second) // successful snapping
          {
            FT visiting_time = time_at_collision;

            // the call to this function will modify 'visiting_time' if the point of snapping is already is visited by 'mc'
            const FT min_visiting_time = time_at_collision - tolerance_;
            const FT max_visiting_time = time_at_collision + tolerance_;
            if(!is_snappable.first->has_motorcycle(mc.id(), min_visiting_time, max_visiting_time, visiting_time))
            {
              // While trying to get the visiting time, if the snapped point is
              // not visited by 'mc', check if it is visited by 'fmc'
              is_snappable.first->has_motorcycle(fmc.id(), min_visiting_time, max_visiting_time, visiting_time);
            }

            // We have snapped so we are igoring times that we had set up as best, but
            // we need to make sure it is still better than the previous one.
            CGAL_assertion(visiting_time <= time_at_closest_collision_memory);

            return tc.treat_potential_collision(is_snappable.first, visiting_time, fmc.id(), visiting_time);
          }
#endif

          // Couldn't snap to anything, the 'collision_location' is definitely a new point
          return tc.treat_potential_collision(collision_location, time_at_collision, fmc.id(), time_at_collision);
        }
      }
    }
  }

  return Collision_information::NO_COLLISION;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_between_tracks(const Motorcycle& mc,
                              const Track_segment& fmc_track,
                              const bool is_fmc_moving_on_track,
                              Collision_information& tc) const
{
  // Use the barycentric coordinate system to compute intersections
  const Point_2 s = geom_traits().construct_point_2_object()(mc.current_location().second[0],
                                                             mc.current_location().second[1]);
  const Point_2 t = geom_traits().construct_point_2_object()(mc.closest_target()->location().second[0],
                                                             mc.closest_target()->location().second[1]);
  const Segment_2 mcs = geom_traits().construct_segment_2_object()(s, t);

  const Motorcycle& fmc = motorcycle(fmc_track.motorcycle_id());
  const Node_ptr fmc_track_source = fmc_track.source();
  const FT time_at_fmc_track_source = fmc_track.time_at_source();
  const Node_ptr fmc_track_target = fmc_track.target();
  const FT time_at_fmc_track_target = fmc_track.time_at_target();

  const Point_2 fs = geom_traits().construct_point_2_object()(fmc_track_source->location().second[0],
                                                              fmc_track_source->location().second[1]);
  const Point_2 ft = geom_traits().construct_point_2_object()(fmc_track_target->location().second[0],
                                                              fmc_track_target->location().second[1]);
  const Segment_2 fmcs = geom_traits().construct_segment_2_object()(fs, ft);

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "¤¤ Checking collision with track of motorcycle #" << fmc.id() << std::endl;
  std::cout << " + source: " << &*fmc_track_source << std::endl << *fmc_track_source << std::endl;
  std::cout << " + target: " << &*fmc_track_target << std::endl << *fmc_track_target << std::endl;
#endif

  CGAL_assertion(mc.current_face() == fmc_track_source->face());

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  if(internal::are_logically_collinear_on_border<Geom_traits>(
                 mc.current_location(), mc.closest_target()->location(),
                 fmc_track_source->location(), fmc_track_target->location()))
  {
    return find_collision_between_collinear_tracks(mc, mcs, fmc, fmc_track, fmcs,
                                                   is_fmc_moving_on_track, tc);
  }
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  bool are_tracks_collinear = internal::are_segments_collinear<Geom_traits>(mcs, fmcs);
#else
  bool are_tracks_collinear = (geom_traits().collinear_2_object()(mcs.source(), mcs.target(), fmcs.source()) &&
                               geom_traits().collinear_2_object()(mcs.source(), mcs.target(), fmcs.target()));
#endif

  // Detect whether the motorcycles share the same supporting line.
  // @todo should this have a tolerance ?
  if(are_tracks_collinear)
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "  /!\\ Tracks are aligned" << std::endl;
#endif
    return find_collision_between_collinear_tracks(mc, mcs, fmc, fmc_track, fmcs,
                                                   is_fmc_moving_on_track, tc);
  }

  // --- From here on, the tracks are not collinear ---
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "Tracks are not collinear" << std::endl;
#endif

  // If either track is degenerate, a collision can only happen if the tracks are collinear
  bool is_fmcs_degenerate = fmc_track.is_degenerate();
  if(is_fmcs_degenerate || mc.is_tentative_track_degenerate())
    return Collision_information::NO_COLLISION;

  // @todo move all the checks below to their own functions

  // Below are a bunch of checks to branch out easily without computing an explicit
  // collision time/position. Note that the tracks are _not_ collinear.
  //
  // #1: Check if the current position of mc is a known intersection with the foreign track
  // #2: Check if the closest target of mc is a known intersection with the foreign track
  // #3/#4: Robustness for intersections on halfedge

  // Check #1: known collision at current_position
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "  Check #1: motorcycle #" << fmc.id() << " between "
            << time_at_fmc_track_source << " " << time_at_fmc_track_target << std::endl;
#endif
  if(mc.current_position()->has_motorcycle(fmc.id(), time_at_fmc_track_source, time_at_fmc_track_target))
  {
    // Ignore this intersection: we are seeking collisions in the tentative track,
    // it means that the position was not blocked
    return Collision_information::NO_COLLISION;
  }

  // Check #2: known collision at closest_target
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "  Check #2: collisition at tentative's track destination ?" << std::endl;
#endif
  FT foreign_visiting_time; // will be filled by 'has_motorcycle' if fmc visits 'closest_target'
  if(mc.closest_target()->has_motorcycle(fmc.id(), time_at_fmc_track_source,
                                         time_at_fmc_track_target, foreign_visiting_time))
  {
    return find_collision_at_tentative_track_destination(mc, fmc, foreign_visiting_time, tc);
  }

  CGAL_assertion(fmc_track_source != mc.current_position());
  CGAL_assertion(fmc_track_source != mc.closest_target());
  CGAL_assertion(fmc_track_target != mc.current_position());
  CGAL_assertion(fmc_track_target != mc.closest_target());

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  // Check #3: collision at the closest target, with foreign track on a halfedge
  //
  // Catch some annoying numerical issue: the configuration of FMCS on a halfedge
  // and the motorcycle destination on the same edge (but do_intersect_2()
  // does not find it due to numerical issues...).
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "  check #3: foreign track and target on the same border" << std::endl;
#endif
  if(internal::are_logically_collinear_on_border<Geom_traits>(
       fmc_track_source->location(), mc.closest_target()->location(), fmc_track_target->location()))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "  foreign track and target are logically collinear on border" << std::endl;
#endif

    if(geom_traits().collinear_are_strictly_ordered_along_line_2_object()(
         fmcs.source(), mcs.target(), fmcs.target()))
    {
      const FT time_at_collision = mc.time_at_closest_target();
      const FT foreign_time_at_collision = time_at_fmc_track_source +
        CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(),
                                          mc.closest_target()->point())) / fmc.speed();

      return tc.treat_potential_collision(mc.closest_target(), time_at_collision, fmc.id(), foreign_time_at_collision);
    }

    return Collision_information::NO_COLLISION;
  }

  // Check #3bis: collision at tentative track's source, with foreign track on a halfedge.
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "  check #3bis: foreign track and source on the same border" << std::endl;
#endif
  if(internal::are_logically_collinear_on_border<Geom_traits>(
       fmc_track_source->location(), mc.current_location(), fmc_track_target->location()))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "  foreign track and source are logically collinear on border" << std::endl;
#endif

    if(geom_traits().collinear_are_strictly_ordered_along_line_2_object()(
         fmcs.source(), mcs.source(), fmcs.target()))
    {
      const FT time_at_collision = mc.current_time();
      const FT foreign_time_at_collision = time_at_fmc_track_source +
        CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(),
                                          mc.current_position()->point())) / fmc.speed();

      return tc.treat_potential_collision(mc.current_position(), time_at_collision, fmc.id(), foreign_time_at_collision);
    }

    return Collision_information::NO_COLLISION;
  }

  // Check #4: collision at foreign track's target, with tentative track on a halfedge.
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "  check #4: tentative track and foreign destination on the same border" << std::endl;
#endif
  if(internal::are_logically_collinear_on_border<Geom_traits>(
      fmc_track_target->location(), mc.closest_target()->location(), mc.current_location()))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "  track and foreign target are logically collinear on border" << std::endl;
#endif

    if(geom_traits().collinear_are_strictly_ordered_along_line_2_object()(
         mcs.source(), fmcs.target(), mcs.target()))
    {
      const FT sqd = CGAL::squared_distance(mc.current_position()->point(),
                                            fmc_track_target->point());
      const FT time_at_collision = mc.current_time() + CGAL::sqrt(sqd) / mc.speed();
      const FT foreign_time_at_collision = time_at_fmc_track_target;

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
      std::cout << "  foreign target in ] track [ " << std::endl;
      std::cout << "  Pts: (" << mc.current_position()->point() << ") -- ("
                              << fmc_track_target->point() << ")" << std::endl;
      std::cout << "  current time: " << mc.current_time() << std::endl;
      std::cout << "  sqd: " << sqd << std::endl;
      std::cout << "  time at collision: " << time_at_collision << std::endl;
#endif

      return tc.treat_potential_collision(fmc_track_target, time_at_collision, fmc.id(), foreign_time_at_collision);
    }

    return Collision_information::NO_COLLISION;
  }

  // Check #4bis: collision at foreign track's source, with tentative track on a halfedge.
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "  check #4: tentative track and foreign source on the same border" << std::endl;
#endif
  if(internal::are_logically_collinear_on_border<Geom_traits>(
      fmc_track_source->location(), mc.closest_target()->location(), mc.current_location()))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "  track and foreign source are logically collinear on border" << std::endl;
#endif

    if(geom_traits().collinear_are_strictly_ordered_along_line_2_object()(
         mcs.source(), fmcs.source(), mcs.target()))
    {
      const FT time_at_collision = mc.current_time() +
        CGAL::sqrt(CGAL::squared_distance(mc.current_position()->point(),
                                          fmc_track_source->point())) / mc.speed();
      const FT foreign_time_at_collision = time_at_fmc_track_source;

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
      std::cout << "  foreign source in ] track [, "
                << "time at collision: " << time_at_collision << std::endl;
#endif

      return tc.treat_potential_collision(fmc_track_source, time_at_collision, fmc.id(), foreign_time_at_collision);
    }

    return Collision_information::NO_COLLISION;
  }
#endif

  // --- The general-est case: the intersection must be computed ---
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "  general case..." << std::endl;
#endif

  if(!geom_traits().do_intersect_2_object()(mcs, fmcs))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "  No intersection (general case)" << std::endl;
#endif
    return Collision_information::NO_COLLISION;
  }

  // Here is what happens next:
  // 1. we compute the intersection
  // 2. we check if this new location is (EXACTLY) an existing node. If it is, we compute
  //    the visiting times and return the appropriate result.
  // 3. if not an existing node, we compute visiting times and check if there are already existing
  //    nodes at these times on the trajectories of 'mc' and 'fmc'. If there is, we use that position.
  //    to the position and return the appropriate result.
  // 4. we check if there is an existing node that is close to the collision point.
  //    If there is, we snap to that position, compute the visiting times and return
  //    the appropriate result.
  // 5. return the new collision point.

  // Step 1: compute the intersection in the barycentric coordinates system
  Point_2 collision = internal::robust_intersection<Geom_traits>(mcs, fmcs, geom_traits());

  // Convert it to a location in the ambiant dimension
  Barycentric_coordinates coords = CGAL::make_array(collision[0], collision[1],
                                                    1. - collision[0] - collision[1]);
  Face_location collision_location = std::make_pair(mc.current_face(), coords);

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  // 1-x-y can result in some nasty "1e-17" imprecisions...
  CGAL::Polygon_mesh_processing::internal::snap_location_to_border<Triangle_mesh>(collision_location, tolerance_);
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "  /!\\ collision between motorcycles #" << mc.id() << " and #" << fmc.id() << std::endl;
  std::cout << "Collision location: " << collision_location.first << " bc: "
            << collision_location.second[0] << " " << collision_location.second[1] << " " << collision_location.second[2] << std::endl;
#endif

  // Step 2: although we might not have known that these two tracks do intersect,
  // their intersection might be a point that has been used previously
  std::pair<Node_ptr, bool> is_already_in_dictionary = nodes().find(collision_location);
  if(is_already_in_dictionary.second)
  {
    Node_ptr collision_point = is_already_in_dictionary.first;
    FT time_at_collision = 0.;

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "Already in the dictionary at: " << &*collision_point << std::endl << *collision_point << std::endl;
#endif

    // Check if 'mc' already visits the known collision point
    if(collision_point == mc.current_position())
    {
      time_at_collision = mc.current_time();
    }
    else if(collision_point == mc.closest_target())
    {
      time_at_collision = mc.time_at_closest_target();
    }
    else // 'collision_point' is a known point but has not (yet) been visited by 'mc'
    {
      // The tentative track of 'mc' can only be intersected at a known point that has 'mc'
      // if that known point is the current position or the closest target.
      CGAL_assertion(!collision_point->has_motorcycle(mc.id(), mc.current_time(), mc.time_at_closest_target()));

      // No choice but to compute the visiting time
      time_at_collision = mc.current_time() +
        CGAL::sqrt(CGAL::squared_distance(mc.current_position()->point(),
                                          collision_point->point())) / mc.speed();

#ifdef CGAL_MOTORCYCLE_GRAPH_SNAPPING_CODE
      // Although we have found an _existing_ point at the location of the intersection,
      // this point was neither the source or the closest target of 'mc'.
      // Global snapping makes sure that points are not too close from one another.
      // Consequently, the times should be different.
      CGAL_assertion(time_at_collision != mc.current_time());
      CGAL_assertion(time_at_collision != mc.time_at_closest_target());
#endif
    }

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "time_at_collision: " << time_at_collision
              << " (closest is: " << tc.time_at_closest_collision << ") " << std::endl;
#endif

    // Partial test of "is_collision_earlier..." to branch out early
    if(time_at_collision <= tc.time_at_closest_collision)
    {
      // Check if 'fmc' already visits the known collision point
      FT foreign_time_at_collision;
      if(collision_point->has_motorcycle(fmc.id(), time_at_fmc_track_source,
                                         time_at_fmc_track_target, foreign_time_at_collision))
      {
        // The collision point is visited by 'fmc' at time 'foreign_time_at_collision'
      }
      else // 'collision_point' is a known point but has not (yet) been visited by 'fmc'
      {
        // No choice but to compute the foreign visiting time
        const FT sqd = CGAL::squared_distance(fmc_track_source->point(),
                                              collision_point->point());
        foreign_time_at_collision = time_at_fmc_track_source + CGAL::sqrt(sqd) / fmc.speed();

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
      std::cout << "  Gotta compute the foreign time " << std::endl;
      std::cout << "  Pts: (" << fmc_track_source->point() << ") -- ("
                              << collision_point->point() << ")" << std::endl;
      std::cout << "  foreign source time: " << time_at_fmc_track_source << std::endl;
      std::cout << "  sqd: " << sqd << std::endl;
      std::cout << "  foreign time at collision: " << foreign_time_at_collision << std::endl;
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_SNAPPING_CODE
        // Although we have found an _existing_ point at the location of the intersection,
        // this point was neither the source or the closest target of 'mc'.
        // Global snapping makes sure that points are not too close from one another.
        // Consequently, the times should be different.
        CGAL_assertion(!fmc.has_target_at_time(foreign_time_at_collision).second);
#endif
      }

      return tc.treat_potential_collision(collision_point, time_at_collision, fmc.id(), foreign_time_at_collision);
    }
  }
  else // The collision location has never been seen before!
  {
    Point collision_point = CGAL::Polygon_mesh_processing::location_to_point(collision_location, mesh());

    FT time_at_collision = mc.current_time() +
      CGAL::sqrt(CGAL::squared_distance(mc.current_position()->point(), collision_point)) / mc.speed();
    FT foreign_time_at_collision = time_at_fmc_track_source +
      CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(), collision_point)) / fmc.speed();

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << "Location never seen before, corresponds to point ("
              << collision_point << ") at time: " << time_at_collision << std::endl;
    std::cout << "time bounds: " << mc.current_time() << " || "
                                 << mc.time_at_closest_target() << std::endl;
    std::cout << "foreign time bounds: " << time_at_fmc_track_source << " || "
                                         << time_at_fmc_track_target << std::endl;
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
    // By construction, the time and foreign_time should be greater
    // than the times at the sources of the tracks (and oppositely for the targets).
    // Some numerical errors can sneak it and, if so, correct the times.
    //
    // It should only be a numerical error, that is a very small error
    if(time_at_collision < mc.current_time())
    {
      CGAL_precondition(time_at_collision + tolerance_ >= mc.current_time());
      time_at_collision = mc.current_time();
    }

    if(time_at_collision > mc.time_at_closest_target())
    {
      CGAL_precondition(time_at_collision - tolerance_ <= mc.time_at_closest_target());
      time_at_collision = mc.time_at_closest_target();
    }

    if(foreign_time_at_collision < time_at_fmc_track_source)
    {
      CGAL_precondition(foreign_time_at_collision + tolerance_ >= time_at_fmc_track_source);
      foreign_time_at_collision = time_at_fmc_track_source;
    }

    if(foreign_time_at_collision > time_at_fmc_track_target)
    {
      CGAL_precondition(foreign_time_at_collision - tolerance_ <= time_at_fmc_track_target);
      foreign_time_at_collision = time_at_fmc_track_target;
    }
#endif

    CGAL_postcondition(time_at_collision >= mc.current_time());
    CGAL_postcondition(time_at_collision <= mc.time_at_closest_target());
    CGAL_postcondition(foreign_time_at_collision >= time_at_fmc_track_source);
    CGAL_postcondition(foreign_time_at_collision <= time_at_fmc_track_target);

    if(tc.compare_collision_time_to_closest(time_at_collision) != Collision_information::LATER_THAN_CURRENT_CLOSEST_TIME)
    {
#ifdef CGAL_MOTORCYCLE_GRAPH_SNAPPING_CODE
      // value used later if we snap times/points
      const FT time_at_closest_collision_memory = tc.time_at_closest_collision;
#endif

      // Step 3:
      // Although there does not exist a point at the location of the collision,
      // this point might be at the same time from the source of the track
      // as another point due to numerical errors.
      std::pair<TPC_iterator, bool> res = mc.has_target_at_time(time_at_collision);
      if(res.second)
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
        std::cout << "Motorcycle #" << mc.id() << " already has a target at time: " << time_at_collision << std::endl;
#endif

        TPC_iterator target_point = res.first;
        CGAL_assertion(target_point->second == time_at_collision);
        Node_ptr alternate_collision = target_point->first;

        // If the times are equal, the points should be very close
        CGAL_assertion(CGAL::squared_distance(alternate_collision->point(), collision_point) < tolerance_);

        // Temporal snap: the collision is now that existing point instead
        return tc.treat_potential_collision(alternate_collision, time_at_collision, fmc.id(), foreign_time_at_collision);
      }

      std::pair<TPC_iterator, bool> fmc_res = fmc.has_target_at_time(foreign_time_at_collision);
      if(fmc_res.second) // there is already a target at that time
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
        std::cout << "Foreign motorcycle #" << fmc.id() << " already has a target at time: " << foreign_time_at_collision << std::endl;
#endif

        TPC_iterator target_point = fmc_res.first;
        Node_ptr alternate_foreign_collision = target_point->first;
        CGAL_assertion(alternate_foreign_collision->face() == fmc.current_face());
        CGAL_assertion(target_point->second == foreign_time_at_collision);

        // Temporal snap: the collision is now that existing point instead
        return tc.treat_potential_collision(alternate_foreign_collision, time_at_collision, fmc.id(), foreign_time_at_collision);
      }

      // At this point, we have a new location at an unknown time...
#ifdef CGAL_MOTORCYCLE_GRAPH_SNAPPING_CODE
      // Step 4:
      // ... but maybe there exists another point that is very close! Check for it,
      // and if needed, snap the new location (and the time) to it.
      std::pair<Node_ptr, bool> is_snappable = find_close_existing_point(collision_location, collision_point);
      if(is_snappable.second) // successful snapping
      {
        FT visiting_time = time_at_collision;
        // the call to this function will modify 'visiting_time' if the point of snapping is already is visited by 'mc'
        is_snappable.first->has_motorcycle(mc.id(), time_at_collision - tolerance_, time_at_collision + tolerance_, visiting_time);

        FT foreign_visiting_time = foreign_time_at_collision;
        // the call to this function will modify 'foreign_visiting_time' if the point of snapping is already is visited by 'fmc'
        is_snappable.first->has_motorcycle(fmc.id(), foreign_time_at_collision - tolerance_,
                                           foreign_time_at_collision + tolerance_, foreign_visiting_time);

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
        std::cout << "successful snapping to " << std::endl << *(is_snappable.first) << std::endl;
        std::cout << "new times: " << visiting_time << " " << time_at_closest_collision_memory << std::endl;
#endif

        // We have snapped so we are igoring times that we had formely set up as best, but
        // we still need to ensure that those times are better than the previous one.
        CGAL_assertion(visiting_time <= time_at_closest_collision_memory);

        return tc.treat_potential_collision(is_snappable.first, visiting_time, fmc.id(), foreign_visiting_time);
      }
#endif

      // Couldn't snap to anything, 'collision_location' is definitely a new point
      return tc.treat_potential_collision(collision_location, time_at_collision, fmc.id(), foreign_time_at_collision);
    }
  }

  return Collision_information::NO_COLLISION;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_complete_track(const Motorcycle& mc,
                                   const Track_segment& fmc_track,
                                   Collision_information& tc)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "-----------------------------" << std::endl;
  std::cout << "¤ Checking for intersection with the complete track of motorcycle #" << fmc_track.motorcycle_id() << std::endl;
#endif

  // 'false' because the motorcycle is not moving on that track
  return find_collision_between_tracks(mc, fmc_track, false /*not moving*/, tc);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_live_motorcycle(Motorcycle& mc,
                                    const Motorcycle& fmc,
                                    Collision_information& tc)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "-----------------------------" << std::endl;
  std::cout << "¤ Checking for intersection with live motorcycle #" << fmc.id() << std::endl;
#endif

  // Uninitialized motorcycles have a degenerate track living in the future)
  if(!fmc.is_initialized())
     return Collision_information::NO_COLLISION;

  if(mc.id() == fmc.id() || // the motorcycles must be different
     mc.current_face() != fmc.current_face() || // the motorcycles must be in the same face
     fmc.status() == Motorcycle::CRASHED) // the foreign motorcycle must be in motion
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
    std::cout << " ignoring fmc..." << std::endl;
    std::cout << "  > motorcycles #" << mc.id() << " and #" << fmc.id() << std::endl;
    std::cout << "  > faces: " << mc.current_face() << " and " << fmc.current_face() << std::endl;
    std::cout << "  > crashed status: " << fmc.status() << std::endl;
#endif
    return Collision_information::NO_COLLISION;
  }

  CGAL_assertion(!fmc.targets().empty());

  Track_segment fmc_track(fmc.id(),
                          fmc.current_position(), fmc.current_time(),
                          fmc.closest_target(), fmc.time_at_closest_target());

  // 'true' because fmc is currently moving on that track
  return find_collision_between_tracks(mc, fmc_track, true /*moving on track*/, tc);
}

// search for a possible collision with another motorcycle between the current
// position of mc and the next target
template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision(Motorcycle& mc, Collision_information& tc)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "~~~~~~~~~X ?" << std::endl;
  std::cout << "Checking for collisions on motorcycle #" << mc.id() << "'s track" << std::endl
            << "Currently on face: " << mc.current_face() << std::endl;
#endif

  // Some sanity checks
  CGAL_precondition(mc.status() != Motorcycle::CRASHED);
  CGAL_precondition(!mc.targets().empty());
  CGAL_precondition(mc.current_face() == mc.closest_target()->face());

  if(mc.is_tentative_track_degenerate())
    return Collision_information::NO_COLLISION;

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "MC tentative track:" << std::endl
            << "source: " << &*(mc.current_position()) << " " << *(mc.current_position()) << std::endl
            << "target: " << &*(mc.closest_target()) << " " << *(mc.closest_target()) << std::endl;
#endif

  // Checking for intersection is done in the following steps:
  // - 1: Check with complete tracks in the face
  // - 2: Check the motorcycles that are currently moving in the face ("live")
  // - 3: Check for intersections at the border of the face with tracks in adjacent faces

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_ COMPLETE TRACKS _-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
#endif

  Collision_return res = Collision_information::NO_COLLISION;

  // Step 1: check complete tracks
  const face_descriptor mc_fd = mc.current_face();
  TFM_iterator it = track_face_map_.find(mc_fd);
  if(it != track_face_map_.end())
  {
    const Track_segment_ptr_container& face_tracks = it->second;

    typename Track_segment_ptr_container::const_iterator ts_it = face_tracks.begin();
    typename Track_segment_ptr_container::const_iterator ts_end = face_tracks.end();
    for(; ts_it!=ts_end; ++ts_it)
    {
      const Track_segment& ts = *(*ts_it);
      Collision_return r = find_collision_with_complete_track(mc, ts, tc);

      if(r == Collision_information::COLLISION) // keep track of the segment as we might have to split it
        tc.foreign_collisions.front().foreign_track_ptr = *ts_it;

      res = res | r;
    }
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_ LIVE MOTOS -_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
#endif

  // Step 2: check incomplete tracks (path of a motorcycle currently moving in the same face)
  MCC_it fmc_it = motorcycles().begin(), fmc_end = motorcycles().end();
  for(; fmc_it!=fmc_end; ++fmc_it)
  {
    Motorcycle& fmc = motorcycle(fmc_it);
    res = res | find_collision_with_live_motorcycle(mc, fmc, tc);
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_COLLISION_VERBOSE
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_- FOREIGNERS _-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
#endif

  // @todo check if an intersection has already been found and AND the track is
  // not on a border:
  // - the source is not on a halfedge
  // - the source and the destination are on different halfedges
  // And, if so, skip below ?
  // with the reasoning that an intersection within the face is necessarily
  // earlier or at the same time as an intersection on the border.

  // Step 3: check adjacent faces
  res = res | find_collision_with_foreign_motorcycles(mc, tc);

  return res;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
generate_enclosing_face()
{
  // generate a bbox that includes all known positions and all crash points
  // 2D only for now @todo
  CGAL_precondition(Geom_traits::dimension() == 2);
  CGAL_assertion(false);

#if 0
  Bbox bbox;

  MCC_it mc_it = motorcycles().begin(), fmc_end = motorcycles().end();
  for(; fmc_it!=fmc_end; ++fmc_it)
  {
    Motorcycle& mc = ++mc_it;
    bbox += mc.input_source().bbox();

    if(mc.input_destination() != boost::none)
      bbox += mc.input_destination()->bbox();

    // this part can be done in O(nlogn) by sorting according to the slopes
    // (farthest intersections happen when the slopes of the motorcycles are close)
    // Use K::Compare_slope_2
    for(motorcycles)
    {
      // segment - segment, segment - ray, or ray-ray intersections @todo
    }
  }

  // Slightly increase the size of the bbox to strictly contain all the (collision) points

  // Manually create the mesh with Euler operations make_rectange()

#endif
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
has_motorcycle_reached_crashing_point(const Motorcycle& mc) const
{
  // A tedious detail:
  // Motorcycles starting from mc's current position at the same time as 'mc' should not block 'mc'...
  // Since such a new motorcycle 'mc2' might have been initialized before 'mc' reached
  // the current position, we might have already blocked 'mc.current_position()' ('== mc2.origin()').
  // If this is the case, we have to ignore the 'blocked' flag...

  bool is_blocked_point = mc.has_reached_blocked_point();

  const FT earliest_time = mc.current_position()->earliest_visiting_time();
  CGAL_assertion(earliest_time <= mc.current_time());

  if(is_blocked_point && earliest_time == mc.current_time())
  {
    typedef std::vector<int>                           Motorcycle_IDS_container;

    std::vector<int> earliest_motorcycles;
    mc.current_position()->earliest_motorcycles(std::back_inserter(earliest_motorcycles));

    Motorcycle_IDS_container::const_iterator id_it = earliest_motorcycles.begin(),
                                             id_end = earliest_motorcycles.end();
    for(; id_it != id_end; ++id_it)
    {
      const int mc2_id = *id_it;
      const Motorcycle& mc2 = motorcycle(mc2_id);

      // The trick here is that if the times are equal a motorcycle that is reaching
      // that point and not departing from it will crash (has_simultaneous_collision being 'true').
      // Thus, if there is a motorcycle that is at this point and not crashed, it's leaving.
      if(earliest_time == mc2.current_time() &&
         mc2.current_position() == mc.current_position() &&
         mc2.status() != Motorcycle::CRASHED)
      {
        is_blocked_point = false;
      }
      else
      {
        // Legitimate point blockade
        is_blocked_point = true;
        break;
      }
    }
  }

  // The function below ignores new starting motorcycles.
  bool at_simultaneous_collision_point = has_simultaneous_collision(mc.current_position());

  return (at_simultaneous_collision_point || is_blocked_point);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
initialize_motorcycle(Motorcycle& mc)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "      _" << std::endl;
  std::cout << "    D/_     _" << std::endl;
  std::cout << "    /(__`=-/" << std::endl;
  std::cout << "  (o)     (o)" << std::endl;
  std::cout << "Initializing motorcycle #" << mc.id() << std::endl;
#endif

  // Add the origin to the node dictionary
  const Point_or_location& input_origin = mc.input_origin();
  add_origin_node(mc, input_origin);

  // Compute the destination if needed, and add it to the node dictionary
  Optional_point_or_location& input_destination = mc.input_destination();
  add_destination_node(mc, input_destination);

  CGAL_assertion(mc.origin() != Node_ptr());
  CGAL_assertion(mc.destination() != Node_ptr());

  mc.add_target(mc.origin(), mc.time_at_origin());
  if(mc.origin() != mc.destination())
    mc.add_target(mc.destination(), mc.time_at_destination());

  // Some sanity checks
  CGAL_postcondition(mc.origin() != Node_ptr());
  CGAL_postcondition(mc.destination() != Node_ptr());
  CGAL_postcondition(mc.time_at_destination() >= mc.time_at_origin());
  CGAL_postcondition(PMP::is_in_face(mc.origin()->location(), mesh()));
  CGAL_postcondition(PMP::is_in_face(mc.destination()->location(), mesh()));
  CGAL_postcondition(mc.origin()->face() == mc.destination()->face());

  CGAL_postcondition(mc.has_target(mc.origin()).second);
  CGAL_postcondition(mc.has_target(mc.destination()).second);
  CGAL_postcondition(mc.origin()->has_motorcycle(mc.id(), mc.time_at_origin()));
  CGAL_postcondition(mc.destination()->has_motorcycle(mc.id(), mc.time_at_destination()));

#ifdef CGAL_MOTORCYCLE_GRAPH_OUTPUT
  mc.output_origin_and_destination();
#endif
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
initialize_tracing()
{
  // if no mesh has been given in input, generate a mesh made of a single quad face
  // that contains all the interesting motorcycle interactions (i.e. crashes)
  if(!is_mesh_provided)
  {
    // @todo reset all motorcycles and existing trace
    generate_enclosing_face();
  }

  motorcycle_pq_.initialize(motorcycles_);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Face_location
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
locate(const Point& p) const
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  CGAL_precondition(!aabb_tree().empty());

  // An AABB tree is a 3D structure, so we need to convert the point to a Point_3.
  // If the point is already a Point_3, this doesn't do anything.
  Point_to_Point_3 to_p3;
  const typename Point_to_Point_3::Point_3& p3 = to_p3(p);

  CGAL_static_assertion((boost::is_same<typename Point_to_Point_3::Point_3,
                                        typename AABB_tree::AABB_traits::Point_3>::value));

  Face_location loc = PMP::locate_with_AABB_tree(p3, aabb_tree(), mesh(),
                                                 parameters::vertex_point_map(aabb_tree_vpm_));

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  PMP::internal::snap_location_to_border<Triangle_mesh>(loc, tolerance_);
#endif

  return loc;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
trace_graph()
{
  initialize_tracing();

  while(!motorcycle_pq_.empty())
  {
    // Get the earliest available event
    Motorcycle_PQE pqe = motorcycle_pq_.top();
    Motorcycle& mc = pqe.motorcycle();
    const FT next_time = mc.time_at_closest_target();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Driving priority queue:" << std::endl << motorcycle_pq_ << std::endl;
#endif
    std::cout << "Driving priority queue size: " << motorcycle_pq_.size();
    std::cout << " closest time: " << next_time << ")" << std::endl << std::endl;

    if(next_time > latest_event_time_)
      latest_event_time_ = next_time;

    CGAL_precondition(mc.status() == Motorcycle::IN_MOTION);

    // Motorcycles are only initialized as late as possible to not wrongly block other traces
    if(!mc.is_initialized())
      initialize_motorcycle(mc);

    // Move the motorcycle to the closest target, which becomes its confirmed position
    drive_to_closest_target(mc);

    // Process where we are and what will be done next
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "At point:" << std::endl
              << " - blocked: " << mc.has_reached_blocked_point() << std::endl
              << " - simultaneous collision: " << has_simultaneous_collision(mc.current_position()) << std::endl;
#endif

    // Note: motorcycles cannot be blocked at their origin (if the origin is a former destination,
    // the intersection would have been detected when the origin was a destination).
    //
    // Exception: collinear directions, but this is checked when we search
    // for collisions on the tentative track
    if(mc.has_left_starting_position() &&
       has_motorcycle_reached_crashing_point(mc))
    {
      crash_motorcycle(mc);
    }
    else if(mc.current_position() == mc.destination())
    {
      if(mc.is_destination_final())
      {
        crash_motorcycle(mc);
      }
      else // Reached the destination, but not done driving yet
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "Reached destination: " << mc.destination()->point();
        std::cout << " Now computing motorcycle's next path..." << std::endl;
#endif
        // Clear any unnecessary targets that might have been built
        mc.clear_targets();

        // New destination to be filled below
      }
    }

    if(mc.status() != Motorcycle::CRASHED && mc.targets().empty())
    {
      if(!initialize_next_path(mc))
      {
        // Couldn't find a new destination
        crash_motorcycle(mc, true /* kill_regardless_of_stop_predicate */);
      }
    }

    if(mc.status() != Motorcycle::CRASHED) // The motorcycle continues!
    {
      // A bunch of output parameters are regrouped into the 'Collision_information' struct,
      // which describes the best (i.e. closest to mc.current_position()) collision.
      Collision_information tc(*this, mc.time_at_closest_target() /*maximum allowed time*/);

      // Check for potential collisions in the tentative track of 'mc'
      Collision_return res = find_collision(mc, tc);

      if(res == Collision_information::COLLISION)
      {
        CGAL_assertion(tc.found_collision());
        treat_collision(mc, tc);
      }
      else
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        CGAL_assertion(res == Collision_information::NO_COLLISION);
        std::cout << " No collision was found!" << std::endl;
#endif
      }

      motorcycle_pq_.update(mc);
    }

    // Block the point that we have just reached
    mc.current_position()->block();
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Finished tracing" << std::endl;
#endif
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
template<typename VertexNodeMap, typename EdgeTrackMap>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
assemble_graph(VertexNodeMap& vnmap,
               EdgeTrackMap& etmap,
               const bool construct_faces)
{
  // below calls 'operator()' and constructs the graph
  internal::Motorcycle_graph_builder<Self>(*this)(vnmap, etmap, construct_faces);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
treat_collision(Motorcycle& mc, const Collision_information& collision_info)
{
  CGAL_precondition(collision_info.closest_collision != boost::none);
  const Node_ptr_or_Face_location& closest_collision = *(collision_info.closest_collision);

  // Insert the collision point in the dictionary, if needed.
  Node_ptr collision;
  if(const Face_location* loc = boost::get<Face_location>(&closest_collision))
  {
    // Motorcycle info will be added later.
    std::pair<Node_ptr, bool> entry = nodes().insert(*loc, mesh());
    collision = entry.first;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    if(!entry.second)
    {
      std::cerr << "Warning: collision location actually already existed in the dictionary:"
                << std::endl << *(entry.first) << std::endl;
    }
#endif
  }
  else
  {
    collision = boost::get<Node_ptr>(closest_collision);
  }

  const face_descriptor fd = mc.current_face();
  Node_ptr collision_in_fd = nodes().get_sibling(collision, fd);
  const FT time_at_collision = collision_info.time_at_closest_collision;

  bool is_collision_at_current_position = (collision_in_fd == mc.current_position());
  bool is_collision_strictly_within_tentative_track = (!is_collision_at_current_position &&
                                                       collision_in_fd != mc.closest_target());

  if(collision_info.must_crash)
  {
    crash_motorcycle(mc);
  }
  else if(is_collision_strictly_within_tentative_track)
  {
    treat_collision(mc, collision_in_fd, time_at_collision);
  }
  else if(is_collision_at_current_position)
  {
    // @fixme a bit hackish...
    // We found a collision at the current point. If it were known collision,
    // it would have been ignored. Since we are here, it means that it is new information.
    // We can't however drive forward just yet because we have to look for collisions
    // in the tentative track beyond the current position.
    // The collision at the current position is now known, and it will be ignored.
    mc.add_target(mc.current_position(), mc.current_time());
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << std::endl << "[[ Post treat_collision() ... ]]" << std::endl;
  std::cout << "Motorcycle involved:" << std::endl << mc << std::endl;
  std::cout << "Collision point:" << std::endl << *collision_in_fd << std::endl;
  std::cout << collision_info.foreign_collisions.size() << " foreign motorcycles" << std::endl;
#endif

  // Now add the collision node to the foreign motorcycles too
  FCC_cit fccit = collision_info.foreign_collisions.begin(),
          fccend = collision_info.foreign_collisions.end();
  for(; fccit!=fccend; ++fccit)
  {
    const Foreign_collision_information& foreign_collision = *fccit;
    const int foreign_motorcycle_id = foreign_collision.fmc_id;

    Motorcycle& fmc = motorcycle(foreign_motorcycle_id);
    const FT foreign_time_at_collision = foreign_collision.foreign_time_at_closest_collision;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Foreign_motorcycle:" << std::endl << fmc << std::endl
              << "Foreign time at collision: " << foreign_time_at_collision << std::endl;
#endif
    const Track_segment_ptr foreign_track = foreign_collision.foreign_track_ptr;
    const bool is_foreign_motorcycle_moving = (foreign_track == Track_segment_ptr());

    if(foreign_collision.must_crash)
    {
      CGAL_assertion(is_foreign_motorcycle_moving); // must be a live foreign motorcycle
      CGAL_assertion(time_at_collision == mc.time_at_origin());
      CGAL_assertion(foreign_time_at_collision == time_at_collision);

      crash_motorcycle(fmc);
      continue;
    }

    const face_descriptor ffd = is_foreign_motorcycle_moving ? fmc.current_face()
                                                             : foreign_track->source()->face();
    CGAL_assertion(ffd != boost::graph_traits<Triangle_mesh>::null_face());

    Node_ptr collision_in_ffd = nodes().get_sibling(collision, ffd);

    // Some sanity tests
    CGAL_postcondition(collision_in_fd->face() == fd);
    CGAL_postcondition(collision_in_ffd->face() == ffd);
    CGAL_postcondition_code(if(fd != ffd) {)
    CGAL_postcondition(collision_in_ffd->is_sibling(collision_in_fd->location()));
    CGAL_postcondition(collision_in_fd->is_sibling(collision_in_ffd->location())); // just to be extra sure
    CGAL_postcondition_code(})

    if(!is_foreign_motorcycle_moving)
    {
      CGAL_assertion(foreign_track != Track_segment_ptr());
      CGAL_assertion(foreign_track->motorcycle_id() == foreign_motorcycle_id);
      CGAL_assertion(foreign_track->time_at_source() <= foreign_time_at_collision);
      CGAL_assertion(foreign_track->time_at_target() >= foreign_time_at_collision);
    }
    // end of sanity tests

    treat_foreign_collision(fmc, collision_in_ffd, foreign_time_at_collision, foreign_track);

    // The target list of the foreign motorcycle has been modified and the queue must be updated
    if(fmc.status() == Motorcycle::IN_MOTION)
      motorcycle_pq_.update(fmc);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << std::endl << "[[ Post [foreign] treat collision... ]]" << std::endl;
    std::cout << "Foreign motorcycle involved:" << std::endl << fmc << std::endl;

    std::cout << "collision point:" << std::endl << *collision << std::endl;
#endif
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
treat_collision(Motorcycle& mc, Node_ptr collision, const FT time_at_collision)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << std::endl << "+++++++++++ Treat collision +++++++++++" << std::endl;
  std::cout << " - motorcycle:" << std::endl << mc << std::endl
            << " - collision point: " << &*collision << std::endl << *(collision) << std::endl
            << " - time at collision: " << time_at_collision << std::endl;
#endif

  CGAL_precondition(collision != Node_ptr());
  CGAL_precondition(collision->face() == mc.current_face());

  CGAL_precondition(time_at_collision > mc.current_time());
  CGAL_precondition(time_at_collision < mc.time_at_closest_target());

  CGAL_assertion(!collision->has_motorcycle(mc.id(), time_at_collision));

  collision->add_motorcycle(mc.id(), time_at_collision);
  mc.add_target(collision, time_at_collision);
  CGAL_postcondition(mc.has_target_at_time(collision, time_at_collision));

  // Call the halving structure to create a new point
  std::pair<Node_ptr, FT> halving_entity =
      compute_halving_point(mc, mc.current_position(), mc.current_time(),
                            collision, time_at_collision);
  Node_ptr halving_point = halving_entity.first;
  const FT time_at_halving_point = halving_entity.second;

  // Degeneracies should have been caught before
  CGAL_postcondition(halving_point != mc.current_position() && halving_point != collision);

  halving_point->add_motorcycle(mc.id(), time_at_halving_point);
  mc.add_target(halving_point, time_at_halving_point);
  CGAL_postcondition(mc.has_target_at_time(halving_point, time_at_halving_point));
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
treat_foreign_collision(Motorcycle& fmc, Node_ptr foreign_collision,
                        const FT foreign_time_at_collision, Track_segment_ptr foreign_track)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "+++++++++++ Treat FOREIGN collision +++++++++++" << std::endl;
  std::cout << " - foreign motorcycle:" << std::endl << fmc << std::endl
            << " - foreign collision point: " << &*foreign_collision << std::endl << *(foreign_collision) << std::endl
            << " - foreign time at collision: " << foreign_time_at_collision << std::endl;
#endif

  // Some sanity checks
  CGAL_precondition(foreign_collision != Node_ptr());

  // If 'fmc' is not moving, ignore collisions after its final destination
  if(fmc.status() != Motorcycle::IN_MOTION && foreign_time_at_collision > fmc.current_time())
    return;

  // Nothing to do if it is already a known collision point
  if(foreign_collision->has_motorcycle(fmc.id(), foreign_time_at_collision))
  {
    // If it's in the tentative track or beyond, it should be a target of 'fmc'
    CGAL_assertion_code(if(fmc.status() == Motorcycle::IN_MOTION && foreign_time_at_collision > fmc.current_time()))
    CGAL_assertion((fmc.has_target_at_time(foreign_time_at_collision)).second);

    return;
  }

  // It is useful to know that the collision point is on the foreign track,
  // even if the collision point is on the confirmed part of the track.
  foreign_collision->add_motorcycle(fmc.id(), foreign_time_at_collision);

  // If the collision point is not on the confirmed track, add it as target
  if(foreign_time_at_collision > fmc.current_time())
  {
    CGAL_assertion(fmc.status() == Motorcycle::IN_MOTION); // we filtered the other cases above

    fmc.add_target(foreign_collision, foreign_time_at_collision);
    CGAL_postcondition(fmc.has_target_at_time(foreign_collision, foreign_time_at_collision));

    // Call the halving structure to create a new point
    std::pair<Node_ptr, FT> foreign_halving_entity =
      compute_halving_point(fmc, fmc.current_position(), fmc.current_time(),
                            foreign_collision, foreign_time_at_collision);
    Node_ptr foreign_halving_point = foreign_halving_entity.first;
    const FT foreign_time_at_halving_point = foreign_halving_entity.second;

    // Degeneracies should have been caught before
    CGAL_postcondition(foreign_halving_point != fmc.current_position() &&
                       foreign_halving_point != foreign_collision);

    if(!foreign_halving_point->has_motorcycle(fmc.id()))
    {
      foreign_halving_point->add_motorcycle(fmc.id(), foreign_time_at_halving_point);
      fmc.add_target(foreign_halving_point, foreign_time_at_halving_point);
      CGAL_postcondition(fmc.has_target_at_time(foreign_halving_point, foreign_time_at_halving_point));
    }
  }
  // else, foreign time <= current_time --> collision point is on the confirmed foreign track
  else
  {
    // It is on the confirmed track, so it is in the past --> we can already block it
    foreign_collision->block();

    // Add it to the track of the foreign motorcycle
    CGAL_assertion(foreign_track != Track_segment_ptr());
    Track_segment_ptr ts_ptr = fmc.track().split_track_segment(foreign_track,
                                                               foreign_collision,
                                                               foreign_time_at_collision);

    // add the new track segment to the face map
    add_track_segment_to_map(foreign_track->face(), ts_ptr);
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
template<typename VertexNodeMap, typename EdgeTrackMap>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
construct_motorcycle_graph(VertexNodeMap& vnmap,
                           EdgeTrackMap& etmap,
                           const bool construct_faces)
{
  trace_graph();

#ifdef CGAL_MOTORCYCLE_GRAPH_OUTPUT
  output_tracks();
  output_all_points();
#endif

  assemble_graph(vnmap, etmap, construct_faces);

#ifdef CGAL_MOTORCYCLE_GRAPH_OUTPUT
  if(construct_faces)
  {
    std::ofstream out("motorcycle_graph.off");
    out.precision(20);
    CGAL::write_off(out, *graph_);
    out.close();
  }
  else
  {
    std::ofstream out("motorcycle_graph.polylines.txt");
    print_motorcycle_graph(out);
  }
#endif

  // @tmp
//  CGAL_postcondition(internal::is_valid(*this));
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
std::pair<typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Node_ptr, bool>
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_close_existing_point(const Face_location& location, const Point& p,
                          const bool allow_same_location) const
{
  const face_descriptor fd = location.first;

  // @todo Brute force for now, need a kd tree
  Node_ptr dit = nodes().begin(), end = nodes().end();
  for(; dit!=end; ++dit)
  {
    // Filter the points from different faces (not a problem for points on face borders
    // since they have a representative in each incident face)
    if(dit->face() != fd)
      continue;

    // - Ignore the same location if it is found
    if(!allow_same_location && dit->location() == location)
      continue;

    if(CGAL::squared_distance(dit->point(), p) <= tolerance_)
    {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "!!! new point: (" << p << ") is close enough to existing point: " << std::endl << *dit << std::endl;
#endif
      return std::make_pair(dit, true);
    }
  }

  return std::make_pair(Node_ptr(), false);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
snap_to_existing_point(Node_ptr& node, const Node_ptr& t)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Snapping: " << std::endl << *node << std::endl
            << " to existing point: " << std::endl << *t << std::endl;
#endif

  if(!node->has_motorcycles())
  {
    CGAL_expensive_assertion_code(MCC_cit mc_it = motorcycles().begin();)
        CGAL_expensive_assertion_code(MCC_cit mc_end = motorcycles().end();)
        CGAL_expensive_assertion_code(for(; mc_it!=mc_end; ++mc_it))
        CGAL_expensive_assertion(!(motorcycle(mc_it).has_target(node)).second);
    nodes().erase(node);
  }
  else
  {
    // Don't want to snap something that is already being used at another point.
    // If we are here, it means there is an issue in the creation of one
    // of the two nodes close to one another.
    CGAL_assertion(false);
  }

  node = t;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
output_all_points() const
{
  typename Nodes::Node_container::const_iterator dit = nodes().begin();
  typename Nodes::Node_container::const_iterator end = nodes().end();

  std::stringstream oss;
  oss << "results_" << geom_traits().dimension() << "/dictionary_points.xyz" << std::ends;
  std::ofstream os(oss.str().c_str());
  if(!os.good())
    return;

  os.precision(20);

  for(; dit!=end; ++dit)
  {
    os << dit->point();
    if(geom_traits().dimension() == 2) // The '.xyz' format expects 3D points
      os << " 0";
    os << '\n';
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
output_tracks() const
{
  MCC_cit mc_it = motorcycles().begin(), mc_end = motorcycles().end();
  for(; mc_it!=mc_end; ++mc_it)
  {
    std::stringstream oss;
    oss << "motorcycle_track_" << mc_it->id() << ".polylines.txt" << std::ends;
    std::ofstream out(oss.str().c_str());
    out.precision(17);

    typename Track::const_iterator tscit = mc_it->track().begin();
    typename Track::const_iterator tsend = mc_it->track().end();
    for(; tscit!=tsend; ++tscit)
      out << "2 " << tscit->source()->point() << " " << tscit->target()->point() << '\n';

    out.close();
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
print_motorcycle_graph(std::ofstream& out) const
{
  out.precision(17);

  typedef typename property_map_selector<Face_graph, CGAL::vertex_point_t>::const_type VPMap;
  VPMap vpm = get_const_property_map(boost::vertex_point, graph());

  typename boost::graph_traits<Face_graph>::edge_iterator eit, eend;
  boost::tie(eit, eend) = edges(graph());

  std::cout << num_vertices(graph()) << " vertices and " << num_edges(graph()) << " edges" << std::endl;

  while(eit != eend)
  {
    out << "2 " << get(vpm, source(*eit, graph())) << " ";
    if(Geom_traits::dimension() == 2)
      out << "0 ";

    out << get(vpm, target(*eit, graph())) << " ";
    if(Geom_traits::dimension() == 2)
      out << "0 ";

    out << '\n';

    ++eit;
  }

  out.close();
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H
