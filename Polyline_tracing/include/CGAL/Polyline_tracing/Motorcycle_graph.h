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

#define CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE // @tmp

#include <CGAL/Polyline_tracing/Motorcycle.h>
#include <CGAL/Polyline_tracing/Motorcycle_graph_node_dictionary.h>
#include <CGAL/Polyline_tracing/Motorcycle_priority_queue.h>
#include <CGAL/Polyline_tracing/Track.h>
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

#include <boost/foreach.hpp> // @fixme CGAL_foreach everywhere
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

namespace internal {

// This struct regroups all useful information on a potential intersection
// @todo completely rework this
template<typename MotorcycleGraph>
struct Collision_information
{
  typedef typename MotorcycleGraph::Triangle_mesh             Triangle_mesh;
  typedef typename MotorcycleGraph::FT                        FT;
  typedef typename MotorcycleGraph::Node_ptr                  Node_ptr;
  typedef typename MotorcycleGraph::face_descriptor           face_descriptor;
  typedef typename MotorcycleGraph::Face_location             Face_location;
  typedef typename MotorcycleGraph::Barycentric_coordinates   Barycentric_coordinates;

  typedef typename MotorcycleGraph::Track_segment_ptr         Track_segment_ptr;

  // Constructor
  Collision_information(const FT max_time_at_collision)
    :
      maximum_time_at_collision(max_time_at_collision),

      // information related to the current face
      is_closest_collision_already_in_dictionary(false),
      closest_collision(),
      closest_collision_location(std::make_pair(boost::graph_traits<Triangle_mesh>::null_face(),
                                                Barycentric_coordinates())),
      time_at_closest_collision(std::numeric_limits<FT>::max()),

      // information related to the neighboring foreign face
      fmc_id(-1),
      is_foreign_motorcycle_moving_on_track(false),
      is_foreign_motorcycle_in_different_face(false),
      foreign_motorcycle_face(boost::graph_traits<Triangle_mesh>::null_face()),
      foreign_time_at_closest_collision(std::numeric_limits<FT>::max()),
      foreign_track()
  { }

  // Functions
  bool found_collision() const
  {
    // Either a node is provided, or the location should be provided
    return (is_closest_collision_already_in_dictionary ||
            closest_collision_location.first != boost::graph_traits<Triangle_mesh>::null_face());
  }

  // Check if the times provided passed in arguments correspond to a collision
  // earlier than the current best, bounded by the maximum time (time at closest target)
  bool is_collision_earlier_than_current_best(const FT time_at_collision,
                                              const FT foreign_time_at_collision,
                                              const bool is_fmc_moving) const
  {
    if(time_at_collision > maximum_time_at_collision)
      return false;

    bool is_equal_collision_time_but_given_priority = false;
    if(time_at_collision == time_at_closest_collision)
    {
      // Note that everything here correspond to the _same_ collision time for
      // two foreign motorcycles.
      //
      // Below serves to handle intersecting multiple tracks at the same point:
      // we want to make sure all the motorcycles correctly visit the common point.
      // Thus:
      // - If we have to choose between two _active_ foreign motorcycles, prefer the foreign
      // motorcycle that would reach first. This is because that motorcycle has by definition
      // a higher chance of reaching and blocking the collision point first.
      // - If one of the motorcycle is moving, give it priority regardless of foreign time.
      // This is because the active motorcycle's tentative track will then be cut
      // in half and we will simply cut tentative tracks till there's only a fixed track left.
      // - If both are inactive, it doesn't matter so keep the current one.
      if(is_fmc_moving && is_foreign_motorcycle_moving_on_track)
        is_equal_collision_time_but_given_priority = (foreign_time_at_collision < foreign_time_at_closest_collision);
      if(is_fmc_moving && !is_foreign_motorcycle_moving_on_track)
        is_equal_collision_time_but_given_priority = true;
      if(!is_fmc_moving && is_foreign_motorcycle_moving_on_track)
        is_equal_collision_time_but_given_priority = false;
      if(!is_fmc_moving && !is_foreign_motorcycle_moving_on_track)
        is_equal_collision_time_but_given_priority = false;
    }

    const bool is_collision_earlier = (time_at_collision < time_at_closest_collision);
    const bool is_better = is_collision_earlier || is_equal_collision_time_but_given_priority;

    if(is_better)
    {
      CGAL_assertion(time_at_collision <= time_at_closest_collision);
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "New earliest collision times: " << time_at_collision << " || "
                                                    << foreign_time_at_collision;
      std::cout << " [previously: " << time_at_closest_collision << " || "
                                    << foreign_time_at_closest_collision << "]" << std::endl;
#endif
    }

    return is_better;
  }

  void reset()
  {
    // information relative to the current face
    is_closest_collision_already_in_dictionary = false;
    closest_collision = Node_ptr();
    closest_collision_location = std::make_pair(boost::graph_traits<Triangle_mesh>::null_face(),
                                                Barycentric_coordinates());
    time_at_closest_collision = std::numeric_limits<FT>::max();

    // information relative to the neighboring foreign face
    fmc_id = -1;
    is_foreign_motorcycle_moving_on_track = false;
    is_foreign_motorcycle_in_different_face = false,
    foreign_motorcycle_face = boost::graph_traits<Triangle_mesh>::null_face();
    foreign_time_at_closest_collision = std::numeric_limits<FT>::max();
    foreign_track = Track_segment_ptr();
  }

public:
  const FT maximum_time_at_collision;

  bool is_closest_collision_already_in_dictionary;
  Node_ptr closest_collision;
  Face_location closest_collision_location;
  FT time_at_closest_collision;

  std::size_t fmc_id;
  bool is_foreign_motorcycle_moving_on_track;
  bool is_foreign_motorcycle_in_different_face;
  face_descriptor foreign_motorcycle_face;
  FT foreign_time_at_closest_collision;
  Track_segment_ptr foreign_track;
};

} // namespace internal

// @todo allow multiple waves of motorcycles (for features)
// @todo handle 'drive till time is T'
// @todo handle degenerate faces / tracks
// @todo snap input points? + Crash a motorcycle if new dest = current_pos ? (currently being done)
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
  typedef typename Geom_traits::Halfedge_graph                              Halfedge_graph;

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
  typedef typename Nodes::Node_ptr                                          Node_ptr;

  typedef typename Geom_traits::Barycentric_coordinates                     Barycentric_coordinates;
  typedef typename Geom_traits::Face_location                               Face_location;

  // Motorcycles
  //@todo slist instead ? what to do with id
  typedef MotorcycleType                                                    Motorcycle;
  typedef std::deque<Motorcycle>                                            Motorcycle_container;
  typedef typename Motorcycle_container::iterator                           MCC_it;
  typedef typename Motorcycle_container::const_iterator                     MCC_cit;

  typedef Motorcycle_priority_queue<Self>                                   Motorcycle_PQ;
  typedef Motorcycle_priority_queue_entry<Self>                             Motorcycle_PQE;

  // Location-related types
  typedef boost::variant<Point, Face_location>                              Point_or_location;
  typedef boost::optional<Point_or_location>                                Optional_point_or_location;

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

  typedef Motorcycle_track_segment<Geom_traits>                             Track_segment;
  typedef Motorcycle_track<Geom_traits>                                     Track;

  typedef typename Track::iterator                                          Track_segment_ptr;
  typedef boost::container::slist<Track_segment_ptr>                        Track_segment_ptr_container;

  // To check collisions with adjacent faces, we need to know for a face the given track
  // track segments in this face. Pointers because tracks live in their respective
  // motorcycle and there is no need to duplicate information.
  typedef boost::unordered_map<face_descriptor,
                               Track_segment_ptr_container>                 Track_face_map;
  typedef typename Track_face_map::iterator                                 TFM_iterator;
  typedef typename Track_face_map::const_iterator                           TFM_const_iterator;

  // Constructor mesh/graph building functions
  template<typename Halfedge_graph_ptr>
  void initialize_graph(Halfedge_graph_ptr graph) { graph_ = graph; }
  void initialize_graph(const boost::param_not_found);

  template<typename Triangle_mesh_ptr, typename NamedParameters>
  void initialize_mesh_and_graph(Triangle_mesh_ptr _mesh, const NamedParameters& np);
  template<typename NamedParameters>
  void initialize_mesh_and_graph(const boost::param_not_found, const NamedParameters& np);

  // Constructor
  template<typename NamedParameters>
  Motorcycle_graph(const NamedParameters& np);

  // Destructor
  ~Motorcycle_graph();

  // Main function to construct the graph
  template<typename VertexNodeMap, typename EdgeTrackMap>
  void construct_motorcycle_graph(VertexNodeMap& vnmap, EdgeTrackMap& etmap);

  void construct_motorcycle_graph()
  {
    boost::dummy_property_map dummy;
    return construct_motorcycle_graph(dummy, dummy);
  }

  // Function to add motorcycles, forwards the arguments to the constructor of 'Motorcycle'
  template<typename ... Args>
  std::size_t add_motorcycle(const Args& ... args);

  // @todo move that to private ?
  void trace_graph();

  // Access
  const Geom_traits& geom_traits() const { return gt_; }
  Triangle_mesh& mesh() { return *mesh_; }
  const Triangle_mesh& mesh() const { return *mesh_; }
  Halfedge_graph& graph() { return *graph_; }
  const Halfedge_graph& graph() const { return *graph_; }
  Nodes& nodes() { return nodes_; }
  const Nodes& nodes() const { return nodes_; }
  Motorcycle_container& motorcycles() { return motorcycles_; }
  const Motorcycle_container& motorcycles() const { return motorcycles_; }

  Motorcycle& motorcycle(const std::size_t id) {
    CGAL_precondition(id >= 0 && id < motorcycles_.size());
    return motorcycles_[id];
  }
  const Motorcycle& motorcycle(const std::size_t id) const {
    CGAL_precondition(id >= 0 && id < motorcycles_.size());
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

  // Validity & ouput
  bool is_valid() const;

  void output_all_points() const;
  void output_motorcycles_origins_and_destinations() const;
  void print_motorcycle_graph() const;

private:
  /// \param fd face in which the segment belongs
  /// \param id the id of the motorcycle
  /// \param s, t the source and target of the oriented segment
  ///
  /// \return iterator in the tracking map
  TFM_iterator add_track_segment_to_map(face_descriptor fd, const Track_segment_ptr ts);

  bool add_origin_node(Motorcycle& mc, const Point_or_location& input_origin,
                       const AABB_tree& tree, const AABB_tree_VPM& vpm);
  bool add_destination_node(Motorcycle& mc, const Optional_point_or_location& input_destination);

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
  void crash_motorcycle(Motorcycle& mc); // @todo add visitors/callbacks ?
  void crash_motorcycles_with_same_origins_and_directions();
  void drive_to_closest_target(Motorcycle& mc);

  // 0, 1, and 3 because we will use this enum in bitwise operations (see operator below)
  enum Collision_return
  {
    NO_COLLISION = 0,
    COLLISION = 1,
    SNAPPED_COLLISION_TO_EXISTING_POINT = 3
  };

  friend inline Collision_return operator|(Collision_return cr1, Collision_return cr2) {
    return static_cast<Collision_return>(static_cast<int>(cr1) | static_cast<int>(cr2));
  }

  // Snapping functions
  std::pair<Node_ptr, bool> find_close_existing_point(const Face_location& location,
                                                      const Point& p,
                                                      const bool allow_same_location = true) const;
  bool try_to_snap_to_close_existing_point(Node_ptr& e);

  // Collisions between two motorcycles in different faces
  Collision_return find_collision_with_foreign_motorcycles(Motorcycle& mc, Collision_information& tc);

  // Below, only the target of the tentative track is on a border
  // ---------------------------------------------------------------------------------
  // collect the different faces in which we seek a collision depending on the location 'dv'
  Collision_return find_collision_with_tentative_track_target_on_border(const Motorcycle& mc,
                                                                        const descriptor_variant dv,
                                                                        Collision_information& tc) const;

  // collect the motorcycles and tracks that we need to seek collisions with in the face 'ffd'
  Collision_return find_collision_with_tentative_track_target_on_border(const Motorcycle& mc,
                                                                        const descriptor_variant dv,
                                                                        const face_descriptor ffd,
                                                                        Collision_information& tc) const;
  // triage based on the validity of 'fmc' and build the foreign track
  Collision_return find_collision_with_tentative_track_target_on_border_with_live_motorcycle_on_foreign_face(const Motorcycle& mc,
                                                                                                             const descriptor_variant dv,
                                                                                                             const face_descriptor ffd,
                                                                                                             const Motorcycle& fmc,
                                                                                                             Collision_information& tc) const;
  // try to find a collision between the tentative tracks's target and the foreign track
  Collision_return find_collision_with_tentative_track_target_on_border_with_track_on_foreign_face(const Motorcycle& mc,
                                                                                                   const descriptor_variant ct_dv,
                                                                                                   const Track_segment& fmc_track,
                                                                                                   const bool is_fmc_moving_on_track,
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
                                                               const bool is_fmc_moving_on_track,
                                                               Collision_information& tc) const;
  // ---------------------------------------------------------------------------------

  // Below, find collisions in a common face
  // ---------------------------------------------------------------------------------
  // collisions between two motorcycles in the same face
  Collision_return find_collision_at_tentative_track_destination(const Motorcycle& mc,
                                                                 const Motorcycle& fmc,
                                                                 const FT fmc_visiting_time,
                                                                 const bool is_fmc_on_foreign_track,
                                                                 const bool is_fmc_moving_on_track,
                                                                 Collision_information& tc) const;
  Collision_return find_collision_between_collinear_tracks(const Motorcycle& mc,
                                                           const Segment_2& mcs,
                                                           const Motorcycle& fmc,
                                                           const Track_segment& fmc_track,
                                                           const Segment_2& fmcs,
                                                           const bool is_fmc_moving_on_track,
                                                           Collision_information& tc) const;
  Collision_return find_collision_between_tracks(const Motorcycle& mc,
                                                 const Segment_2& mcs,
                                                 const Motorcycle& fmc,
                                                 const Track_segment& fmc_track,
                                                 const bool is_fmc_moving_on_track,
                                                 Collision_information& tc) const;
  Collision_return find_collision_with_complete_track(const Motorcycle& mc,
                                                      const Segment_2& mcs,
                                                      const Track_segment& fmc_track,
                                                      Collision_information& tc);
  Collision_return find_collision_with_live_motorcycle(Motorcycle& mc,
                                                       const Segment_2& mcs,
                                                       const Motorcycle& fmc,
                                                       Collision_information& tc);

  // \return collision point (if any), time at the collision for `mc`, id of
  //         the foreign intersecting motorcycle, time at the collision for
  //         the foreign motorcycle.
  Collision_return find_collision(Motorcycle& mc, Collision_information& tc);

  // \brief If in 2D and no mesh is passed, generate a triangle that encloses
  //        all the points and all the collisions
  void generate_enclosing_face();
  bool has_motorcycle_reached_crashing_point(const Motorcycle& mc) const;
  bool has_motorcycle_reached_final_destination(const Motorcycle& mc) const;
  void initialize_motorcycles();
  void initialize_tracing();
  bool is_AABB_tree_needed() const;
  bool is_motorcycle_position_blocked(const Motorcycle& mc) const;
  Face_location locate(const Point& p, const AABB_tree& tree, const AABB_tree_VPM vpm) const;

  void treat_collision(Motorcycle& mc, const Collision_information& collision);
  void treat_collision(Motorcycle& mc, Node_ptr collision_point, const FT time_at_collision,
                       Motorcycle& fmc, Node_ptr foreign_collision_point, const FT foreign_time_at_collision,
                       Track_segment_ptr foreign_track);

  void visit_point(Motorcycle& mc, Collision_information& tc);

private:
  Geom_traits gt_;

  Nodes nodes_; // points that will be used throughout the algorithm
  Motorcycle_container motorcycles_;
  Motorcycle_PQ motorcycle_pq_; // motorcycle priority queue

  Triangle_mesh* mesh_;
  Halfedge_graph* graph_;
  bool is_mesh_provided, is_graph_provided;

  // map to store the completed tracks of the motorcycles for each face of the mesh
  Track_face_map track_face_map_;

  const FT tolerance_ = 1e-13;
};

// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
initialize_graph(const boost::param_not_found)
{
  graph_ = new Halfedge_graph();
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
    mesh_(),
    graph_(),
    is_mesh_provided(true),
    is_graph_provided(true),
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
std::size_t
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
add_motorcycle(const Args& ... args)
{
  std::size_t new_id = number_of_motorcycles();
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
  CGAL_precondition(ts->motorcycle_id() >= 0 && ts->motorcycle_id() < number_of_motorcycles());

  Track_segment_ptr_container l;
  l.push_front(ts);

  std::pair<TFM_iterator, bool> is_insert_success = track_face_map_.insert(std::make_pair(fd, l));

  if(!is_insert_success.second)
    is_insert_success.first->second.push_front(ts);

  return is_insert_success.first;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
add_origin_node(Motorcycle& mc,
                const Point_or_location& input_origin,
                const AABB_tree& tree, const AABB_tree_VPM& vpm)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  Node_ptr origin;
  Face_location origin_location;
  Point origin_point;

  if(const Point* origin_point_ptr = boost::get<Point>(&input_origin))
  {
    origin_point = *origin_point_ptr;
    origin_location = locate(origin_point, tree, vpm);
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

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  // Try to find an existing point close to that location
  std::pair<Node_ptr, bool> is_snappable = find_close_existing_point(origin_location, origin_point);
  if(is_snappable.second)
  {
    // @todo should the time be changed? If so, how?
    origin = is_snappable.first;
  }
  else
#endif
  {
    std::pair<Node_ptr, bool> is_insert_successful = nodes().insert(origin_location, origin_point, mesh());
    CGAL_assertion(is_insert_successful.second);
    origin = is_insert_successful.first;
  }

  CGAL_postcondition(origin != Node_ptr());
  mc.origin() = origin;
  mc.current_position() = mc.origin();

  mc.origin()->add_motorcycle(mc.id(), mc.current_time());

  return true;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
add_destination_node(Motorcycle& mc,
                     const Optional_point_or_location& input_destination)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  // At the start of this function, mc.origin() must already be initialized
  CGAL_precondition(mc.origin() != Node_ptr());

  if(input_destination == boost::none) // A destination was not provided
    return compute_and_set_next_destination(mc);

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
      PMP::locate_in_common_face(destination_point, origin_location, destination_location, mesh());
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

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  // Try to find an existing point close to that location
  std::pair<Node_ptr, bool> is_snappable = find_close_existing_point(destination_location, destination_point);

  if(is_snappable.second)
  {
    // @todo should the time be changed? If so, how?
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
  mc.destination()->add_motorcycle(mc.id(), mc.time_at_destination());

  return true;
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
  if(p->face() != q->face())
  {
    std::cerr << "Error: middle point computation with different faces" << std::endl;
    CGAL_assertion(false);
  }

  const Barycentric_coordinates& p_coords = p->location().second;
  const Barycentric_coordinates& q_coords = q->location().second;

  Barycentric_coordinates middle_coords = CGAL::make_array(0.5*(p_coords[0] + q_coords[0]),
                                                           0.5*(p_coords[1] + q_coords[1]),
                                                           0.5*(p_coords[2] + q_coords[2]));
  Face_location middle_loc = std::make_pair(p->face(), middle_coords);
  const FT time_at_r = 0.5 * (p_time + q_time);
  std::pair<Node_ptr, bool> entry = nodes().insert(middle_loc, mesh());

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  New middle point: (" << entry.first->point()
                                       << ") at time: " << time_at_r << std::endl;
  std::cout << "Location: " << p->face()
            << " bc: " << middle_coords[0] << " "
                       << middle_coords[1] << " "
                       << middle_coords[2] << std::endl;
#endif

  return std::make_pair(entry.first, time_at_r);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
compute_and_set_next_destination(Motorcycle& mc)
{
  Tracer_result res = mc.compute_next_destination(nodes(), mesh());

  if(!res.template get<0>()) // couldn't find a next path
    return false;

  const Node_ptr& next_origin = res.template get<1>();
  Node_ptr next_destination = res.template get<2>();
  const FT time_at_next_destination = res.template get<3>();
  const bool is_destination_final = res.template get<4>();

  // @todo try to snap the coordinates to a possibly-nearby border here.

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  try_to_snap_to_close_existing_point(next_destination);
#endif

  mc.destination() = next_destination;
  mc.time_at_destination() = time_at_next_destination;
  mc.set_destination_finality(is_destination_final);

  if(mc.destination() == mc.current_position())
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    if(!is_destination_final)
      std::cerr << "Warning: new destination is the current position but 'is_final' is set to 'no'!" << std::endl;
#endif

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
crash_motorcycle(Motorcycle& mc)
{
  if(mc.is_crashed())
    return;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~X" << std::endl;
  std::cout << "Crashing " << mc;
#endif

  mc.current_position()->block();
  mc.clear_targets();
  mc.crash();
  motorcycle_pq_.erase(mc);
}

// @todo crash null speed too (in another function)
template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
crash_motorcycles_with_same_origins_and_directions()
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Checking for motorcycles with same origins and directions" << std::endl;
#endif

  // @todo handle motorcycles starting from the same point & with same directions
  // but not on the same face...

  // brute force, for now
  // A smarter version is to sort motorcycles by direction (slope),
  // and check for consecutive entries @todo (and for surfaces...)
  MCC_it mc_it = motorcycles().begin(), mc_end = motorcycles().end();
  for(; mc_it!=mc_end; ++mc_it)
  {
    Motorcycle& mc = motorcycle(mc_it);

    if(mc.origin() == mc.destination() || mc.is_crashed())
      continue;

    MCC_it fmc_it = motorcycles().begin();
    for(; fmc_it!=mc_end; ++fmc_it)
    {
      Motorcycle& fmc = motorcycle(fmc_it);

      // Note: not ignoring crashed motorcycles in case of > 2 motorcycles with
      // same origin and destination

      if(fmc.id() == mc.id() ||
         fmc.origin() == fmc.destination() || // a degenerate track does not block anything
         mc.origin() != fmc.origin()) // must have identical origins
        continue;

      CGAL_assertion(mc.current_face() == fmc.current_face());

      Point_2 mc_s(mc.origin()->location().second[0], mc.origin()->location().second[1]);
      Point_2 mc_d(mc.destination()->location().second[0], mc.destination()->location().second[1]);
      Point_2 fmc_d(fmc.destination()->location().second[0], fmc.destination()->location().second[1]);
      Point_2 fmc_s(fmc.origin()->location().second[0], fmc.origin()->location().second[1]);

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
      // Add some tolerance to the definition of "collinearity"
      Vector_2 mc_v(mc_s, mc_d);
      Vector_2 fmc_v(fmc_s, fmc_d);

      FT mc_v_n = mc_v * mc_v;
      FT fmc_v_n = fmc_v * fmc_v;

      FT sp = geom_traits().compute_scalar_product_2_object()(mc_v, fmc_v);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "SProduct: " << sp << std::endl;
      std::cout << "SProduct normalized " << sp * sp / (fmc_v_n * mc_v_n ) << std::endl;
#endif

      if(CGAL::abs( 1 - sp * sp / (fmc_v_n * mc_v_n) ) < tolerance_)
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "Crashing degenerate motorcycles: " << mc.id() << " and " << fmc.id() << std::endl;
#endif
        // drive the motorcycles to create their (degenerate tracks)
        drive_to_closest_target(mc);
        crash_motorcycle(mc);

        drive_to_closest_target(fmc);
        crash_motorcycle(fmc);

        break;
      }
#else // CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
      // only aligned tracks block one another
      if(!geom_traits().collinear_2_object()(mc_s /* == fmc.origin()->point() */, mc_d, fmc_d))
        continue;

      std::cout << "Collinear tracks with the same origin" << std::endl;

      // Moving away from each other from the same point is allowed.
      if(geom_traits().angle_2_object()(mc_s, mc_d, fmc_s, fmc_d) == CGAL::ACUTE)
      {
        std::cout << "Crashing degenerate motorcycles: " << mc.id() << " and " << fmc.id() << std::endl;
        drive_to_closest_target(mc);
        crash_motorcycle(mc);

        drive_to_closest_target(fmc);
        crash_motorcycle(fmc);

        break;
      }
#endif
    }
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
    const Track_segment_ptr ts_ptr = --(mc.track().end());
    CGAL_postcondition(mc.current_position() == ts_ptr->target());

    add_track_segment_to_map(mc.current_face(), ts_ptr);

    // If we have just added a second track segment, we can remove the first degenerate one
    if(mc.track().size() == 2)
    {
      // @todo
    }
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_foreign_motorcycles(Motorcycle& mc, Collision_information& tc)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "~~~~~~~~~X ?" << std::endl;
  std::cout << "Checking for collisions on motorcycle #" << mc.id() << "'s track"
            << " with foreign faces" << std::endl;
#endif

  // We can look only at collisions with the closest target, except if the whole
  // segment "position -- closest_target" is on the same border halfedge.

  descriptor_variant target_dv = PMP::get_descriptor_from_location(mc.closest_target()->location(), mesh());
  if(const face_descriptor* fd_ptr = boost::get<face_descriptor>(&target_dv))
  {
    // The target is not on the border, thus there's simply nothing to do because
    // we don't care about intersections at the origin.

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " Tentative track's target is not on border" << std::endl;
#endif
    return NO_COLLISION;
  }

  descriptor_variant source_dv = PMP::get_descriptor_from_location(mc.current_location(), mesh());
  if(const face_descriptor* fd_ptr = boost::get<face_descriptor>(&source_dv))
  {
    // Tentative track's source is not on a border.

    // Small skip: if we have already found an intersection strictly within the face,
    // there's no point to check adjacent faces, since the intersection will be
    // at a later time.
    if(tc.time_at_closest_collision < mc.time_at_closest_target())
      return NO_COLLISION;

    return find_collision_with_tentative_track_target_on_border(mc, target_dv, tc);
  }
  else // tentative track's source and closest target are on a border
  {
    // check if source and targets lie on the same halfedge
    halfedge_descriptor hd = halfedge(mc.current_face(), mesh()), done(hd);
    bool are_on_same_halfedge = false;

    do
    {
      if(PMP::is_on_halfedge(mc.current_position()->location(), hd, mesh()) &&
         PMP::is_on_halfedge(mc.closest_target()->location(), hd, mesh()))
      {
        are_on_same_halfedge = true;
        break;
      }

      hd = next(hd, mesh());
    } while(hd != done);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Tentative track on the same halfedge: " << are_on_same_halfedge << std::endl;
#endif

    if(are_on_same_halfedge)
    {
      // same halfedge, means that we must consider the full segment and look
      // for intersections in the opposite face
      return find_foreign_collision_with_tentative_track_on_border(mc, hd, tc);

      if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&target_dv))
      {
        // closest target is on a vertex => need to check the faces incident to 'vd'
        return find_collision_with_tentative_track_target_on_border(mc, target_dv, tc);
      }
    }
    else // not on the same halfedge, only look at the destination
    {
      return find_collision_with_tentative_track_target_on_border(mc, target_dv, tc);
    }
  }
}

// Below, only the target of the tentative track is on a border
// ---------------------------------------------------------------------------------
template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_tentative_track_target_on_border(const Motorcycle& mc,
                                                     const descriptor_variant dv,
                                                     Collision_information& tc) const
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "¤ Find collision with tentative track target of motorcycle #" << mc.id() << " on border" << std::endl;
#endif

  CGAL_expensive_precondition(dv == CGAL::Polygon_mesh_processing::
                              get_descriptor_from_location(mc.closest_target()->location(), mesh()));

  if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&dv))
  {
    const vertex_descriptor vd = *vd_ptr;

    // check all incident faces at 'vd' and intersections at vd
    const halfedge_descriptor hd = halfedge(vd, mesh());
    BOOST_FOREACH(face_descriptor ffd, CGAL::faces_around_target(hd, mesh()))
    {
      if(ffd == mc.current_face() || ffd == boost::graph_traits<Triangle_mesh>::null_face())
        continue;

      return find_collision_with_tentative_track_target_on_border(mc, dv, ffd, tc);
    }
  }
  else // mc's closest target is on a halfedge
  {
    const halfedge_descriptor hd = boost::get<halfedge_descriptor>(dv);

    if(is_border(edge(hd, mesh()), mesh()))
      return NO_COLLISION;

    // check opposite face for intersection at the mc.closest_target()
    const face_descriptor ffd = face(opposite(hd, mesh()), mesh());
    return find_collision_with_tentative_track_target_on_border(mc, dv, ffd, tc);
  }

  return NO_COLLISION;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_tentative_track_target_on_border(const Motorcycle& mc,
                                                     const descriptor_variant dv,
                                                     const face_descriptor ffd,
                                                     Collision_information& tc) const
{
  CGAL_precondition(ffd != boost::graph_traits<Triangle_mesh>::null_face());
  CGAL_precondition(mc.current_face() != ffd);

  Collision_return res = NO_COLLISION;

  // Step 1: check complete tracks
  TFM_const_iterator it = track_face_map_.find(ffd);
  if(it != track_face_map_.end())
  {
    const Track_segment_ptr_container& face_tracks = it->second;

    typename Track_segment_ptr_container::const_iterator tl_it = face_tracks.begin();
    typename Track_segment_ptr_container::const_iterator tl_end = face_tracks.end();
    for(; tl_it!=tl_end; ++tl_it)
    {
      const Track_segment& ts = *(*tl_it);

      Collision_return r = find_collision_with_tentative_track_target_on_border_with_track_on_foreign_face(mc, dv, ts, false /*fmc is not moving*/, tc);

      // Need to keep the foreign track in memory to add a new point on the confirmed track...
      if(r == COLLISION)
      {
        CGAL_assertion(!tc.is_foreign_motorcycle_moving_on_track);
        tc.foreign_track = *tl_it;
      }

      res = res | r;

      if(res == SNAPPED_COLLISION_TO_EXISTING_POINT)
        return res;
    }
  }

  // Step 2: check incomplete tracks (path of a motorcycle currently moving in the same face)
  MCC_cit fmc_it = motorcycles().begin(), fmc_end = motorcycles().end();
  for(; fmc_it!=fmc_end; ++fmc_it)
  {
    const Motorcycle& fmc = motorcycle(fmc_it);
    res = res | find_collision_with_tentative_track_target_on_border_with_live_motorcycle_on_foreign_face(mc, dv, ffd, fmc, tc);

    if(res == SNAPPED_COLLISION_TO_EXISTING_POINT)
      return res;
  }

  return res;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_tentative_track_target_on_border_with_live_motorcycle_on_foreign_face(const Motorcycle& mc,
                                                                                          const descriptor_variant dv,
                                                                                          const face_descriptor ffd,
                                                                                          const Motorcycle& fmc,
                                                                                          Collision_information& tc) const
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE_PLUS
  std::cout << "¤ Checking for foreign intersection with live motorcycle #" << fmc.id() << std::endl;
#endif

  CGAL_precondition(ffd != boost::graph_traits<Triangle_mesh>::null_halfedge());
  CGAL_precondition(mc.current_face() != ffd);

  if(// the foreign motorcycle must be in the foreign face 'ffd'
     fmc.current_face() != ffd ||
     // the foreign motorcycle must be in motion
     fmc.is_crashed())
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE_PLUS
    std::cout << " ignoring 'fmc' in foreign face..." << std::endl;
    std::cout << "  > motorcycles #" << mc.id() << " and #" << fmc.id() << std::endl;
    std::cout << "  > faces: " << fmc.current_face() << " and " << fmc.current_face() << std::endl;
    std::cout << "  > crashed status: " << fmc.is_crashed() << std::endl;
#endif
    return NO_COLLISION;
  }

  CGAL_assertion(fmc.id() != mc.id());

  Track_segment fmc_track(fmc.id(), fmc.current_position(), fmc.current_time(),
                          fmc.closest_target(), fmc.time_at_closest_target());

  return find_collision_with_tentative_track_target_on_border_with_track_on_foreign_face(mc, dv, fmc_track, true /*fmc is not moving*/, tc);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_tentative_track_target_on_border_with_track_on_foreign_face(const Motorcycle& mc,
                                                                                const descriptor_variant ct_dv,
                                                                                const Track_segment& fmc_track,
                                                                                const bool is_fmc_moving_on_track,
                                                                                Collision_information& tc) const
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  const std::size_t fmc_id = fmc_track.motorcycle_id();

  const Motorcycle& fmc = motorcycle(fmc_id);
  const Node_ptr fmc_track_source = fmc_track.source();
  const Node_ptr fmc_track_destination = fmc_track.target();

  const bool is_fmcs_degenerate = (fmc_track_source == fmc_track_destination);

  const face_descriptor ffd = fmc_track_source->face();
  CGAL_assertion(ffd == fmc_track_destination->face());

  const Node_ptr ct = mc.closest_target();
  const Node_ptr ct_in_ffd = ct->sibling(ffd);
  CGAL_postcondition(ffd == ct_in_ffd->face());

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "-----------------------------" << std::endl;
  std::cout << "¤¤ Checking collision with single point on border "
            << "of foreign motorcycle #" << fmc_id << std::endl;
  std::cout << "moving on track: " << is_fmc_moving_on_track << std::endl;
  std::cout << " + closest target: " << &*ct << std::endl << *ct << std::endl;
  std::cout << " + location in foreign face: " << " " << ct_in_ffd->face() << " bc: "
                                               << ct_in_ffd->barycentric_coordinate(0) << " "
                                               << ct_in_ffd->barycentric_coordinate(1) << " "
                                               << ct_in_ffd->barycentric_coordinate(2) << std::endl;
  std::cout << " + source: " << &*fmc_track_source << std::endl << *fmc_track_source << std::endl;
  std::cout << " + target: " << &*fmc_track_destination << std::endl << *fmc_track_destination << std::endl;
#endif

  const FT time_at_collision = mc.time_at_closest_target();
  const FT time_at_fmc_track_source = fmc_track.time_at_source();
  const FT time_at_fmc_track_destination = fmc_track.time_at_target();

  if(is_fmcs_degenerate)
  {
    // the only possible intersection is if mc's target is the same point
    // as the degenerate foreign track
    if(mc.closest_target() == fmc_track_source)
    {
      return find_collision_at_tentative_track_destination(mc, fmc, time_at_fmc_track_source,
                                                           true /*on foreign face*/, is_fmc_moving_on_track, tc);
    }

    return NO_COLLISION;
  }

  FT foreign_visiting_time;
  if(ct->has_motorcycle(fmc.id(), time_at_fmc_track_source,
                        time_at_fmc_track_destination, foreign_visiting_time))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "  /!\\ Tentative path collides with track on foreign face of motorcycle #: " << fmc.id()
              << " at the closest target. Time: " << time_at_collision << std::endl;
#endif

    if(tc.is_collision_earlier_than_current_best(time_at_collision, foreign_visiting_time, is_fmc_moving_on_track))
    {
      tc.reset();
      tc.is_closest_collision_already_in_dictionary = true;
      tc.closest_collision = mc.closest_target();
      tc.time_at_closest_collision = time_at_collision;

      tc.fmc_id = fmc.id();
      tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
      tc.is_foreign_motorcycle_in_different_face = true;
      tc.foreign_motorcycle_face = ffd;
      tc.foreign_time_at_closest_collision = foreign_visiting_time;

      return COLLISION;
    }
  }
  else if(const vertex_descriptor* vd_ptr = boost::get<vertex_descriptor>(&ct_dv))
  {
    // If the closest target is on a vertex_descriptor, then the only possible
    // intersection is with 'fmc_track_source' or 'fmc_track_destination'
    // and it will (should) have been found with the check above if it exists.
    return NO_COLLISION;
  }
  else if(const halfedge_descriptor* hd_ptr = boost::get<halfedge_descriptor>(&ct_dv))
  {
    halfedge_descriptor hd = *hd_ptr;
    // Need to check that the track [fmc_track_source, fmc_track_destination]
    // does not contain mc.closest_target()

    // If the extremities of the foreign track are not on a border halfedge,
    // then there can't be an intersection with a point on the border (except
    // for source or target, which have been checked above)

    // Check if source and targets lie on the same halfedge
    halfedge_descriptor cfhd = halfedge(ffd, mesh()), done(cfhd);
    bool are_on_same_halfedge = false;

    do
    {
      if(PMP::is_on_halfedge(fmc_track_source->location(), cfhd, mesh()) &&
         PMP::is_on_halfedge(fmc_track_destination->location(), cfhd, mesh()))
      {
        are_on_same_halfedge = true;
        break;
      }

      cfhd = next(cfhd, mesh());
    } while(cfhd != done);

    if(!are_on_same_halfedge)
      return NO_COLLISION;

    // 'hd' is in the non-foreign face, and we want the halfedge in the foreign face
    halfedge_descriptor opp_hd = opposite(hd, mesh());

    if(cfhd != opp_hd)
      return NO_COLLISION;

    // We are now in the configuration of 'mc' having a single point on a halfedge,
    // and the foreign track is on the opposite halfedge

    const Point_2 s = geom_traits().construct_point_2_object()(fmc_track_source->location().second[0],
                                                               fmc_track_source->location().second[1]);
    const Point_2 t = geom_traits().construct_point_2_object()(fmc_track_destination->location().second[0],
                                                               fmc_track_destination->location().second[1]);
    const Point_2 ct2 = geom_traits().construct_point_2_object()(ct_in_ffd->barycentric_coordinate(0),
                                                                 ct_in_ffd->barycentric_coordinate(1));

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "s-ct2-t: " << s << " || " << ct2 << " || " << t << std::endl;
#endif

    CGAL_assertion(s != ct2 && t != ct2);

    // Below might fail due to numerical errors, but it is supposed to be 'true'
#ifdef CGAL_POLYLINE_TRACING_ENABLE_RIGOROUS_PRECONDITIONS
    CGAL_assertion(geom_traits().collinear_2_object()(s, ct2, t));
#endif

    // Check if the closest target is in between the source and the target
    if(!geom_traits().collinear_are_strictly_ordered_along_line_2_object()(s, ct2, t))
      return NO_COLLISION;

    // From here on, 'ct2' is strictly in between 's' and 't'

    // No choice but to compute the foreign time
    const FT time_at_fmc_track_source = fmc_track.time_at_source();
    const FT foreign_time_at_collision = time_at_fmc_track_source +
      CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(),
                                         mc.closest_target()->point())) / fmc.speed();

    if(tc.is_collision_earlier_than_current_best(time_at_collision, foreign_time_at_collision, is_fmc_moving_on_track))
    {
      tc.reset();
      tc.is_closest_collision_already_in_dictionary = true;
      tc.closest_collision = mc.closest_target();
      tc.time_at_closest_collision = time_at_collision;

      tc.fmc_id = fmc_id;
      tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
      tc.is_foreign_motorcycle_in_different_face = true;
      tc.foreign_motorcycle_face = ffd;
      tc.foreign_time_at_closest_collision = foreign_time_at_collision;

      return COLLISION;
    }
  }
  else
  {
    // Motorcycle is not moving on a border and we shouldn't be here
    CGAL_assertion(false);
  }

  return NO_COLLISION;
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

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "¤ Checking collision with tentative track on border" << std::endl;
#endif

  const halfedge_descriptor opp_hd = opposite(hd, mesh());
  if(is_border(opp_hd, mesh()))
    return NO_COLLISION;

  Collision_return res = NO_COLLISION;
  const face_descriptor ffd = face(opp_hd, mesh());

  // Step 1: check complete tracks
  TFM_const_iterator it = track_face_map_.find(ffd);
  if(it != track_face_map_.end())
  {
    const Track_segment_ptr_container& face_tracks = it->second;

    typename Track_segment_ptr_container::const_iterator tl_it = face_tracks.begin();
    typename Track_segment_ptr_container::const_iterator tl_end = face_tracks.end();
    for(; tl_it!=tl_end; ++tl_it)
    {
      const Track_segment& ts = *(*tl_it);

      Collision_return r = find_collision_with_track_on_foreign_face(mc, hd, ts, false /*is_fmc_moving_on_track*/, tc);

      // Need to keep the foreign track in memory to add a new point on the confirmed track...
      if(r == COLLISION)
      {
        CGAL_assertion(!tc.is_foreign_motorcycle_moving_on_track);
        tc.foreign_track = *tl_it;
      }

      res = res | r;

      if(res == SNAPPED_COLLISION_TO_EXISTING_POINT)
        return res;
    }
  }

  // Step 2: check incomplete tracks (path of a motorcycle currently moving in the same face)
  MCC_it fmc_it = motorcycles().begin(), fmc_end = motorcycles().end();
  for(; fmc_it!=fmc_end; ++fmc_it)
  {
    Motorcycle& fmc = motorcycle(fmc_it);
    res = res | find_collision_with_live_motorcycle_on_foreign_face(mc, hd, fmc, tc);

    if(res == SNAPPED_COLLISION_TO_EXISTING_POINT)
      return res;
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

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "¤ Checking for foreign intersection with live motorcycle #" << fmc.id()
            << " in foreign face: " << ffd << std::endl;
#endif

  CGAL_precondition(!is_border(edge(hd, mesh()), mesh()));
  CGAL_precondition(mc.current_face() != ffd);

  if(// the foreign motorcycle must be in the foreign face 'ffd'
     fmc.current_face() != ffd ||
     // the foreign motorcycle must be in motion
     fmc.is_crashed())
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << " ignoring 'fmc' in foreign face... " << std::endl;
    std::cout << "  > motorcycles #" << mc.id() << " and #" << fmc.id() << std::endl;
    std::cout << "  > faces: " << mc.current_face() << " and " << fmc.current_face() << std::endl;
    std::cout << "  > crashed status: " << fmc.is_crashed() << std::endl;
#endif
    return NO_COLLISION;
  }

  CGAL_assertion(fmc.id() != mc.id());

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

  const std::size_t fmc_id = fmc_track.motorcycle_id();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "¤¤ Checking collision with tentative track on border "
            << "and foreign motorcycle #" << fmc_id << std::endl;
#endif

  CGAL_precondition(!is_border(edge(hd, mesh()), mesh()));

  const Motorcycle& fmc = motorcycle(fmc_id);
  const Node_ptr fmc_track_source = fmc_track.source();
  const FT time_at_fmc_track_source = fmc_track.time_at_source();
  const Node_ptr fmc_track_destination = fmc_track.target();
  const FT time_at_fmc_track_destination = fmc_track.time_at_target();

  const halfedge_descriptor opp_hd = opposite(hd, mesh());

  bool is_fts_on_opp_hd = PMP::is_on_halfedge(fmc_track_source->location(), opp_hd, mesh());
  bool is_ftd_on_opp_hd = PMP::is_on_halfedge(fmc_track_destination->location(), opp_hd, mesh());

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
                                                         time_at_fmc_track_source,
                                                         is_fmc_moving_on_track, tc);
    }
  }
  else if(is_ftd_on_opp_hd) // !is_fts_on_opp_hd && is_ftd_on_opp_hd
  {
    // only possible intersection is at the destination
    return find_collision_with_foreign_track_extremity(mc, hd, fmc, fmc_track_destination,
                                                       time_at_fmc_track_destination,
                                                       is_fmc_moving_on_track, tc);
  }

  return NO_COLLISION;
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

  const std::size_t fmc_id = fmc_track.motorcycle_id();
  const Motorcycle& fmc = motorcycle(fmc_id);
  const Node_ptr fmc_track_source = fmc_track.source();
  const Node_ptr fmc_track_destination = fmc_track.target();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "¤¤¤ Find collision between collinear tracks of motorcycles #"
            << mc.id() << " and #" << fmc.id() << std::endl;
  std::cout << "   foreign track: " << std::endl << *fmc_track_source << std::endl
                                                 << *fmc_track_destination << std::endl;
#endif

  CGAL_precondition(PMP::is_on_halfedge(mc.current_position()->location(), hd, mesh()));
  CGAL_precondition(PMP::is_on_halfedge(mc.closest_target()->location(), hd, mesh()));

  const halfedge_descriptor opp_hd = opposite(hd, mesh());
  CGAL_precondition(!is_border(opp_hd, mesh()));
  const face_descriptor ffd = face(opp_hd, mesh());

  const Face_location& cp_in_ffd = mc.current_position()->sibling(ffd)->location();
  const Face_location& ct_in_ffd = mc.closest_target()->sibling(ffd)->location();

  const Point_2 s = geom_traits().construct_point_2_object()(cp_in_ffd.second[0], cp_in_ffd.second[1]);
  const Point_2 t = geom_traits().construct_point_2_object()(ct_in_ffd.second[0], ct_in_ffd.second[1]);
  const Segment_2 mcs = geom_traits().construct_segment_2_object()(s, t);

  const Point_2 fs = geom_traits().construct_point_2_object()(fmc_track_source->location().second[0],
                                                              fmc_track_source->location().second[1]);
  const Point_2 ft = geom_traits().construct_point_2_object()(fmc_track_destination->location().second[0],
                                                              fmc_track_destination->location().second[1]);
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
                                            const bool is_fmc_moving_on_track,
                                            Collision_information& tc) const
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  // this is the case of 'mc' tentative track being on a border, and a foreign
  // track with a single point on this same border

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "¤¤¤ Checking collision with tentative track on border"
            << " with foreign motorcycle " << fmc.id()
            << " and single foreign point on border: " << std::endl;
#endif

  // mc's track is non-degenerate
  CGAL_precondition(mc.current_position() != mc.closest_target());
  // mc's track in on the halfedge
  CGAL_precondition(PMP::is_on_halfedge(mc.current_position()->location(), hd, mesh()));
  CGAL_precondition(PMP::is_on_halfedge(mc.closest_target()->location(), hd, mesh()));
  // the foreign extremity is on a halfedge
  CGAL_precondition(PMP::is_on_face_border(foreign_extremity->location(), mesh()));

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
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
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "s == e" << std::endl;
#endif
    return NO_COLLISION;
  }
  else if(t == e) // intersection at mc's closest target
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "t == e" << std::endl;
#endif
    const FT time_at_collision = mc.time_at_closest_target();

    // Compare to current tentative collision to keep the closest intersection
    if(tc.is_collision_earlier_than_current_best(time_at_collision, foreign_time_at_collision, is_fmc_moving_on_track))
    {
      tc.reset();
      tc.is_closest_collision_already_in_dictionary = true;
      tc.closest_collision = mc.closest_target();
      tc.time_at_closest_collision = time_at_collision;

      tc.fmc_id = fmc.id();
      tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
      tc.is_foreign_motorcycle_in_different_face = true;
      tc.foreign_motorcycle_face = ffd;
      tc.foreign_time_at_closest_collision = foreign_time_at_collision;

      return COLLISION;
    }
  }
  else // general case
  {
    // The assertion below might fail due to numerical errors, but it is, logically,
    // a correct statement (case of three points on the same halfedge)
#ifdef CGAL_POLYLINE_TRACING_ENABLE_RIGOROUS_PRECONDITIONS
      CGAL_assertion(geom_traits().collinear_2_object()(s, e, t));
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "Intersection not at source or target" << std::endl;
#endif

    if(!geom_traits().collinear_are_strictly_ordered_along_line_2_object()(s, e, t))
      return NO_COLLISION;

    // From here on, e is on ]s;t[
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "e is on ]s;t[" << std::endl;
#endif

    const Point collision_point = foreign_extremity->point();
    const FT time_at_collision = mc.current_time() +
      CGAL::sqrt(CGAL::squared_distance(mc.current_position()->point(),
                                        collision_point)) / mc.speed();

    if(tc.is_collision_earlier_than_current_best(time_at_collision, foreign_time_at_collision, is_fmc_moving_on_track))
    {
      tc.reset();
      tc.is_closest_collision_already_in_dictionary = true;
      tc.closest_collision = foreign_extremity;
      tc.time_at_closest_collision = time_at_collision;

      tc.fmc_id = fmc.id();
      tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
      tc.is_foreign_motorcycle_in_different_face = true;
      tc.foreign_motorcycle_face = ffd;
      tc.foreign_time_at_closest_collision = foreign_time_at_collision;

      return COLLISION;
    }
  }

  return NO_COLLISION;
}

// collisions between two motorcycles in the same face
template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_at_tentative_track_destination(const Motorcycle& mc,
                                              const Motorcycle& fmc,
                                              const FT fmc_visiting_time,
                                              const bool is_fmc_on_foreign_face,
                                              const bool is_fmc_moving_on_track,
                                              Collision_information& tc) const
{
  FT time_at_collision = mc.time_at_closest_target();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  /!\\ Tentative path collides with track: " << fmc.id()
            << " at the closest target. Time: " << time_at_collision << std::endl;
#endif

  if(tc.is_collision_earlier_than_current_best(time_at_collision, fmc_visiting_time, is_fmc_moving_on_track))
  {
    tc.reset();
    tc.is_closest_collision_already_in_dictionary = true;
    tc.closest_collision = mc.closest_target();
    tc.time_at_closest_collision = time_at_collision;

    tc.fmc_id = fmc.id();
    tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
    tc.is_foreign_motorcycle_in_different_face = is_fmc_on_foreign_face;
    tc.foreign_time_at_closest_collision = fmc_visiting_time;

    return COLLISION;
  }

  return NO_COLLISION;
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
  // Below might fail due to numerical errors, but we are treating here the
  // case of two collinear tracks, possibly on two different faces incident to
  // the same edge.
#ifdef CGAL_POLYLINE_TRACING_ENABLE_RIGOROUS_PRECONDITIONS
  CGAL_precondition(geom_traits().collinear_2_object()(mcs.source(), fmcs.source(), mcs.target()));
  CGAL_precondition(geom_traits().collinear_2_object()(mcs.source(), fmcs.target(), mcs.target()));
#endif

  // Many different configurations exist, e.g. (_S is for source, _T for target):
  //  MC_S  ---- FMC_S ---- FMC_T ---- MC_T
  //  FMC_T ---- MC_S  ---- FMC_S ---- MC_T
  // etc.
  // If, on the ray MC_S->MC_T,
  // - FMC_S is "before" MC_S, then it doesn't matter for MC whichever respective
  //   direction the motorcycles are moving in.
  // - FMC_S is MC_S, then it only matters if they are moving in the same direction
  //   but this already treated before the algorithm starts, in the function
  //   'crash_motorcycles_with_same_sources_and_directions()'
  // - FMC_S is "after" MC_S, then it depends on the motorcycles' directions.

  if(mcs.source() == fmcs.source())
    return NO_COLLISION;

  bool is_fmcs_degenerate = geom_traits().is_degenerate_2_object()(fmcs);

  // Compute the respective direction of the two motorcycles:
  CGAL_precondition(mcs.source() != mcs.target());
  bool are_motorcycles_moving_in_the_same_direction =
    (is_fmcs_degenerate ||
     geom_traits().angle_2_object()(mcs.source(), mcs.target(), fmcs.source(), fmcs.target()) == CGAL::ACUTE);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  is degen: " << is_fmcs_degenerate << std::endl;
  std::cout << "  angle: " << geom_traits().angle_2_object()(mcs.source(), mcs.target(),
                                                             fmcs.source(), fmcs.target()) << std::endl;
  std::cout << "  are motorcycles moving in the same direction: "
            << are_motorcycles_moving_in_the_same_direction << std::endl;
#endif

  FT time_at_collision = 0.;
  const Node_ptr fmc_track_source = fmc_track.source();
  const FT time_at_fmc_track_source = fmc_track.time_at_source();
  const Node_ptr fmc_track_destination = fmc_track.target();
  const FT time_at_fmc_track_destination = fmc_track.time_at_target();
  const face_descriptor ffd = fmc_track_source->face();
  CGAL_assertion(ffd == fmc_track_destination->face());

  const bool are_motorcycles_on_the_same_face =
    (mc.current_face() == fmc_track_source->face());

  // Some sanity checks -----
  CGAL_assertion(fmc_track_source->face() == fmc_track_destination->face());
  CGAL_assertion(time_at_fmc_track_source <= time_at_fmc_track_destination);

  if(!are_motorcycles_on_the_same_face)
  {
    // Check that all the track points are indeed on the same halfedge
    CGAL_precondition_code
    (
      boost::optional<halfedge_descriptor> hd =
        CGAL::Polygon_mesh_processing::internal::common_halfedge(fmc_track_source->face(),
                                                                 mc.current_face(),
                                                                 mesh());
    )
    CGAL_precondition(bool(hd));
    CGAL_precondition_code(halfedge_descriptor opp_hd = opposite(*hd, mesh());)
    CGAL_precondition(CGAL::Polygon_mesh_processing::is_on_halfedge(fmc_track_source->location(), *hd, mesh()));
    CGAL_precondition(CGAL::Polygon_mesh_processing::is_on_halfedge(fmc_track_destination->location(), *hd, mesh()));
    CGAL_precondition(CGAL::Polygon_mesh_processing::is_on_halfedge(mc.current_position()->location(), opp_hd, mesh()));
    CGAL_precondition(CGAL::Polygon_mesh_processing::is_on_halfedge(mc.closest_target()->location(), opp_hd, mesh()));
  }
  // end of sanity checks -----

  // The motorcycles move in the same direction
  if(are_motorcycles_moving_in_the_same_direction)
  {
    // If there's an intersection, 'mc' will impact fmcs' source.

    // The weird configuration of both motorcycles moving in the same direction
    // AND with the same source is handled by crashing motorcycles at the very
    // beginning, see function: 'crash_motorcycles_with_same_sources_and_directions()'
    CGAL_assertion(is_fmcs_degenerate || mcs.source() != fmcs.source());

    if(mcs.target() == fmcs.source())
    {
      time_at_collision = mc.time_at_closest_target();
    }
    // Note that here, we know that fmcs.source() != mcs.source() and mcs.target()
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
      return NO_COLLISION;
    }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "  Motorcycles #" << mc.id() << " crashes into the source of Motorcycle #"
                                   << fmc.id() << " at time: " << time_at_collision << std::endl;
#endif

    if(tc.is_collision_earlier_than_current_best(time_at_collision, time_at_fmc_track_source, is_fmc_moving_on_track))
    {
      tc.reset();
      tc.is_closest_collision_already_in_dictionary = true;
      tc.closest_collision = fmc_track_source;
      tc.time_at_closest_collision = time_at_collision;

      tc.fmc_id = fmc.id();
      tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
      tc.is_foreign_motorcycle_in_different_face = !are_motorcycles_on_the_same_face;
      tc.foreign_motorcycle_face = ffd;
      tc.foreign_time_at_closest_collision = time_at_fmc_track_source;
      CGAL_assertion(!tc.is_foreign_motorcycle_in_different_face || mc.current_face() != ffd);

      return COLLISION;
    }
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
      return NO_COLLISION;

    // If mc's target is in [mcs, fmcs], then there is no intersection
    if(mcs.target() != fmcs.source() && // to be able to use strictly (and there is an intersection if 'true')
       mcs.target() != fmcs.target() &&
       geom_traits().collinear_are_strictly_ordered_along_line_2_object()(mcs.target(),
                                                                          fmcs.target(),
                                                                          fmcs.source()))
      return NO_COLLISION;

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
                                            fmc_track_destination->point())) / mc.speed();

        // @todo time snapping here ?

        CGAL_assertion(!mc.has_target_at_time(time_at_collision).second);
      }
      else
      {
        // fmcs.target() can't be 'before' mcs.source() because 'not_moving' means
        // that we are on a confirmed track and if fmcs.target() is 'after' mcs.target(),
        // then there is no intersection.
        return NO_COLLISION;
      }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "  Motorcycles #" << mc.id() << " crashes into the final position of Motorcycle #"
              << fmc.id() << " at time: " << time_at_collision << std::endl;
#endif

      if(tc.is_collision_earlier_than_current_best(time_at_collision, time_at_fmc_track_destination, is_fmc_moving_on_track))
      {
        tc.reset();
        tc.is_closest_collision_already_in_dictionary = true;
        tc.closest_collision = fmc_track_destination;
        tc.time_at_closest_collision = time_at_collision;

        tc.fmc_id = fmc.id();
        tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
        tc.foreign_time_at_closest_collision = time_at_fmc_track_destination;
        tc.is_foreign_motorcycle_in_different_face = !are_motorcycles_on_the_same_face;
        tc.foreign_motorcycle_face = ffd;
        CGAL_assertion(!tc.is_foreign_motorcycle_in_different_face ||
                       mc.current_face() != ffd);

        return COLLISION;
      }
    }
    else // The foreign motorcycle is (also) moving
    {
      // The collision is at the middle point and both motorcycles reach it at the same time.
      // Note that this point might not actually be reached by either motorcycle,
      // e.g. if a motorcycle crashes before reaching it.

      // @todo, if speeds are ever allowed to change while tracing, the speed
      // of fmc here must be changed to the speed on the track segment 'fmc_track'
      const FT sqd = CGAL::squared_distance(mc.current_position()->point(),
                                            fmc_track_source->point());
      time_at_collision = mc.current_time() +
        (CGAL::sqrt(sqd) - fmc.speed() * (mc.current_time() - time_at_fmc_track_source)) / (mc.speed() + fmc.speed());

      // @todo time snapping here

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "  sqd: " << sqd << std::endl;
      std::cout << "  speeds: " << mc.speed() << " " << fmc.speed() << std::endl;
      std::cout << "  current times: " << mc.current_time() << " " << time_at_fmc_track_source << std::endl;
      std::cout << "  final time: " << time_at_collision << std::endl;
      std::cout << "  § mc and fmc would meet at time: " << time_at_collision << std::endl;
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
      // By construction, the time and foreign_time should be greater
      // than the times at the sources of the tracks. Some numerical errors
      // can sneak it and, if so, snap the time.
      //
      // It should only be a numerical error, that is a very small error
      if(time_at_collision < mc.current_time())
      {
        CGAL_precondition(time_at_collision + tolerance_ >= mc.current_time());
        time_at_collision = mc.current_time();
      }

      if(time_at_collision < time_at_fmc_track_source)
      {
        CGAL_precondition(time_at_collision + tolerance_ >= time_at_fmc_track_source);
        time_at_collision = time_at_fmc_track_source;
      }
#endif

      CGAL_postcondition(time_at_collision >= time_at_fmc_track_source);
      CGAL_postcondition(time_at_collision >= mc.current_time());

      if(tc.is_collision_earlier_than_current_best(time_at_collision, time_at_collision, is_fmc_moving_on_track))
      {
        // both values are used later when we snap times/points
        const FT time_at_closest_collision_memory = tc.time_at_closest_collision;
        const FT foreign_time_at_closest_collision_memory = tc.foreign_time_at_closest_collision;

        tc.reset();
        tc.time_at_closest_collision = time_at_collision;

        tc.fmc_id = fmc.id();
        tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
        tc.is_foreign_motorcycle_in_different_face = !are_motorcycles_on_the_same_face;
        tc.foreign_motorcycle_face = ffd;
        tc.foreign_time_at_closest_collision = time_at_collision;

        // Temporal snapping ---------------------------------------------------
        // Try to find the collision point by checking if any of the motorcycles
        // has a point at that time.
        std::pair<TPC_iterator, bool> mc_res = mc.has_target_at_time(time_at_collision);
        if(mc_res.second) // there is already a target at that time
        {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
          std::cout << "Motorcycle #" << mc.id() << " already has a target at time: " << time_at_collision << std::endl;
#endif

          TPC_iterator target_point = mc_res.first;
          CGAL_assertion(target_point->second == time_at_collision);
          Node_ptr alternate_collision = target_point->first;

          tc.is_closest_collision_already_in_dictionary = true;
          tc.closest_collision = alternate_collision;

          return COLLISION;
        }

        // Same check, but with the foreign time at collision
        std::pair<TPC_iterator, bool> fmc_res = fmc.has_target_at_time(time_at_collision);
        if(fmc_res.second) // there is already a target at that time
        {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
          std::cout << "Foreign motorcycle #" << fmc.id() << " already has a target at time: " << time_at_collision << std::endl;
#endif

          TPC_iterator target_point = fmc_res.first;
          Node_ptr alternate_foreign_collision = target_point->first;
          CGAL_assertion(alternate_foreign_collision->face() == fmc.current_face());
          CGAL_assertion(target_point->second == time_at_collision);

          tc.is_closest_collision_already_in_dictionary = true;
          tc.closest_collision = alternate_foreign_collision;

          return COLLISION;
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

          tc.is_closest_collision_already_in_dictionary = true;
          tc.closest_collision = collision_entry.first;

          // We previously searched by time but couldn't find anything but the
          // point existed. Check if that point is visited by either 'mc' or 'fmc';
          // if it's the case, we need to repare the time to be that of the existing
          // point.

          // Add a small tolerance on the time since we previously didn't find any target at the exact time
          FT visiting_time; // will be filled by the call to 'has_motorcycle'
          if(collision_entry.first->has_motorcycle(mc.id(), time_at_collision - tolerance_,
                                                   time_at_collision + tolerance_, visiting_time))
          {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
            std::cout << "Motorcycle #" << mc.id() << " already has a target at time: " << visiting_time << std::endl;
#endif

            // Assert that we are still the closest collision (not sure what to do otherwise)
            CGAL_assertion(visiting_time < time_at_closest_collision_memory);

            tc.time_at_closest_collision = visiting_time;
            tc.foreign_time_at_closest_collision = visiting_time; // times are equal in this configuration

            return COLLISION;
          }

          // Try with 'fmc'
          FT foreign_visiting_time;
          if(collision_entry.first->has_motorcycle(fmc.id(), time_at_collision - tolerance_,
                                                   time_at_collision + tolerance_, foreign_visiting_time))
          {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
            std::cout << "Foreign motorcycle #" << fmc.id() << " already has a target at time: " << foreign_visiting_time << std::endl;
#endif

            // Assert that we are still the closest collision (not sure what to do otherwise)
            CGAL_assertion_code(if(tc.time_at_closest_collision == time_at_closest_collision_memory))
                CGAL_assertion(foreign_visiting_time < foreign_time_at_closest_collision_memory);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
            std::cout << "found: fmc.id(): " << fmc.id() << " in pt: " << std::endl << *(collision_entry.first) << std::endl;
            std::cout << "foreign_visiting_time: " << foreign_visiting_time << std::endl;
#endif
            tc.foreign_time_at_closest_collision = foreign_visiting_time;
            tc.time_at_closest_collision = foreign_visiting_time; // times are equal in this configuration

            return COLLISION;
          }
        }
        else
        {
          // At this point, we have a new location at an unknown time...
#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
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
            CGAL_assertion(visiting_time < time_at_closest_collision_memory ||
                           visiting_time < foreign_time_at_closest_collision_memory);

            tc.is_closest_collision_already_in_dictionary = true;
            tc.closest_collision = is_snappable.first;
            tc.time_at_closest_collision = visiting_time;
            tc.foreign_time_at_closest_collision = visiting_time;

            return SNAPPED_COLLISION_TO_EXISTING_POINT;
          }
#endif

          // Couldn't snap to anything, 'collision_location' is definitely a new point
          tc.is_closest_collision_already_in_dictionary = false;
          tc.closest_collision_location = collision_location;

          return COLLISION;
        }
      }
    }
  }

  return NO_COLLISION;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_between_tracks(const Motorcycle& mc,
                              const Segment_2& mcs,
                              const Motorcycle& fmc,
                              const Track_segment& fmc_track,
                              const bool is_fmc_moving_on_track,
                              Collision_information& tc) const
{
  // Non degenerate mc segment
  CGAL_precondition(mc.current_position() != mc.closest_target());
  CGAL_precondition(mcs.source() != mcs.target());

  const Node_ptr fmc_track_source = fmc_track.source();
  const FT time_at_fmc_track_source = fmc_track.time_at_source();
  const Node_ptr fmc_track_destination = fmc_track.target();
  const FT time_at_fmc_track_destination = fmc_track.time_at_target();

  // Both tracks must be on the same face
  CGAL_precondition(fmc_track_source->face() == fmc_track_destination->face());

  const Point_2 s = geom_traits().construct_point_2_object()(fmc_track_source->location().second[0],
                                                             fmc_track_source->location().second[1]);
  const Point_2 t = geom_traits().construct_point_2_object()(fmc_track_destination->location().second[0],
                                                             fmc_track_destination->location().second[1]);
  const Segment_2 fmcs = geom_traits().construct_segment_2_object()(s, t);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "¤¤ Checking collision with track of motorcycle #" << fmc.id() << std::endl;
  std::cout << " + source: " << &*fmc_track_source << std::endl << *fmc_track_source << std::endl;
  std::cout << " + target: " << &*fmc_track_destination << std::endl << *fmc_track_destination << std::endl;
#endif

  // Ignore the case of a degenerate fmc track starting at the same source as mc's
  bool is_fmcs_degenerate = geom_traits().is_degenerate_2_object()(fmcs);
  if(is_fmcs_degenerate)
  {
    if(mcs.source() == fmcs.source())
    {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "degenerate fmc and mcs.source() == fmcs.source()" << std::endl;
#endif
      return NO_COLLISION;
    }
#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
    else if(internal::are_logically_collinear_on_border<Geom_traits>(
              mc.current_position()->location(), fmc_track_source->location(), mc.closest_target()->location()))
    {
      return find_collision_between_collinear_tracks(mc, mcs, fmc, fmc_track, fmcs,
                                                     is_fmc_moving_on_track, tc);
    }
#endif
  }

  // Detect whether the motorcycles share the same supporting line.
  // Note that we know that 'mcs' is not degenerate.
  // @todo should this have a tolerance ?
  if(geom_traits().collinear_2_object()(mcs.source(), mcs.target(), fmcs.source()) &&
     geom_traits().collinear_2_object()(mcs.source(), mcs.target(), fmcs.target()))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "  /!\\ Tracks are aligned" << std::endl;
#endif
    return find_collision_between_collinear_tracks(mc, mcs, fmc, fmc_track, fmcs,
                                                   is_fmc_moving_on_track, tc);
  }

  // --- From here on, the tracks are not collinear ---

  // @todo below to another function
  // @todo move all easy exit checks to their own functions

  // Below are a bunch of checks to branch out easily without computing an explicit
  // intersection.
  // - #1: Check if the current position of mc is a known intersection with the foreign track
  // - #2: Check if the closest target of mc is a known intersection with the foreign track
  // - #3: Robustness for intersections on halfedge

  // Check #1: known collision at current_position
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  check #1: motorcycle #" << fmc.id() << " between "
            << time_at_fmc_track_source << " " << time_at_fmc_track_destination << std::endl;
#endif
  if(mc.current_position()->has_motorcycle(fmc.id(), time_at_fmc_track_source, time_at_fmc_track_destination))
  {
    // Ignore this intersection: since we are seeking collision in the tentative track,
    // it means that the position was not blocked
    return NO_COLLISION;
  }

  // Check #2: known collision at closest_target
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  check #2: collisition at tentative's track destination ?" << std::endl;
#endif
  FT foreign_visiting_time; // will be filled by 'has_motorcycle' if fmc visits 'closest_target'
  if(mc.closest_target()->has_motorcycle(fmc.id(), time_at_fmc_track_source,
                                         time_at_fmc_track_destination, foreign_visiting_time))
  {
    return find_collision_at_tentative_track_destination(mc, fmc, foreign_visiting_time,
                                                         false /*same face*/, is_fmc_moving_on_track, tc);
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  // Check #3: collision at destination, with foreign track on an edge
  // Catch some annoying numerical issue: the configuration of FMCS on a halfedge
  // and the motorcycle destination on the same edge (but somehow, do_intersect_2()
  // does not find it...).
  // Only doing it for the closest_target because we don't care about the source.
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  check #3: foreign track and target on the same border" << std::endl;
#endif
  CGAL_assertion(fmc_track_source != mc.closest_target());
  CGAL_assertion(fmc_track_destination != mc.closest_target());

  if(!is_fmcs_degenerate &&
     internal::are_logically_collinear_on_border<Geom_traits>(
      fmc_track_source->location(), mc.closest_target()->location(), fmc_track_destination->location()))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "  foreign track and target are logically collinear on border" << std::endl;
#endif

    if(geom_traits().collinear_are_strictly_ordered_along_line_2_object()(s, mcs.target(), t))
    {
      const FT time_at_collision = mc.time_at_closest_target();
      const FT foreign_time_at_collision = time_at_fmc_track_source +
        CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(),
                                          mc.closest_target()->point())) / fmc.speed();

      if(tc.is_collision_earlier_than_current_best(time_at_collision, foreign_time_at_collision, is_fmc_moving_on_track))
      {
        tc.reset();
        tc.is_closest_collision_already_in_dictionary = true;
        tc.closest_collision = mc.closest_target();
        tc.time_at_closest_collision = time_at_collision;

        tc.fmc_id = fmc.id();
        tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
        tc.foreign_time_at_closest_collision = foreign_time_at_collision;

        return COLLISION;
      }
    }

    return NO_COLLISION;
  }

  // Check #4: collision at foreign destination, with track and foreign destination
  // on the same halfedge.
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  check #4: track and foreign destination on a same halfedge" << std::endl;
#endif
  CGAL_assertion(fmc_track_destination != mc.current_position());
  CGAL_assertion(fmc_track_destination != mc.closest_target());
  if(internal::are_logically_collinear_on_border<Geom_traits>(
      fmc_track_destination->location(), mc.closest_target()->location(), mc.current_location()))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "  track and foreign target are logically collinear on border" << std::endl;
#endif

    if(geom_traits().collinear_are_strictly_ordered_along_line_2_object()(mcs.source(), t, mcs.target()))
    {
      const FT sqd = CGAL::squared_distance(mc.current_position()->point(),
                                            fmc_track_destination->point());
      const FT time_at_collision = mc.current_time() + CGAL::sqrt(sqd) / mc.speed();
      const FT foreign_time_at_collision = time_at_fmc_track_destination;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "  foreign target in ] track [ " << std::endl;
      std::cout << "  Pts: (" << mc.current_position()->point() << ") -- ("
                              << fmc_track_destination->point() << ")" << std::endl;
      std::cout << "  current time: " << mc.current_time() << std::endl;
      std::cout << "  sqd: " << sqd << std::endl;
      std::cout << "  time at collision: " << time_at_collision << std::endl;
#endif

      if(tc.is_collision_earlier_than_current_best(time_at_collision, foreign_time_at_collision, is_fmc_moving_on_track))
      {
        tc.reset();
        tc.is_closest_collision_already_in_dictionary = true;
        tc.closest_collision = fmc_track_destination;
        tc.time_at_closest_collision = time_at_collision;

        tc.fmc_id = fmc.id();
        tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
        tc.foreign_time_at_closest_collision = foreign_time_at_collision;

        return COLLISION;
      }
    }

    return NO_COLLISION;
  }

  // Check #4bis: collision at foreign source, with track and foreign source
  // on the same halfedge.
  CGAL_assertion(fmc_track_source != mc.current_position());
  CGAL_assertion(fmc_track_source != mc.closest_target());
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  check #4: track and foreign source on a same halfedge" << std::endl;
#endif
  if(internal::are_logically_collinear_on_border<Geom_traits>(
      fmc_track_source->location(), mc.closest_target()->location(), mc.current_location()))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "  track and foreign source are logically collinear on border" << std::endl;
#endif

    if(geom_traits().collinear_are_strictly_ordered_along_line_2_object()(mcs.source(), s, mcs.target()))
    {
      const FT time_at_collision = mc.current_time() +
        CGAL::sqrt(CGAL::squared_distance(mc.current_position()->point(),
                                          fmc_track_source->point())) / mc.speed();
      const FT foreign_time_at_collision = time_at_fmc_track_source;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "  foreign source in ] track [, "
                << "time at collision: " << time_at_collision << std::endl;
#endif

      if(tc.is_collision_earlier_than_current_best(time_at_collision, foreign_time_at_collision, is_fmc_moving_on_track))
      {
        tc.reset();
        tc.is_closest_collision_already_in_dictionary = true;
        tc.closest_collision = fmc_track_source;
        tc.time_at_closest_collision = time_at_collision;

        tc.fmc_id = fmc.id();
        tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
        tc.foreign_time_at_closest_collision = foreign_time_at_collision;

        return COLLISION;
      }
    }

    return NO_COLLISION;
  }
#endif

  // --- The general case: the intersection must be computed ---
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  general case..." << std::endl;
#endif

  // Ignoring the case of a degenerate fmcs because if there is an intersection,
  // it will have been caught by the first part of that function,
  // branching: "collinear > moving in the same direction"
  if(is_fmcs_degenerate)
  {
    CGAL_assertion(!geom_traits().do_intersect_2_object()(mcs, fmcs));
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "  No intersection with degenerate fmcs track" << std::endl;
#endif
    return NO_COLLISION;
  }

  if(!geom_traits().do_intersect_2_object()(mcs, fmcs))
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "  No intersection (general case)" << std::endl;
#endif
    return NO_COLLISION;
  }

  // Here is what happens next:
  // 1. we compute the intersection
  // 2. we check if that new location is (EXACTLY) an existing node. If it is, we compute
  //    the visiting times and return the appropriate result (COLLISION, MUST_VISIT, etc.).
  // 3. we compute visiting times and check if there are already existing nodes at
  //    these times on the trajectories of 'mc' and 'fmc'. If there is, we snap
  //    to the position and return the appropriate result.
  // 4. we check if there is an existing node that is close to the collision point.
  //    If there is, we snap to that position, compute the visiting times and return
  //    the appropriate result.
  // 5. return the new collision point

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

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "  /!\\ collision between motorcycles #" << mc.id() << " and #" << fmc.id() << std::endl;
  std::cout << "Collision location: " << collision_location.first << " bc: "
            << collision_location.second[0] << " " << collision_location.second[1] << " " << collision_location.second[2] << std::endl;
#endif

  // Step 2: although we might not have known that these two tracks do intersect,
  // their intersection might be a point that has already been used
  std::pair<Node_ptr, bool> is_already_in_dictionary = nodes().find(collision_location);
  if(is_already_in_dictionary.second)
  {
    Node_ptr collision_point = is_already_in_dictionary.first;
    FT time_at_collision = 0.;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Already in the dictionary at: " << &*collision_point << std::endl << *collision_point << std::endl;
#endif

    // Check if 'mc' already visits the known collision point
    if(collision_point == mc.closest_target())
    {
      time_at_collision = mc.time_at_closest_target();
    }
    else if(collision_point == mc.current_position())
    {
      time_at_collision = mc.current_time();
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

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
      // Although we have found an _existing_ point at the location of the intersection,
      // this point was neither the source or the closest target of 'mc'.
      // Global snapping makes sure that points are not too close from one another.
      // Consequently, the times should be different.
      CGAL_assertion(time_at_collision != mc.current_time());
      CGAL_assertion(time_at_collision != mc.time_at_closest_target());
#endif
    }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "time_at_collision: " << time_at_collision
              << " (closest is: " << tc.time_at_closest_collision << ") " << std::endl;
#endif

    // Partial test of "is_collision_earlier..." to branch out early
    if(time_at_collision <= tc.time_at_closest_collision)
    {
      // Check if 'fmc' already visits the known collision point
      FT foreign_time_at_collision;
      if(collision_point->has_motorcycle(fmc.id(), time_at_fmc_track_source,
                                         time_at_fmc_track_destination, foreign_time_at_collision))
      {
        // The collision point is visited by 'fmc' at time 'foreign_time_at_collision'
      }
      else // 'collision_point' is a known point but has not (yet) been visited by 'fmc'
      {
        // No choice but to compute the foreign visiting time
        const FT sqd = CGAL::squared_distance(fmc_track_source->point(),
                                              collision_point->point());
        foreign_time_at_collision = time_at_fmc_track_source + CGAL::sqrt(sqd) / fmc.speed();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "  Gotta compute the foreign time " << std::endl;
      std::cout << "  Pts: (" << fmc_track_source->point() << ") -- ("
                              << collision_point->point() << ")" << std::endl;
      std::cout << "  foreign source time: " << time_at_fmc_track_source << std::endl;
      std::cout << "  sqd: " << sqd << std::endl;
      std::cout << "  foreign time at collision: " << foreign_time_at_collision << std::endl;
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
        // Although we have found an _existing_ point at the location of the intersection,
        // this point was neither the source or the closest target of 'mc'.
        // Global snapping makes sure that points are not too close from one another.
        // Consequently, the times should be different.
        CGAL_assertion(!fmc.has_target_at_time(tc.foreign_time_at_closest_collision).second);
#endif
      }

      if(tc.is_collision_earlier_than_current_best(time_at_collision, foreign_time_at_collision, is_fmc_moving_on_track))
      {
        tc.reset();
        tc.is_closest_collision_already_in_dictionary = true;
        tc.closest_collision = collision_point;
        tc.time_at_closest_collision = time_at_collision;

        tc.fmc_id = fmc.id();
        tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
        tc.foreign_time_at_closest_collision = foreign_time_at_collision;

        return COLLISION;
      }
    }
  }
  else // The collision location has never been seen before!
  {
    Point collision_point = CGAL::Polygon_mesh_processing::location_to_point(collision_location, mesh());

    FT time_at_collision = mc.current_time() +
      CGAL::sqrt(CGAL::squared_distance(mc.current_position()->point(), collision_point)) / mc.speed();
    FT foreign_time_at_collision = time_at_fmc_track_source +
      CGAL::sqrt(CGAL::squared_distance(fmc_track_source->point(), collision_point)) / fmc.speed();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Location never seen before, corresponds to point ("
              << collision_point << ") at time: " << time_at_collision << std::endl;
    std::cout << "time bounds: " << mc.current_time() << " || "
                                 << mc.time_at_closest_target() << std::endl;
    std::cout << "foreign time bounds: " << time_at_fmc_track_source << " || "
                                         << time_at_fmc_track_destination << std::endl;
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
    // By construction, the time and foreign_time should be greater
    // than the times at the sources of the tracks (and oppositely for the targets).
    // Some numerical errors can sneak it and, if so, snap the times.
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

    if(foreign_time_at_collision > time_at_fmc_track_destination)
    {
      CGAL_precondition(foreign_time_at_collision - tolerance_ <= time_at_fmc_track_destination);
      foreign_time_at_collision = time_at_fmc_track_destination;
    }
#endif

    CGAL_postcondition(time_at_collision >= mc.current_time());
    CGAL_postcondition(time_at_collision <= mc.time_at_closest_target());
    CGAL_postcondition(foreign_time_at_collision >= time_at_fmc_track_source);
    CGAL_postcondition(foreign_time_at_collision <= time_at_fmc_track_destination);

    if(tc.is_collision_earlier_than_current_best(time_at_collision, foreign_time_at_collision, is_fmc_moving_on_track))
    {
      // both values are used later when we snap times/points
      const FT time_at_closest_collision_memory = tc.time_at_closest_collision;
      const FT foreign_time_at_closest_collision_memory = tc.foreign_time_at_closest_collision;

      tc.reset();
      tc.time_at_closest_collision = time_at_collision;
      tc.fmc_id = fmc.id();
      tc.is_foreign_motorcycle_moving_on_track = is_fmc_moving_on_track;
      tc.foreign_time_at_closest_collision = foreign_time_at_collision;

      // Step 3:
      // Although there does not exist a point at the location of the collision,
      // this point might be at the same time from the source of the track
      // as another point due to numerical errors.
      std::pair<TPC_iterator, bool> res = mc.has_target_at_time(time_at_collision);
      if(res.second)
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "Motorcycle #" << mc.id() << " already has a target at time: " << time_at_collision << std::endl;
#endif

        TPC_iterator target_point = res.first;
        CGAL_assertion(target_point->second == time_at_collision);
        Node_ptr alternate_collision = target_point->first;

        // If the times are equal, the points should be very close
        CGAL_assertion(CGAL::squared_distance(alternate_collision->point(), collision_point) < tolerance_);

        // Temporal snap: the collision is now that existing point instead
        tc.is_closest_collision_already_in_dictionary = true;
        tc.closest_collision = alternate_collision;

        return COLLISION;
      }

      std::pair<TPC_iterator, bool> fmc_res = fmc.has_target_at_time(foreign_time_at_collision);
      if(fmc_res.second) // there is already a target at that time
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "Foreign motorcycle #" << fmc.id() << " already has a target at time: " << foreign_time_at_collision << std::endl;
#endif

        TPC_iterator target_point = fmc_res.first;
        Node_ptr alternate_foreign_collision = target_point->first;
        CGAL_assertion(alternate_foreign_collision->face() == fmc.current_face());
        CGAL_assertion(target_point->second == foreign_time_at_collision);

        tc.is_closest_collision_already_in_dictionary = true;
        tc.closest_collision = alternate_foreign_collision;

        return COLLISION;
      }

      // At this point, we have a new location at an unknown time...
#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
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

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "successful snapping to " << std::endl << *(is_snappable.first) << std::endl;
        std::cout << "new times: " << visiting_time << " " << time_at_closest_collision_memory << std::endl;
#endif

        // We have snapped so we are igoring times that we had formely set up as best, but
        // we still need to ensure that those times are better than the previous one.
        CGAL_assertion(visiting_time <= time_at_closest_collision_memory);
        CGAL_assertion(visiting_time < time_at_closest_collision_memory ||
                       foreign_visiting_time < foreign_time_at_closest_collision_memory);

        tc.is_closest_collision_already_in_dictionary = true;
        tc.closest_collision = is_snappable.first;
        tc.time_at_closest_collision = visiting_time;
        tc.foreign_time_at_closest_collision = foreign_visiting_time;

          return SNAPPED_COLLISION_TO_EXISTING_POINT;
      }
#endif

      // Couldn't snap to anything, 'collision_location' is definitely a new point
      tc.is_closest_collision_already_in_dictionary = false;
      tc.closest_collision_location = collision_location;

      return COLLISION;
    }
  }

  return NO_COLLISION;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_complete_track(const Motorcycle& mc, const Segment_2& mcs,
                                   const Track_segment& fmc_track,
                                   // below are out parameters
                                   Collision_information& tc)
{
  const Motorcycle& fmc = motorcycle(fmc_track.motorcycle_id());

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "-----------------------------" << std::endl;
  std::cout << "¤ Checking for intersection with the complete track of motorcycle #" << fmc.id() << std::endl;
#endif

  // 'false' because the motorcycle is not moving on that track // @todo replace with named enums
  return find_collision_between_tracks(mc, mcs, fmc, fmc_track, false /*not moving*/, tc);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Collision_return
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_collision_with_live_motorcycle(Motorcycle& mc, const Segment_2& mcs,
                                    const Motorcycle& fmc,
                                    // below are out parameters
                                    Collision_information& tc)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE_PLUS
  std::cout << "-----------------------------" << std::endl;
  std::cout << "¤ Checking for intersection with live motorcycle #" << fmc.id() << std::endl;
#endif

  if(// the motorcycles must be different
     mc.id() == fmc.id() ||
     // the motorcycles must be in the same face
     mc.current_face() != fmc.current_face() ||
     // the foreign motorcycle must be in motion
     fmc.is_crashed())
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE_PLUS
    std::cout << " ignoring fmc..." << std::endl;
    std::cout << "  > motorcycles #" << mc.id() << " and #" << fmc.id() << std::endl;
    std::cout << "  > faces: " << mc.current_face() << " and " << fmc.current_face() << std::endl;
    std::cout << "  > crashed status: " << fmc.is_crashed() << std::endl;
#endif
    return NO_COLLISION;
  }

  Track_segment fmc_track(fmc.id(), fmc.current_position(), fmc.current_time(),
                          fmc.closest_target(), fmc.time_at_closest_target());

  // 'true' because fmc is currently moving on that track
  return find_collision_between_tracks(mc, mcs, fmc, fmc_track, true, tc);
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

  CGAL_precondition(!mc.is_crashed());
  CGAL_precondition(!mc.targets().empty());

  // The motorcycles must be on the same face
  CGAL_precondition(mc.current_face() == mc.closest_target()->face());

  // Use the barycentric coordinate systems to compute intersections
  const Point_2 s = geom_traits().construct_point_2_object()(mc.current_location().second[0],
                                                             mc.current_location().second[1]);
  const Point_2 t = geom_traits().construct_point_2_object()(mc.closest_target()->location().second[0],
                                                             mc.closest_target()->location().second[1]);
  const Segment_2 mc_tentative_track = geom_traits().construct_segment_2_object()(s, t);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "MC tentative track:" << std::endl
            << "source: " << &*(mc.current_position()) << " " << *(mc.current_position()) << std::endl
            << "target: " << &*(mc.closest_target()) << " " << *(mc.closest_target()) << std::endl;
#endif

  // A degenerate tentative track has no interesting collisions
  if(mc_tentative_track.is_degenerate())
    return NO_COLLISION;

  // Checking for intersection is done in the following steps:
  // - 1: Check with complete tracks in the face
  // - 2: Check the motorcycles that are currently moving in the face ("live")
  // - 3: Check for intersections at the border of the face with tracks in adjacent faces

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_ COMPLETE TRACKS _-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
  std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << std::endl;
#endif

  Collision_return res = NO_COLLISION;

  // Step 1: check complete tracks
  const face_descriptor mc_fd = mc.current_face();
  TFM_iterator it = track_face_map_.find(mc_fd);
  if(it != track_face_map_.end())
  {
    const Track_segment_ptr_container& face_tracks = it->second;

    typename Track_segment_ptr_container::const_iterator tl_it = face_tracks.begin();
    typename Track_segment_ptr_container::const_iterator tl_end = face_tracks.end();
    for(; tl_it!=tl_end; ++tl_it)
    {
      const Track_segment& ts = *(*tl_it);

      Collision_return r = find_collision_with_complete_track(mc, mc_tentative_track, ts, tc);

      // Need to keep the foreign track in memory to add a new point on the confirmed track...
      if(r == COLLISION)
      {
        CGAL_assertion(!tc.is_foreign_motorcycle_moving_on_track);
        tc.foreign_track = *tl_it;
      }

      res = res | r;

      if(res == SNAPPED_COLLISION_TO_EXISTING_POINT)
        return res;
    }
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
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
    res = res | find_collision_with_live_motorcycle(mc, mc_tentative_track, fmc, tc);

    if(res == SNAPPED_COLLISION_TO_EXISTING_POINT)
      return res;
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
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

    // this part is brute force for now, but can be done in O(nlogn) by sorting
    // according to the slopes (farthest intersections happen when the slopes
    // of the motorcycles are close)
    for(std::size_t fmc_id = 0; fmc_id<nm; ++fmc_id)
    {
      // segment - segment, segment - ray, or ray-ray intersections @todo
    }
  }

  // Slightly increase the size of the bbox to strictly contain all points

  // Manually create the mesh with Euler operations
#endif
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
has_motorcycle_reached_crashing_point(const Motorcycle& mc) const
{
  return (// multiple motorcycles will reach mc's current position at the same time
          mc.has_reached_simultaneous_collision_point() ||
          // the current position might be blocked (possibly in a representation in another face)
          is_motorcycle_position_blocked(mc));
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
has_motorcycle_reached_final_destination(const Motorcycle& mc) const
{
  return mc.is_destination_final();
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
initialize_motorcycles()
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Initialize motorcycles" << std::endl;
#endif
  namespace PMP = CGAL::Polygon_mesh_processing;

  AABB_tree tree;
  AABB_tree_VPM vpm(mesh());

  // if no mesh has been given in input, generate a mesh made of a single quad face
  // that contains all the interesting motorcycle interactions (i.e. crashes)
  if(!is_mesh_provided)
    generate_enclosing_face();

  if(is_AABB_tree_needed())
    PMP::build_AABB_tree(mesh(), tree, parameters::vertex_point_map(vpm));

  MCC_it mc_it = motorcycles().begin(), mc_end = motorcycles().end();
  for(; mc_it!=mc_end; ++mc_it)
  {
    Motorcycle& mc = motorcycle(mc_it);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "     _" << std::endl;
    std::cout << "   D/_" << std::endl;
    std::cout << "   /(__`=-/" << std::endl;
    std::cout << "  (o)     (o)" << std::endl;
    std::cout << "Initializing motorcycle #" << mc.id() << std::endl;
#endif

    // Add the origin to the node dictionary
    const Point_or_location& input_origin = mc.input_origin();
    add_origin_node(mc, input_origin, tree, vpm);

    // Compute the destination if needed, and add it to the node dictionary
    Optional_point_or_location& input_destination = mc.input_destination();
    add_destination_node(mc, input_destination);

    // Initialize the motorcycle targets
    mc.targets().insert(std::make_pair(mc.origin(), mc.time_at_origin()));

    if(mc.origin() != mc.destination())
      mc.targets().insert(std::make_pair(mc.destination(), mc.time_at_destination()));

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
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
initialize_tracing()
{
  initialize_motorcycles();
  motorcycle_pq_.initialize(motorcycles_);

  // Two motorcycles leaving from the same point is allowed, as long as their
  // directions are different. At any later stage, a point being reached
  // by two motorcycles will mean a crash.
  crash_motorcycles_with_same_origins_and_directions();

#ifdef CGAL_MOTORCYCLE_GRAPH_OUTPUT
  output_motorcycles_origins_and_destinations();
#endif
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
is_AABB_tree_needed() const
{
  // an AABB tree must be built if at least one origin is given as a geometric point
  MCC_cit mc_it = motorcycles().begin(), mc_end = motorcycles().end();
  for(; mc_it!=mc_end; ++mc_it)
  {
    const Motorcycle& mc = motorcycle(mc_it);
    if(const Point* p = boost::get<Point>(&(mc.input_origin()))) // input was given as a 'Point'
      return true;
  }

  return false;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
is_motorcycle_position_blocked(const Motorcycle& mc) const
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  if(mc.has_reached_blocked_point())
    return true;

  // @todo check the correctness of below
  // to avoid self blocking while crossing mesh edges
  Node_ptr position = mc.current_position();
  if(position->earliest_motorcycle()->first < mc.current_time())
    return true;

  return false;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Face_location
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
locate(const Point& p, const AABB_tree& tree, const AABB_tree_VPM vpm) const
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  // An AABB tree is a 3D structure, so we need to convert the point to a Point_3.
  // If the point is already a Point_3, this doesn't do anything.
  Point_to_Point_3 to_p3;
  const typename Point_to_Point_3::Point_3& p3 = to_p3(p);

  CGAL_static_assertion((boost::is_same<typename Point_to_Point_3::Point_3,
                                        typename AABB_tree::AABB_traits::Point_3>::value));

  Face_location loc = PMP::locate_with_AABB_tree(p3, tree, mesh(), parameters::vertex_point_map(vpm));

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
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "---" << std::endl;
    std::cout << "Driving priority queue:" << std::endl << motorcycle_pq_ << std::endl;
#endif

    // get the earliest available event
    Motorcycle_PQE pqe = motorcycle_pq_.top();
    Motorcycle& mc = pqe.motorcycle();

    std::cout << "Driving priority queue size: " << motorcycle_pq_.size();
    std::cout << " (closest time: " << mc.time_at_closest_target() << ")" << std::endl << std::endl;

    // move the motorcycle to the closest target, which becomes its confirmed position
    drive_to_closest_target(mc);

    if(mc.current_position() == mc.destination())
    {
      if(has_motorcycle_reached_final_destination(mc) ||
         has_motorcycle_reached_crashing_point(mc))
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "Reached motorcycle's crashing point:" << std::endl
                  << " - final destination: " << has_motorcycle_reached_final_destination(mc) << std::endl
                  << " - blocked: " << is_motorcycle_position_blocked(mc) << std::endl
                  << " - simultaneous collision: " << mc.has_reached_simultaneous_collision_point() << std::endl;
#endif
        crash_motorcycle(mc);
      }
      else // Not crashing yet, try to compute the next path
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cout << "Reached destination: " << mc.destination()->point();
        std::cout << " Now computing motorcycle's next path..." << std::endl;
#endif
        // Clear any unnecessary targets that might have been built
        mc.clear_targets();

        if(initialize_next_path(mc))
        {
          // A new path was found and set up, update the queue and continue
          motorcycle_pq_.update(mc);

          // Note that we have intentionally not blocked the point in this case!!
          continue;
        }
        else
        {
          // Couldn't find a next destination, crash the motorcycle
          crash_motorcycle(mc);
        }
      }
    }
    // The motorcycle has not reached its destination, but still might be crashing
    else if(has_motorcycle_reached_crashing_point(mc) &&
            // hackish to prevent multiple motorcycles starting from the same origin
            // (but with different directions) from blocking each other
            //
            // @fixme: can specify starting time ==> '0' is not acceptable
            // number of elements in the track == 0 ?
            mc.current_time() != 0.)
    {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "Reached crashing point:" << std::endl
                << " - blocked: " << is_motorcycle_position_blocked(mc) << std::endl
                << " - simultaneous collision: " << mc.has_reached_simultaneous_collision_point() << std::endl;
#endif
      crash_motorcycle(mc);
    }
    // the motorcycle can continue without issue towards its destination
    else
    {
      // A bunch of output parameters are regrouped into the 'Collision_information' struct,
      // which describes the best (closest to mc.current_position()) tentative collision.
      Collision_information tc(mc.time_at_closest_target() /*maximum allowed time*/);

      // Check for potential collisions within the face for the next move of 'mc'
      Collision_return res = find_collision(mc, tc);

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE // @todo use another macro for snapping stuff
      if(res == SNAPPED_COLLISION_TO_EXISTING_POINT)
      {
        // This handles the configuration where the collision computed is close
        // to an existing point. We only add combinatorial information using this
        // existing point and re-add the current positions as targets.
        // On the next pass, combinatorial checks will select the existing point
        // as valid intersection directly.
        visit_point(mc, tc);

        // - Re-add the current positions of the motorcycles to re-evaluate
        //   potential intersections in the tentative path
        // - Update the priority queue for the two motorcycles
        if(tc.closest_collision != mc.current_position())
          mc.add_target(mc.current_position(), mc.current_time());
        motorcycle_pq_.update(mc);

        Motorcycle& fmc = motorcycle(tc.fmc_id);
        if(!fmc.is_crashed())
        {
          // Check that it is a target to know if it makes sense to update 'fmc'
          if(fmc.has_target(tc.closest_collision).second && tc.closest_collision != fmc.current_position())
            fmc.add_target(fmc.current_position(), fmc.current_time());

          motorcycle_pq_.update(fmc);
        }

        continue;
      }
      else
#endif
      if(res == COLLISION)
      {
        CGAL_assertion(tc.found_collision());
        treat_collision(mc, tc);
      }
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      else
      {
        CGAL_assertion(res == NO_COLLISION);
        std::cout << " No collision was found!" << std::endl;
      }
#endif

      // The target list of 'mc' was modified and the PQ must be updated.
      // The PQ entry of 'fmc' is modified in 'treat_collision()', if needed.
      CGAL_assertion(!mc.is_crashed());
      motorcycle_pq_.update(mc);

      // Block the current position
      mc.current_position()->block();
    }
  }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Finished tracing" << std::endl;
#endif
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
treat_collision(Motorcycle& mc, const Collision_information& collision_info)
{
  const std::size_t foreign_motorcycle_id = collision_info.fmc_id;
  Motorcycle& fmc = motorcycle(foreign_motorcycle_id);
  const FT time_at_collision = collision_info.time_at_closest_collision;
  const FT time_at_foreign_collision = collision_info.foreign_time_at_closest_collision;

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << std::endl << "+++++++++++ Treat collision [PREPROCESSING] +++++++++++" << std::endl;
  std::cout << " - motorcycle:" << std::endl << mc << std::endl
            << " - foreign_motorcycle:" << std::endl << fmc << std::endl
            << " - moving: " << collision_info.is_foreign_motorcycle_moving_on_track << std::endl
            << " - time_at_collision: " << time_at_collision << std::endl
            << " - time_at_foreign_collision: " << time_at_foreign_collision << std::endl;
#endif

  const face_descriptor fd = mc.current_face();
  const face_descriptor ffd =
    collision_info.is_foreign_motorcycle_in_different_face ? collision_info.foreign_motorcycle_face
                                                           : fd;

  CGAL_assertion(ffd != boost::graph_traits<Triangle_mesh>::null_face());
  CGAL_assertion((collision_info.is_foreign_motorcycle_in_different_face && fd != ffd) || fd == ffd);


  // Insert the collision point in the dictionary, if needed.
  Node_ptr collision;
  if(!collision_info.is_closest_collision_already_in_dictionary)
  {
    // Motorcycle info will be added later.
    std::pair<Node_ptr, bool> entry = nodes().insert(collision_info.closest_collision_location, mesh());
    collision = entry.first;

    if(!entry.second)
    {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cerr << "Warning: collision location actually already existed in the dictionary:"
                << std::endl << *(entry.first) << std::endl;
#endif
    }
  }
  else
  {
    collision = collision_info.closest_collision;
  }

  Node_ptr collision_in_fd = nodes().get_sibling(collision, fd);
  Node_ptr collision_in_ffd = nodes().get_sibling(collision, ffd);

  // Some sanity tests
  CGAL_postcondition(collision_in_fd->face() == fd);
  CGAL_postcondition(collision_in_ffd->face() == ffd);
  CGAL_postcondition_code(if(fd != ffd) {)
  CGAL_postcondition(collision_in_ffd->is_sibling(collision_in_fd->location()));
  CGAL_postcondition(collision_in_fd->is_sibling(collision_in_ffd->location())); // just to be extra sure
  CGAL_postcondition_code(})

  Track_segment_ptr foreign_track = collision_info.foreign_track;

  if(!collision_info.is_foreign_motorcycle_moving_on_track)
  {
    CGAL_assertion(foreign_track != Track_segment_ptr());
    CGAL_assertion(foreign_track->motorcycle_id() == foreign_motorcycle_id);
    CGAL_assertion(foreign_track->time_at_source() <= time_at_foreign_collision);
    CGAL_assertion(foreign_track->time_at_target() >= time_at_foreign_collision);
  }
  else
  {
    // @todo is that even needed  --> just check vs default constructed iterator?
    foreign_track = fmc.track().end();
  }

  // Treat_collision handles all types of collisions
  treat_collision(mc, collision_in_fd, time_at_collision,
                  fmc, collision_in_ffd, time_at_foreign_collision,
                  foreign_track);
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
treat_collision(Motorcycle& mc, Node_ptr collision_point, const FT time_at_collision,
                Motorcycle& fmc, Node_ptr foreign_collision_point, const FT foreign_time_at_collision,
                Track_segment_ptr foreign_track)
{
  // @todo factorize this a bit so we don't have multiple calls to "has_motorcycle"
  // followed by "add_motorcycle": this is all in log(n) complexity (admittedly,
  // log(n) on something that is unlikely to contain more than 2 elements, but still)

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "+++++++++++ Treat collision +++++++++++" << std::endl;
  std::cout << " - collision_point: " << &*collision_point << std::endl << *(collision_point) << std::endl;
  if(collision_point != foreign_collision_point)
    std::cout << " - foreign collision_point: " << &*foreign_collision_point << std::endl << *(foreign_collision_point) << std::endl;
#endif

  // Some sanity checks
  CGAL_precondition(mc.id() != std::size_t(-1));
  CGAL_precondition(fmc.id() != std::size_t(-1));

  CGAL_precondition(collision_point != Node_ptr());
  CGAL_precondition(foreign_collision_point != Node_ptr());
  CGAL_precondition(collision_point->point() == foreign_collision_point->point());
  CGAL_precondition(collision_point->face() == mc.current_face());

  // Can't give an upper bound on the (foreign_)time_at_collision due to front collisions
  CGAL_precondition(time_at_collision >= mc.current_time());

  if(// the impact is closer than the next target
     time_at_collision <= mc.time_at_closest_target() &&
     // the collision is not the next target of 'mc' or the foreign track
     // does not know the collision point yet
     (collision_point != mc.closest_target() ||
      !collision_point->has_motorcycle(fmc.id(), foreign_time_at_collision)))
  {
    if(!collision_point->has_motorcycle(mc.id(), time_at_collision))
    {
      // Call the halving structure to create a new point
      std::pair<Node_ptr, FT> halving_entity =
        compute_halving_point(mc, mc.current_position(), mc.current_time(),
                              collision_point, time_at_collision);
      Node_ptr halving_point = halving_entity.first;
      const FT time_at_halving_point = halving_entity.second;

      // Degeneracies should have been caught before
      CGAL_postcondition(halving_point != mc.current_position() &&
                         halving_point != collision_point);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "Adding collision point: " << &*collision_point
                << " and halving point: " << &*halving_point
                << " to motorcycle #" << mc.id() << std::endl;
#endif

      mc.add_target(collision_point, time_at_collision);
      mc.add_target(halving_point, time_at_halving_point);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "Adding motorcycle #" << mc.id()
                << " to collision point: " << &*collision_point
                << " and halving point: " << &*halving_point << std::endl;
#endif

      halving_point->add_motorcycle(mc.id(), time_at_halving_point);
      collision_point->add_motorcycle(mc.id(), time_at_collision);

      CGAL_postcondition(mc.has_target_at_time(collision_point, time_at_collision));
    }
    // If we have snapped the collision point to the current position, re-add it to the targets.
    // Note that we won't find the same intersection again because 'collision_point'
    // (which is 'mc.current_position') now combinatorially knows that there is
    // an intersection with the foreign motorcycle at 'mc.current_time' (and it will be ignored).
    // See "Check #1: known collision at current_position"
    else if(collision_point == mc.current_position())
    {
      mc.add_target(collision_point, time_at_collision);
    }

    // Now, do the same for the foreign motorcycle
    if(!foreign_collision_point->has_motorcycle(fmc.id(), foreign_time_at_collision) &&
       // If the motorcycle is crashed, ignore points after the final destination
       ((fmc.is_crashed() && foreign_time_at_collision < fmc.current_time()) ||
       // Otherwise, ignore points that are farther than the current closest point for the foreign track
       // (otherwise you can get nasty stuff like halving points == existing points, etc.)
        foreign_time_at_collision < fmc.time_at_closest_target()))
    {
      // It is useful to know that the collision point is on the foreign track,
      // even if the collision point is on the confirmed part of the track.
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
      std::cout << "Adding foreign motorcycle #" << fmc.id()
                << " to foreign collision point: " << &*foreign_collision_point << std::endl;
#endif
      foreign_collision_point->add_motorcycle(fmc.id(), foreign_time_at_collision);

      if(// the collision point is not on the confirmed track for the foreign mc
         foreign_time_at_collision > fmc.current_time())
      {
        if(!fmc.is_crashed()) // do not add targets to a crashed motorcycle
        {
          // Call the halving structure to create a new point
          std::pair<Node_ptr, FT> foreign_halving_entity =
            compute_halving_point(fmc, fmc.current_position(), fmc.current_time(),
                                  foreign_collision_point, foreign_time_at_collision);
          Node_ptr foreign_halving_point = foreign_halving_entity.first;
          const FT foreign_time_at_halving_point = foreign_halving_entity.second;

          // Degeneracies should have been caught before
          CGAL_postcondition(foreign_halving_point != fmc.current_position() &&
                             foreign_halving_point != foreign_collision_point);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
          std::cout << "Adding foreign collision point: " << &*foreign_collision_point
                    << " and halving point: " << &*foreign_halving_point
                    << " to motorcycle #" << fmc.id() << std::endl;
#endif
          fmc.add_target(foreign_collision_point, foreign_time_at_collision);
          fmc.add_target(foreign_halving_point, foreign_time_at_halving_point);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
          std::cout << "Adding foreign motorcycle #" << fmc.id()
                    << " to halving point: " << &*foreign_halving_point << std::endl;
#endif
          foreign_halving_point->add_motorcycle(fmc.id(), foreign_time_at_halving_point);
          CGAL_postcondition(fmc.has_target_at_time(foreign_collision_point, foreign_time_at_collision));

          // The target list of the foreign motorcycle was modified and the queue must be updated
          motorcycle_pq_.update(fmc);
        }
      }
      else // foreign time <= current_time --> collision point is on the confirmed foreign track
      {
        foreign_collision_point->block();

        // Add it to the track of the foreign motorcycle
        CGAL_assertion(foreign_track != fmc.track().end());
        Track_segment_ptr ts_ptr = fmc.track().split_track_segment(foreign_track,
                                                                   foreign_collision_point,
                                                                   foreign_time_at_collision);

        // add the new track segment to the face map
        add_track_segment_to_map(foreign_track->face(), ts_ptr);
      }
    }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << std::endl << "[[ Post-treatment... ]]" << std::endl;
    std::cout << "Motorcycles involved: " << std::endl << mc << std::endl << fmc << std::endl;

    std::cout << "collision point:" << std::endl << *collision_point << std::endl;
    if(collision_point != foreign_collision_point)
      std::cout << "foreign collision point:" << std::endl << *foreign_collision_point << std::endl;
#endif
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
std::pair<typename Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::Node_ptr, bool>
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
find_close_existing_point(const Face_location& location, const Point& p,
                          const bool allow_same_location) const
{
  // The collision must not already be in the dictionary
  CGAL_expensive_precondition_code(!nodes().find(location).second);

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
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
try_to_snap_to_close_existing_point(Node_ptr& node)
{
  // 'false' because we want to ignore the same location
  std::pair<Node_ptr, bool> is_snappable = find_close_existing_point(node->location(), node->point(), false);
  if(is_snappable.second)
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Snapping: " << std::endl << *node << std::endl
              << " to existing point: " << std::endl << *(is_snappable.first) << std::endl;
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
      // of the two close nodes.
      CGAL_assertion(false);
    }

    node = is_snappable.first;
    return true;
  }

  return false;
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
visit_point(Motorcycle& mc, Collision_information& tc)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " ---- Visiting point " << &*(tc.closest_collision) << std::endl
                                       << *(tc.closest_collision) << std::endl;
  std::cout << "with mc: " << std::endl << mc << std::endl;
  std::cout << "with fmc: " << std::endl << motorcycle(tc.fmc_id) << std::endl;
  std::cout << "times: " << tc.time_at_closest_collision << " || " << tc.foreign_time_at_closest_collision << std::endl;
#endif

  const bool are_times_equal =
      (CGAL::abs(tc.foreign_time_at_closest_collision - tc.time_at_closest_collision) <= tolerance_);

  Motorcycle& fmc = motorcycle(tc.fmc_id);

  // First, check if 'mc' is known at the collision point
  FT min_visiting_time = tc.time_at_closest_collision - tolerance_;
  FT max_visiting_time = tc.time_at_closest_collision + tolerance_;
  FT visiting_time;
  bool is_visited_by_mc = tc.closest_collision->has_motorcycle(
                            mc.id(), min_visiting_time, max_visiting_time, visiting_time);

  FT min_foreign_visiting_time = tc.foreign_time_at_closest_collision - tolerance_;
  FT max_foreign_visiting_time = tc.foreign_time_at_closest_collision + tolerance_;
  FT foreign_visiting_time;
  bool is_visited_by_fmc = tc.closest_collision->has_motorcycle(
                             tc.fmc_id, min_foreign_visiting_time,
                             max_foreign_visiting_time, foreign_visiting_time);

  // Make sure that even if we snap, the times stay equal
  if(is_visited_by_fmc)
  {
    if(is_visited_by_mc)
    {
      if(are_times_equal) // visited by both, equal times
      {
        // Can't change times otherwise we create inconsistencies
        CGAL_assertion(CGAL::abs(tc.foreign_time_at_closest_collision - tc.time_at_closest_collision) <= tolerance_);
      }
    }
    else if(are_times_equal)// only visited by fmc, equal times
      tc.time_at_closest_collision = foreign_visiting_time;
  }
  else if(is_visited_by_mc && are_times_equal) // only visited by mc, equal times
    tc.foreign_time_at_closest_collision = visiting_time;

  // Need to get the representants in the correct faces
  face_descriptor mc_fd = mc.current_face();
  face_descriptor fmc_fd = tc.is_foreign_motorcycle_in_different_face ? tc.foreign_motorcycle_face
                                                                      : mc_fd;

  Node_ptr collision_in_mc_face = nodes().get_sibling(tc.closest_collision, mc_fd);
  Node_ptr collision_in_fmc_face = nodes().get_sibling(tc.closest_collision, fmc_fd);

  if(is_visited_by_mc)
  {
    // Consistency: there can't be a target that is ever so slightly before the current point
    CGAL_assertion(visiting_time >= mc.current_time());

    CGAL_assertion_code(if(visiting_time == mc.current_time()))
    CGAL_assertion(collision_in_mc_face == mc.current_position());
    CGAL_assertion_code(if(visiting_time == mc.time_at_closest_target()))
    CGAL_assertion(collision_in_mc_face == mc.closest_target());

    if(!mc.has_target(collision_in_mc_face).second) // nothing to do if it's already a target
      mc.add_target(collision_in_mc_face, visiting_time);
  }
  else // the snapping collision point is not yet visited by 'mc'
  {
    // Consistency: there can't be a target around that time: snapping should create safe zones around points
    CGAL_assertion(!mc.has_target_at_time(min_visiting_time, max_visiting_time).second);
    CGAL_assertion(tc.time_at_closest_collision > mc.current_time());
    CGAL_assertion(tc.time_at_closest_collision < mc.time_at_closest_target());

    mc.add_target(collision_in_mc_face, tc.time_at_closest_collision);
    collision_in_mc_face->add_motorcycle(mc.id(), tc.time_at_closest_collision);
  }

  // Same, for the foreign motorcycle 'fmc'
  if(is_visited_by_fmc)
  {
    if(// can't add targets if the motorcycle is crashed
       !fmc.is_crashed() &&
       // nothing to do if it's already a target
       !fmc.has_target(collision_in_fmc_face).second &&
       // can't be a target if it's younger than the current time
       foreign_visiting_time >= fmc.current_time())
    {
      fmc.add_target(collision_in_fmc_face, foreign_visiting_time);
    }
  }
  else // the snapping collision point is not yet visited by 'fmc'
  {
    CGAL_assertion(!fmc.has_target_at_time(min_foreign_visiting_time, max_foreign_visiting_time).second);
    CGAL_assertion(tc.foreign_time_at_closest_collision <= fmc.time_at_destination());

    if(!fmc.is_crashed() && tc.foreign_time_at_closest_collision >= fmc.current_time())
      fmc.add_target(collision_in_fmc_face, tc.foreign_time_at_closest_collision);

    collision_in_fmc_face->add_motorcycle(fmc.id(), tc.foreign_time_at_closest_collision);

    // Check if the point is on the confirmed part of the foreign motorcycle's track
    if(tc.foreign_time_at_closest_collision <= fmc.current_time())
      collision_in_fmc_face->block();
  }

  CGAL_postcondition(mc.has_target(collision_in_mc_face).second);
  CGAL_postcondition(fmc.is_crashed() ||
                     tc.foreign_time_at_closest_collision < fmc.current_time() ||
                     fmc.has_target(collision_in_fmc_face).second);
  CGAL_postcondition(tc.closest_collision->has_motorcycle(mc.id()));
  CGAL_postcondition(tc.closest_collision->has_motorcycle(fmc.id()));

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "Post-visit: " << std::endl;
  std::cout << "collision: " << std::endl << *(tc.closest_collision) << std::endl;
  std::cout << "mc: " << std::endl << mc << std::endl;
  std::cout << "fmc: " << std::endl << fmc << std::endl;
#endif
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
bool
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
is_valid() const
{
  std::cout << "Checking trace validity..." << std::endl;

  // brute force track validity check @todo make it not as brute force
  MCC_cit mc_it = motorcycles().begin(), mc_end = motorcycles().end();
  for(; mc_it!=mc_end; ++mc_it)
  {
    const Motorcycle& mc = motorcycle(mc_it);
    const Track& mc_track = mc.track();
    std::size_t mc_track_size = mc_track.size();

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "Track of motorcycle: " << mc.id() << std::endl << mc_track << std::endl;
#endif

    if(mc_track_size == 0)
      continue;

    const Track_segment& first_ts = mc_track.front();
    Node_ptr previous_segment_target = first_ts.target();
    FT time_at_previous_segment_target = first_ts.time_at_target();

    typename Track::const_iterator mc_tcit = mc_track.begin(), end = mc_track.end();
    for(; mc_tcit!=end; ++mc_tcit)
    {
      const Track_segment& mc_ts = *mc_tcit;

      const Node_ptr mc_track_source = mc_ts.source();
      const FT time_at_mc_track_source = mc_ts.time_at_source();
      const Node_ptr mc_track_target = mc_ts.target();
      const FT time_at_mc_track_target = mc_ts.time_at_target();

      if(mc_track_source->face() != mc_track_target->face())
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cerr << "Error: track is crossing faces (";
        std::cerr << mc_track_source->point() << " || " << mc_track_target->point() << ")" << std::endl;
        std::cerr << "motorcycle id: " << mc.id() << std::endl;
#endif
        return false;
      }

      // check the visiting times consistency
      if(time_at_mc_track_source > time_at_mc_track_target)
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cerr << "Error: source must be before target (" << time_at_mc_track_source
                  << " vs " << time_at_mc_track_target << ")" << std::endl;
        std::cerr << "motorcycle id: " << mc.id() << std::endl;
#endif
        return false;
      }

      // Ignore degenerate segments, except if the whole track is a single degenerate segment
      const bool is_mc_degenerate = (mc_track_source == mc_track_target);
      if(is_mc_degenerate && mc_track_size > 1)
        continue;

      // check the continuity between two consecutive tracks
      if(mc_track_source->point() != previous_segment_target->point() ||
         time_at_mc_track_source != time_at_previous_segment_target)
      {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
        std::cerr << "Error in track continuity: points ("
                  << mc_track_source->point() << ") and (" << previous_segment_target->point() << ") "
                  << "times " << time_at_mc_track_source << " & " << time_at_previous_segment_target << std::endl;
        std::cerr << "motorcycle id: " << mc.id() << std::endl;
#endif
        return false;
      }

      previous_segment_target = mc_track_target;
      time_at_previous_segment_target = time_at_mc_track_target;

      // Here comes the mega brute force: check each track segment against ALL the other track segments
      const Point_2 ts = geom_traits().construct_point_2_object()(mc_track_source->barycentric_coordinate(0),
                                                                  mc_track_source->barycentric_coordinate(1));
      const Point_2 tt = geom_traits().construct_point_2_object()(mc_track_target->barycentric_coordinate(0),
                                                                  mc_track_target->barycentric_coordinate(1));
      Segment_2 s = geom_traits().construct_segment_2_object()(ts, tt);

      MCC_cit fmc_it = motorcycles().begin(), fmc_end = motorcycles().end();
      for(; fmc_it!=fmc_end; ++fmc_it)
      {
        const Motorcycle& fmc = motorcycle(fmc_it);
        const Track& fmc_track = fmc.track();
        std::size_t fmc_track_size = fmc_track.size();

        if(fmc_track_size == 0)
          continue;

        typename Track::const_iterator fmc_tcit = fmc_track.begin(), fend = fmc_track.end();
        for(; fmc_tcit!=fend; ++fmc_tcit)
        {
          const Track_segment& fmc_ts = *fmc_tcit;

          const Node_ptr fmc_track_source = fmc_ts.source();
          const FT time_at_fmc_track_source = fmc_ts.time_at_source();
          const Node_ptr fmc_track_target = fmc_ts.target();
          const FT time_at_fmc_track_target = fmc_ts.time_at_target();

          if(mc_ts == fmc_ts)
            continue;

          // @fixme there can be valid intersections on the edge
          if(mc_track_source->face() != fmc_track_source->face())
            continue;

          // Ignore degenerate segments, except if the whole track is a single degenerate segment
          const bool is_fmc_degenerate = (fmc_track_source == fmc_track_target);
          if(is_fmc_degenerate && fmc_track_size > 1)
            continue;

          const Point_2 fts = geom_traits().construct_point_2_object()(fmc_track_source->barycentric_coordinate(0),
                                                                       fmc_track_source->barycentric_coordinate(1));
          const Point_2 ftt = geom_traits().construct_point_2_object()(fmc_track_target->barycentric_coordinate(0),
                                                                       fmc_track_target->barycentric_coordinate(1));
          Segment_2 fs = geom_traits().construct_segment_2_object()(fts, ftt);

          if(geom_traits().do_intersect_2_object()(s, fs))
          {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
            std::cout << "Intersection ¤~~~~~~~~~~~~~~~~~¤" << std::endl;
            std::cout << "motorcycle #" << mc.id() << " (track size: " << mc_track.size();
            std::cout << ") with motorcycle #" << fmc.id() << " (track size: " << fmc_track.size() << ")" << std::endl;
            std::cout << "DECITs:" << std::endl << *mc_track_source << std::endl << *mc_track_target << std::endl
                                                << *fmc_track_source << std::endl << *fmc_track_target << std::endl;
#endif

            // Check that the only possible intersection is an extremity
            if(is_mc_degenerate)
            {
              if(is_fmc_degenerate) // both tracks are generate
              {
                CGAL_assertion(mc_track_source == fmc_track_source);
              }
              else // only mc is denegerate
              {
                CGAL_assertion((fmc_track_source == mc_track_source && fmc_track_target != mc_track_source) ||
                               (fmc_track_source != mc_track_source && fmc_track_target == mc_track_source));
              }
            }
            else if(is_fmc_degenerate) // degenerate fmc track
            {
              CGAL_assertion((mc_track_source == fmc_track_source && mc_track_target != fmc_track_source) ||
                             (mc_track_source != fmc_track_source && mc_track_target == fmc_track_source));
            }
            else // nothing is degenerate
            {
              CGAL_assertion((mc_track_source == fmc_track_source && mc_track_source != fmc_track_target &&
                              mc_track_target != fmc_track_source && mc_track_target != fmc_track_target) ||
                             (mc_track_source != fmc_track_source && mc_track_source == fmc_track_target &&
                              mc_track_target != fmc_track_source && mc_track_target != fmc_track_target) ||
                             (mc_track_source != fmc_track_source && mc_track_source != fmc_track_target &&
                              mc_track_target == fmc_track_source && mc_track_target != fmc_track_target) ||
                             (mc_track_source != fmc_track_source && mc_track_source != fmc_track_target &&
                              mc_track_target != fmc_track_source && mc_track_target == fmc_track_target));
            }

            // Below is actually a bit more subtle than just "consecutive in the track of the same motorcycle",
            // due to degenerate track segments.
            bool is_fmc_next_track_segment_after_mc = (mc.id() == fmc.id() && mc_track_target == fmc_track_source);
            bool is_mc_next_track_segment_after_fmc = (mc.id() == fmc.id() && mc_track_source == fmc_track_target);

            // An intersection at the target of the track must crash the motorcycle
            // if the motorcycle reaches this collision point at a later time
            // than another motorcycle. In other words, if there is an intersection at
            // 'mc_track_target', 'mc_track_target' must be the last track entry
            // for 'mc' if the visiting time is lower for the foreign motorcycle.
            if(!is_fmc_next_track_segment_after_mc) // otherwise false positives
            {
              if((mc_track_target == fmc_track_source && time_at_fmc_track_source <= time_at_mc_track_target) ||
                 (mc_track_target == fmc_track_target && time_at_fmc_track_target <= time_at_mc_track_target))
              {
                // must be the last item of the track
                if(mc_tcit != --(mc_track.end()))
                {
                  std::cerr << "Motorcycle: " << std::endl << mc << std::endl;
                  std::cerr << "should have been stopped at: " << std::endl << *mc_track_target << std::endl;
                  std::cerr << "by foreign motorcycle : " << std::endl << fmc << std::endl;
                  if(mc_track_target == fmc_track_source)
                    std::cerr << "times: " << time_at_mc_track_target << " vs " << time_at_fmc_track_source << std::endl;
                  else
                    std::cerr << "times: " << time_at_mc_track_target << " vs " << time_at_fmc_track_target << std::endl;

                  CGAL_assertion(false);
                }
              }
            }

            // If the intersection is at mc's source, then mc must pass first...
            // ... except if it is mc's first track segment and then equal time is allowed
            // as long as it is also fmc's first track segment
            if(!is_mc_next_track_segment_after_fmc) // otherwise false positives
            {
              if(mc_track_source == fmc_track_source)
              {
                if(time_at_mc_track_source == time_at_fmc_track_source &&
                   mc_track_source == mc_track.begin()->source() &&
                   fmc_track_source == fmc_track.begin()->source())
                {
                  // This is the valid case of two motorcycles starting from the same position,
                  // but in different directions.
                }
                else
                {
                  if(time_at_mc_track_source >= time_at_fmc_track_source)
                  {
                    std::cerr << "Motorcycle: " << std::endl << mc << std::endl;
                    std::cerr << "should have been stopped at: " << std::endl << *mc_track_source << std::endl;
                    std::cerr << "by foreign motorcycle : " << std::endl << fmc << std::endl;
                    CGAL_assertion(false);
                  }
                }
              }
              else if(mc_track_source == fmc_track_target)
              {
                if(time_at_mc_track_source >= time_at_fmc_track_target)
                {
                  std::cerr << "Motorcycle: " << std::endl << mc << std::endl;
                  std::cerr << "should have been stopped at: " << std::endl << *mc_track_source << std::endl;
                  std::cerr << "by foreign motorcycle : " << std::endl << fmc << std::endl;
                  CGAL_assertion(false);
                }
              }
            }
          }
        }
      }
    }
  }

  std::cout << "Valid trace" << std::endl;
  return true;
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
  CGAL_assertion(os.good());
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
output_motorcycles_origins_and_destinations() const
{
  std::stringstream oss_orig, oss_dest;
  oss_orig << "results_" << geom_traits().dimension() << "/motorcycles_origins.xyz" << std::ends;
  oss_dest << "results_" << geom_traits().dimension() << "/motorcycles_destinations.xyz" << std::ends;
  std::ofstream oss(oss_orig.str().c_str());
  std::ofstream osd(oss_dest.str().c_str());
  oss.precision(17);
  osd.precision(17);

  MCC_cit mc_it = motorcycles().begin(), mc_end = motorcycles().end();
  for(; mc_it!=mc_end; ++mc_it)
  {
    const Motorcycle& mc = motorcycle(mc_it);

    oss << mc.origin()->point();
    osd << mc.destination()->point();

    if(geom_traits().dimension() == 2) // The '.xyz' format expects 3D points
    {
      oss << " 0";
      osd << " 0";
    }

    oss << '\n';
    osd << '\n';
  }
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
template<typename VertexNodeMap, typename EdgeTrackMap>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
construct_motorcycle_graph(VertexNodeMap& vnmap, EdgeTrackMap& etmap)
{
  trace_graph();
  internal::Motorcycle_graph_builder<Self>(*this)(vnmap, etmap);
  print_motorcycle_graph();

#ifdef CGAL_MOTORCYCLE_GRAPH_OUTPUT
  output_all_points();
#endif

  CGAL_postcondition(is_valid());
}

template<typename MotorcycleGraphTraits, typename MotorcycleType>
void
Motorcycle_graph<MotorcycleGraphTraits, MotorcycleType>::
print_motorcycle_graph() const
{
  std::ofstream out("motorcycle_graph.polylines.txt");
  out.precision(17);

  typedef typename property_map_selector<Halfedge_graph, CGAL::vertex_point_t>::const_type VPMap;
  VPMap vpm = get_const_property_map(boost::vertex_point, graph());

  typename boost::graph_traits<Halfedge_graph>::edge_iterator eit, eend;
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
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_H
