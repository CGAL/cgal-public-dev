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

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_NODE_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_NODE_H

#include <CGAL/assertions.h>
#include <CGAL/functional.h>
#include <CGAL/Polygon_mesh_processing/locate.h>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/container/slist.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>

#include <algorithm>
#include <iostream>
#include <list>
#include <ostream>
#include <utility>

namespace CGAL {

namespace Polyline_tracing {

namespace internal {

// @tmp to easily change types without having to do it in all classes
template<typename Node>
struct Node_ptr_type
{
//  typedef const Node*                                  type
  typedef typename std::set<Node>::iterator            type;
};

template<typename NodePtr>
struct Node_ptr_comparer
  : public CGAL::binary_function<const NodePtr&, const NodePtr&, bool>
{
  // The set corresponds to the same point in different faces --> only compare face descriptors
  bool operator()(const NodePtr& e1, const NodePtr& e2) const {
    return (e1->face() < e2->face());
  }
};

template<typename NodePtr, typename FaceDescriptor>
struct Node_ptr_finder
  : public CGAL::binary_function<const NodePtr&, const FaceDescriptor&, bool>
{
  bool operator()(const NodePtr& e, const FaceDescriptor& fd) {
    return (e->face() < fd);
  }
};

} // namespace internal

// -----------------------------------------------------------------------------
//                               Node_base class
//
// The base contains all the information except the (face, barycentric) location.
// It is dissociated from the 'Node' class because we create multiple nodes
// for points on border of mesh faces (as many as there are incident faces).
// However, these different nodes all share a lot of information, regrouped here.
// -----------------------------------------------------------------------------

// CRTP used to access quickly other nodes corresponding to the same geometric position.
template<typename MotorcycleGraphTraits, typename Derived>
class Motorcycle_graph_node_base
{
  typedef Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>     Self;

  typedef Derived                                                        Node;
  typedef typename internal::Node_ptr_type<Node>::type                   Node_ptr;

public:
  typedef MotorcycleGraphTraits                                          Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                            Triangle_mesh;
  typedef typename Geom_traits::Halfedge_graph                           Halfedge_graph;

  typedef typename Geom_traits::FT                                       FT;
  typedef typename Geom_traits::Point_d                                  Point;
  typedef typename Geom_traits::Face_location                            Face_location;

  typedef typename boost::graph_traits<Halfedge_graph>::vertex_descriptor hg_vertex_descriptor;

  // A container of motorcycles that visit this point. We need to efficiently know:
  // - if a motorcycle visits this point,
  // - the ordered visiting times to detect simultaneous collisions.
  //
  // Note that we need multisets:
  // - for the first container, because a motorcycle might visit a point twice
  //   (e.g. if the trajectory makes a loop and self-intersects at the point);
  // - for the second container, because arrival times might not be unique.
  typedef boost::bimap<
            boost::bimaps::multiset_of<std::size_t>, // multi-set of motorcycles
            boost::bimaps::multiset_of<FT> > // multi-set of visiting times
                                                                         Visiting_motorcycles_container;
  typedef typename Visiting_motorcycles_container::iterator              VMC_it;
  typedef typename Visiting_motorcycles_container::value_type            value_type;
  typedef typename Visiting_motorcycles_container::size_type             size_type;

  typedef typename Visiting_motorcycles_container::left_iterator         VMC_left_it;
  typedef typename Visiting_motorcycles_container::left_const_iterator   VMC_left_cit;
  typedef typename Visiting_motorcycles_container::left_value_type       VMC_left_value_type;

  typedef typename Visiting_motorcycles_container::right_const_iterator  VMC_right_cit;

  // Two reasons for not using an unordered set:
  // - The size of the set is the number of incident faces in the mesh, it's likely to be small (<6)
  // - We need to call find with a different key than 'Node_ptr'
  typedef internal::Node_ptr_comparer<Node_ptr>                          Sibling_comparer;
  typedef std::set<Node_ptr, Sibling_comparer>                           Siblings_container;

  // Constructor
  Motorcycle_graph_node_base(const Point& p)
    : position_(p),
      blocked_(false),
      visiting_mcs_(),
      siblings_(Sibling_comparer()),
      out_vd_(boost::graph_traits<Halfedge_graph>::null_vertex())
  { }

  // Access
  const Point& point() const { return position_; }
  bool is_blocked() const { return blocked_; }
  void block() const { blocked_ = true; }
  std::size_t number_of_visiting_motorcycles() const { return visiting_mcs_.left.size(); }
  Visiting_motorcycles_container& visiting_motorcycles() { return visiting_mcs_; }
  const Visiting_motorcycles_container& visiting_motorcycles() const { return visiting_mcs_; }
  Siblings_container& siblings() { return siblings_; }
  const Siblings_container& siblings() const { return siblings_; }
  hg_vertex_descriptor& graph_vertex() { return out_vd_; }
  const hg_vertex_descriptor& graph_vertex() const { return out_vd_; }

  std::pair<VMC_it, bool> add_motorcycle(const std::size_t id, const FT time) const;
  void add_motorcycles(const Visiting_motorcycles_container& foreign_visiting_mcs) const;

  template<typename MotorcycleOutputIterator>
  MotorcycleOutputIterator earliest_motorcycles(MotorcycleOutputIterator out) const;
  FT earliest_visiting_time() const;

  // Checks if the motorcycle 'id' visits.
  VMC_left_it find_motorcycle(const std::size_t id) const;

  // Checks if the motorcycle 'id' visits at time 'visiting_time'.
  VMC_left_it find_motorcycle(const std::size_t id, const FT visiting_time) const;

  // Checks if the motorcycle 'id' visits at time between 'min_' and 'max_visiting_time'.
  // The last parameter is optional and can be used to grab the visiting time
  // 'strictly_at_X' to include the boundary of the interval or not.
  VMC_left_it find_motorcycle(const std::size_t id,
                              const FT min_visiting_time, const FT max_visiting_time,
                              FT& visiting_time /*time at which it is visited*/,
                              const bool strictly_at_min = false,  const bool strictly_at_max = false) const;
  VMC_left_it find_motorcycle(const std::size_t id,
                              const FT min_visiting_time, const FT max_visiting_time,
                              const bool strictly_at_min = false, const bool strictly_at_max = false) const;

  // Checks if a motorcycle visits the point (possibly within a given time interval)
  bool has_motorcycle(const std::size_t id) const;
  bool has_motorcycle(const std::size_t id, const FT visiting_time) const;
  bool has_motorcycle(const std::size_t id,
                      const FT min_visiting_time, const FT max_visiting_time,
                      FT& visiting_time /*time at which it is visited*/,
                      const bool strictly_at_min = false, const bool strictly_at_max = false) const;
  bool has_motorcycle(const std::size_t id,
                      const FT min_visiting_time, const FT max_visiting_time,
                      const bool strictly_at_min = false, const bool strictly_at_max = false) const;
  bool has_motorcycles() const;

  size_type remove_motorcycle(const std::size_t id) const;

  // Output
  friend std::ostream& operator<<(std::ostream& out, const Self& dec)
  {
    out << "  Point: (" << dec.point() << ") blocked: " << dec.is_blocked() << std::endl;
    out << "  Visiting motorcycles:" << std::endl;
    VMC_right_cit vmc_it = dec.visiting_motorcycles().right.begin();
    VMC_right_cit end = dec.visiting_motorcycles().right.end();
    for(; vmc_it!=end; ++vmc_it)
      out << "\t motorcycle #" << vmc_it->second << " time: " << vmc_it->first << std::endl;

    out << "  Siblings:" << std::endl;
    typename Siblings_container::const_iterator smcit = dec.siblings().begin();
    for(; smcit!=dec.siblings().end(); ++smcit)
    {
      std::cout << "\t location fd: " << (*smcit)->face() << " / bc: ["
                                      << (*smcit)->barycentric_coordinate(0) << " "
                                      << (*smcit)->barycentric_coordinate(1) << " "
                                      << (*smcit)->barycentric_coordinate(2) << "]" << std::endl;;
    }

    return out;
  }

private:
  // This class is meant to be an element of a set, which must be const.
  // However, some of its members must be modifiable (e.g. visiting motorcycles)
  // so they are made mutable. It is safe because the set is only ordered
  // with the location, which is never modified.

  // The position of the point (only for information)
  const Point position_;

  // Whether the position is blocked or not
  mutable bool blocked_;

  // The motorcycles visiting this point
  mutable Visiting_motorcycles_container visiting_mcs_;

  // When a point is on the border of a face, multiple 'Face_location' represent
  // the same point. These locations are referred to as 'siblings'.
  // @todo make it optional
  mutable Siblings_container siblings_;

  // Equivalent vertex in the output graph
  mutable hg_vertex_descriptor out_vd_;
};

template<typename MotorcycleGraphTraits, typename Derived>
std::pair<typename Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::VMC_it, bool>
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
add_motorcycle(const std::size_t id, const FT time) const
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << " > Point " << this << " is visited by motorcycle #" << id << " at time: " << time << std::endl;
#endif
  CGAL_expensive_precondition(!has_motorcycle(id, time));
  return visiting_mcs_.insert(value_type(id, time));
}

template<typename MotorcycleGraphTraits, typename Derived>
void
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
add_motorcycles(const Visiting_motorcycles_container& foreign_visiting_mcs) const
{
  VMC_left_cit vmc_it = foreign_visiting_mcs.left.begin();
  VMC_left_cit end = foreign_visiting_mcs.left.end();
  for(; vmc_it!=end; ++vmc_it)
  {
    // Although it's a multiset, we don't want to insert redundant information
    // in the visiting motorcycle container
    std::pair<VMC_left_it, bool> is_insert_successful =
      visiting_mcs_.left.insert(VMC_left_value_type(vmc_it->first, vmc_it->second));
    CGAL_expensive_precondition(is_insert_successful.second);
  }
}

template<typename MotorcycleGraphTraits, typename Derived>
template<typename MotorcycleOutputIterator>
MotorcycleOutputIterator
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
earliest_motorcycles(MotorcycleOutputIterator out) const
{
  if(visiting_motorcycles().empty())
    return out;

  VMC_right_cit mc_it = visiting_motorcycles().right.begin(),
                mc_end = visiting_motorcycles().right.end();
  const FT earliest_visiting_time = mc_it->first;

  do
  {
    *out++ = mc_it->second;
    ++mc_it;
  }
  while(mc_it != mc_end && mc_it->first == earliest_visiting_time);

  return out;
}

template<typename MotorcycleGraphTraits, typename Derived>
typename Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::FT
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
earliest_visiting_time() const
{
  CGAL_precondition(!visiting_motorcycles().empty());
  return visiting_motorcycles().right.begin()->first;
}

template<typename MotorcycleGraphTraits, typename Derived>
typename Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::VMC_left_it
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
find_motorcycle(const std::size_t id) const
{
  // Since 'lower_bound' is used, it returns the first motorcycle fitting this
  VMC_left_it it = visiting_mcs_.left.lower_bound(id);
  if(it == visiting_mcs_.left.end() || it->first != id)
    return visiting_mcs_.left.end();

  return it;
}

template<typename MotorcycleGraphTraits, typename Derived>
typename Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::VMC_left_it
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
find_motorcycle(const std::size_t id, const FT visiting_time) const
{
  VMC_left_it mit = find_motorcycle(id);
  if(mit == visiting_motorcycles().left.end())
    return mit;

  bool is_valid_iterator = true;
  while(is_valid_iterator)
  {
    CGAL_assertion(mit->first == id);
    if(mit->second == visiting_time)
      return mit;

    ++mit;
    is_valid_iterator = (mit != visiting_motorcycles().left.end() && mit->first == id);
  }

  return visiting_mcs_.left.end();
}

template<typename MotorcycleGraphTraits, typename Derived>
typename Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::VMC_left_it
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
find_motorcycle(const std::size_t id,
                const FT min_visiting_time, const FT max_visiting_time,
                FT& visiting_time,
                const bool strictly_at_min, const bool strictly_at_max) const
{
  CGAL_precondition(min_visiting_time <= max_visiting_time);

  VMC_left_it mit = find_motorcycle(id);
  if(mit == visiting_motorcycles().left.end())
    return mit;

  bool is_valid_iterator = true;
  while(is_valid_iterator) // while still considering the motorcycle 'id'
  {
    CGAL_assertion(mit->first == id);
    const FT time = mit->second;

    if((time > min_visiting_time || (!strictly_at_min && time == min_visiting_time)) &&
       (time < max_visiting_time || (!strictly_at_max && time == max_visiting_time)))
    {
      visiting_time = time;
      return mit;
    }

    ++mit;
    is_valid_iterator = (mit != visiting_motorcycles().left.end() && mit->first == id);
  }

  return visiting_mcs_.left.end();
}

template<typename MotorcycleGraphTraits, typename Derived>
typename Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::VMC_left_it
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
find_motorcycle(const std::size_t id,
                const FT min_visiting_time, const FT max_visiting_time,
                const bool strictly_at_min, const bool strictly_at_max) const
{
  FT useless;
  return find_motorcycle(id, min_visiting_time, max_visiting_time, useless,
                         strictly_at_min, strictly_at_max);
}

template<typename MotorcycleGraphTraits, typename Derived>
bool
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
has_motorcycle(const std::size_t id) const
{
  return (find_motorcycle(id) != visiting_motorcycles().left.end());
}

template<typename MotorcycleGraphTraits, typename Derived>
bool
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
has_motorcycle(const std::size_t id, const FT visiting_time) const
{
  return (find_motorcycle(id, visiting_time) != visiting_motorcycles().left.end());
}

template<typename MotorcycleGraphTraits, typename Derived>
bool
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
has_motorcycle(const std::size_t id,
               const FT min_visiting_time, const FT max_visiting_time,
               FT& visiting_time,
               const bool strictly_at_min, const bool strictly_at_max) const
{
  CGAL_precondition(min_visiting_time <= max_visiting_time);
  return (find_motorcycle(id, min_visiting_time, max_visiting_time, visiting_time,
                          strictly_at_min, strictly_at_max) != visiting_motorcycles().left.end());
}

template<typename MotorcycleGraphTraits, typename Derived>
bool
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
has_motorcycle(const std::size_t id,
               const FT min_visiting_time, const FT max_visiting_time,
               const bool strictly_at_min, const bool strictly_at_max) const
{
  CGAL_precondition(min_visiting_time <= max_visiting_time);

  FT useless;
  return (find_motorcycle(id, min_visiting_time, max_visiting_time, useless,
                          strictly_at_min, strictly_at_max) != visiting_motorcycles().left.end());
}

template<typename MotorcycleGraphTraits, typename Derived>
bool
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
has_motorcycles() const
{
  return !visiting_motorcycles().empty();
}

template<typename MotorcycleGraphTraits, typename Derived>
typename Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::size_type
Motorcycle_graph_node_base<MotorcycleGraphTraits, Derived>::
remove_motorcycle(const std::size_t id) const
{
  // Note that this will remove all instances of motorcycle 'id' in the multiset
  return visiting_mcs_.left.erase(id);
}


// -----------------------------------------------------------------------------
//                                    Node class
//
// This class represents a point that is involved in the motorcycle graph algorithm
//
// It is useful for robustness to regroup them in a single dictionary-like structure,
// sorted by location because intersections are computed on locations (instead
// of points)
// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits>
class Motorcycle_graph_node
{
  typedef Motorcycle_graph_node<MotorcycleGraphTraits>                   Self;
  typedef Motorcycle_graph_node_base<MotorcycleGraphTraits, Self>        Base;

  typedef typename internal::Node_ptr_type<Self>::type                   Node_ptr;

  // Node bases container
  typedef boost::container::slist<Base>                                  NB_container;
  typedef typename NB_container::iterator                                NBC_it;

public:
  typedef MotorcycleGraphTraits                                          Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                            Triangle_mesh;
  typedef typename Geom_traits::Halfedge_graph                           Halfedge_graph;

  typedef typename Geom_traits::FT                                       FT;
  typedef typename Geom_traits::Point_d                                  Point;

  typedef typename Geom_traits::Face_location                            Face_location;
  typedef typename Geom_traits::Barycentric_coordinates                  Barycentric_coordinates;

  typedef typename boost::graph_traits<Halfedge_graph>::vertex_descriptor hg_vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor   face_descriptor;

  typedef typename Base::size_type                                       size_type;
  typedef typename Base::Visiting_motorcycles_container                  Visiting_motorcycles_container;
  typedef typename Base::VMC_it                                          VMC_it;
  typedef typename Base::VMC_left_it                                     VMC_left_it;
  typedef typename Base::VMC_right_cit                                   VMC_right_cit;

  typedef typename Base::Siblings_container                              Siblings_container;
  typedef typename Siblings_container::iterator                          SC_it;

  // Access
  const Face_location& location() const { return location_; }
  const face_descriptor face() const { return location_.first; }
  const Barycentric_coordinates& barycentric_coordinates() const { return location_.second; }
  const FT barycentric_coordinate(int i) const { CGAL_precondition(i>=0 && i<3); return location_.second[i]; }

    // the "base" is a pointer to an object of type 'Base' that is stored in a container
    // (because the base is common to multiple locations)
  NBC_it& base() const { CGAL_precondition(base_ != NBC_it()); return base_; }
  void set_base(NBC_it new_base) const { base_ = new_base; }

    // return the point's location in another face
    // \pre the location is on the border of 'fd'
  const Face_location& location(face_descriptor fd) const;

  // Constructor
  Motorcycle_graph_node(const Face_location& location) : location_(location), base_() { }

  // ---------------------------------------------------------------------------
  // Simple wrappers to artifically create a base
  const Point& point() const { return base()->point(); }
  bool is_blocked() const { return base()->is_blocked(); }
  void block() const { base()->block(); }
  std::size_t number_of_visiting_motorcycles() const { return base()->number_of_visiting_motorcycles(); }
  Visiting_motorcycles_container& visiting_motorcycles() { return base()->visiting_motorcycles(); }
  const Visiting_motorcycles_container& visiting_motorcycles() const { return base()->visiting_motorcycles(); }
  Siblings_container& siblings() const { return base()->siblings(); }
  hg_vertex_descriptor& graph_vertex() const { return base()->graph_vertex(); }

  std::pair<VMC_it, bool> add_motorcycle(const std::size_t id, const FT time) const {
    return base()->add_motorcycle(id, time);
  }
  void add_motorcycles(const Visiting_motorcycles_container& foreign_visiting_mcs) const {
    return base()->add_motorcycles(foreign_visiting_mcs);
  }

  template<typename MotorcycleOutputIterator>
 MotorcycleOutputIterator earliest_motorcycles(MotorcycleOutputIterator out) const { return base()->earliest_motorcycles(out); }
  FT earliest_visiting_time() const { return base()->earliest_visiting_time(); }

  VMC_left_it find_motorcycle(const std::size_t id) const { return base()->find_motorcycle(id); }
  // Check if the motorcycle 'id' visits at time 'visiting_time'.
  VMC_left_it find_motorcycle(const std::size_t id, const FT visiting_time) const {
    return base()->find_motorcycle(id, visiting_time);
  }

  VMC_left_it find_motorcycle(const std::size_t id, const FT min_visiting_time, const FT max_visiting_time,
                              FT& visiting_time, const bool strictly_at_min = false, const bool strictly_at_max = false) const {
    return base()->find_motorcycle(id, min_visiting_time, max_visiting_time, visiting_time, strictly_at_min, strictly_at_max);
  }
  VMC_left_it find_motorcycle(const std::size_t id, const FT min_visiting_time, const FT max_visiting_time,
                              const bool strictly_at_min = false, const bool strictly_at_max = false) const {
    return base()->find_motorcycle(id, min_visiting_time, max_visiting_time, strictly_at_min, strictly_at_max);
  }

  bool has_motorcycle(const std::size_t id) const { return base()->has_motorcycle(id); }
  bool has_motorcycle(const std::size_t id, const FT visiting_time) const { return base()->has_motorcycle(id, visiting_time); }
  bool has_motorcycle(const std::size_t id,
                      const FT min_visiting_time, const FT max_visiting_time,
                      FT& visiting_time /*time at which it is visited*/,
                      const bool strictly_at_min = false, const bool strictly_at_max = false) const {
    return base()->has_motorcycle(id, min_visiting_time, max_visiting_time, visiting_time, strictly_at_min, strictly_at_max);
  }
  bool has_motorcycle(const std::size_t id,
                      const FT min_visiting_time, const FT max_visiting_time,
                      const bool strictly_at_min = false, const bool strictly_at_max = false) const {
    return base()->has_motorcycle(id, min_visiting_time, max_visiting_time, strictly_at_min, strictly_at_max);
  }
  bool has_motorcycles() const { return base()->has_motorcycles(); }

  // check if the two earliest motorcycles meet at the same time
  size_type remove_motorcycle(const std::size_t id) const { return base()->remove_motorcycle(id); }
  // ---------------------------------------------------------------------------

  // Functions
  SC_it lower_bound(face_descriptor fd) const
  {
    // Use lower bound and a custom comparer to avoid creating a dummy iterator
    internal::Node_ptr_finder<Node_ptr, face_descriptor> comp;
    return std::lower_bound(siblings().begin(), siblings().end(), fd, comp);
  }

  Node_ptr sibling(face_descriptor fd) const
  {
    CGAL_precondition(!siblings().empty());
    SC_it it = lower_bound(fd);

    // Sibling must be present
    CGAL_assertion(it != siblings().end());
    CGAL_assertion((*it)->face() == fd);

    return *it;
  }

  bool is_sibling(const Node_ptr& other_node) const
  {
    CGAL_expensive_assertion_code(if(other_node->point() == point()))
    CGAL_expensive_assertion(is_sibling(other_node->location()));

    return (other_node->point() == point());
  }

  bool is_sibling(const Face_location& other_location) const
  {
    if(other_location == location())
      return true;

    if(siblings().empty())
      return false;

    SC_it it = lower_bound(other_location.first);
    if(it == siblings().end() || (*it)->face() != other_location.first)
      return false;

    return true;
  }

  // To build a set<Node>
  friend bool operator<(const Self& le, const Self& re)
  {
    if(le.face() == re.face())
    {
      // If the faces are equal, lexicographically compare the barycentric coordinates
      return std::lexicographical_compare(le.barycentric_coordinates().begin(), le.barycentric_coordinates().end(),
                                          re.barycentric_coordinates().begin(), re.barycentric_coordinates().end());
    }

    return le.face() < re.face();
  }

  // To build an unordered_set<Node>
  friend bool operator==(const Self& le, const Self& re)
  {
    return (le.face() == re.face() &&
            le.location().second == re.location().second);
  }

  friend std::size_t hash_value(const Self& dec)
  {
    boost::hash<face_descriptor> face_hasher;
    std::size_t seed = 0;
    boost::hash_combine(seed, face_hasher(dec.face()));
    boost::hash_combine(seed, boost::hash_range(dec.location().second.begin(),
                                                dec.location().second.end()));
    return seed;
  }

  // Output
  friend std::ostream& operator<<(std::ostream& out, const Self& dec)
  {
    out << "  Location: " << dec.face() << " barycentric coordinates: { "
        << dec.barycentric_coordinate(0) << "; " << dec.location().second[1]
        << "; " << dec.location().second[2] << "}" << std::endl;
    out << *(dec.base());
    return out;
  }

private:
  const Face_location location_; // location in the mesh
  mutable NBC_it base_; // everything else (mutable because this class is the value type of a set)
};

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_NODE_H
