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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_DICTIONARY_H
#define CGAL_POLYLINE_TRACING_DICTIONARY_H

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/Polygon_mesh_processing/locate.h>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/functional/hash.hpp>
#include <boost/next_prior.hpp>
#include <boost/unordered_set.hpp>

#include <iostream>
#include <list>
#include <map>
#include <ostream>
#include <utility>
#include <vector>

namespace CGAL {

namespace Polyline_tracing {

namespace internal {

template<typename Face_location>
struct Face_location_comparer
{
  bool operator()(const Face_location& fl1, const Face_location& fl2) const
  {
    return fl1.first < fl2.first;
  }
};

} // namespace internal

// -----------------------------------------------------------------------------
//                        Dictionary_entry_base class
// The base contains all the information except the (face, barycentric) location.
// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits_>
class Dictionary_entry_base
{
  typedef Dictionary_entry_base<MotorcycleGraphTraits_>                  Self;

public:
  typedef MotorcycleGraphTraits_                                         Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                            Triangle_mesh;

  typedef typename Geom_traits::FT                                       FT;
  typedef typename Geom_traits::Point_d                                  Point;
  typedef typename Geom_traits::Face_location                            Face_location;

  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor   face_descriptor;

  // A container of motorcycles that visit this point. We need to efficiently know:
  // - if a motorcycle visits this point,
  // - the ordered visiting times to detect simultaneous collisions.
  //
  // Note that we need multisets:
  // - for the first container, because a motorcycle might visit a point twice
  //   (e.g. if the trajectory makes a loop and self-intersects at the point);
  // - for the second container, because arrival times might not be unique.
  typedef boost::bimap<
            boost::bimaps::multiset_of<std::size_t>, // set of motorcycles
            boost::bimaps::multiset_of<FT> > // multi-set of visiting times
                                                                         Visiting_motorcycles_container;
  typedef typename Visiting_motorcycles_container::iterator              VMC_it;
  typedef typename Visiting_motorcycles_container::value_type            value_type;
  typedef typename Visiting_motorcycles_container::size_type             size_type;

  typedef typename Visiting_motorcycles_container::left_iterator         VMC_left_it;
  typedef typename Visiting_motorcycles_container::left_const_iterator   VMC_left_cit;
  typedef typename Visiting_motorcycles_container::left_value_type       VMC_left_value_type;

  typedef typename Visiting_motorcycles_container::right_const_iterator  VMC_right_cit;

  typedef internal::Face_location_comparer<Face_location>                FL_comparer;
  typedef std::set<Face_location, FL_comparer>                           Siblings_container;

  // Constructor
  Dictionary_entry_base(const Point& p)
    : position_(p), blocked_(false), visiting_mcs_(), siblings_(FL_comparer())
  { }

  // Access
  const Point& point() const { return position_; }
  bool is_blocked() const { return blocked_; }
  void block() const { blocked_ = true; }
  Visiting_motorcycles_container& visiting_motorcycles() { return visiting_mcs_; }
  const Visiting_motorcycles_container& visiting_motorcycles() const { return visiting_mcs_; }
  Siblings_container& siblings() { return siblings_; }
  const Siblings_container& siblings() const { return siblings_; }

  std::pair<VMC_it, bool> add_motorcycle(const std::size_t id, const FT time) const;
  void add_motorcycles(const Visiting_motorcycles_container& foreign_visiting_mcs) const;

  // Returns an iterator to the motorcycle that first visits the point.
  VMC_right_cit earliest_motorcycle() const;

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

  // Checks if the two earliest motorcycles meet at the same time
  bool has_simultaneous_collision() const;

  size_type remove_motorcycle(const std::size_t id) const;

  // Output
  friend std::ostream& operator<<(std::ostream& out, const Self& dec)
  {
    out << "Point: (" << dec.point() << ") blocked: " << dec.is_blocked() << std::endl;
    out << "  Visiting motorcycles:" << std::endl;
    VMC_left_cit vmc_it = dec.visiting_motorcycles().left.begin();
    VMC_left_cit end = dec.visiting_motorcycles().left.end();
    for(; vmc_it!=end; ++vmc_it)
      out << "\t motorcycle #" << vmc_it->first << " time: " << vmc_it->second << std::endl;

    out << "  Siblings:" << std::endl;
    typename Siblings_container::const_iterator smcit = dec.siblings().begin();
    for(; smcit!=dec.siblings().end(); ++smcit)
    {
      std::cout << "\t location is fd: " << smcit->first << " / bc: ["
                                         << smcit->second[0] << " "
                                         << smcit->second[1] << " "
                                         << smcit->second[2] << "]" << std::endl;;
    }

    return out;
  }

private:
  // This class is meant to be an element of a set entry, which is const. However,
  // its members must still be modifiable so they are made mutable.
  // It is safe because the set is only ordered with the location, which is never modified.

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
};

template<typename MotorcycleGraphTraits_>
std::pair<typename Dictionary_entry_base<MotorcycleGraphTraits_>::VMC_it, bool>
Dictionary_entry_base<MotorcycleGraphTraits_>::
add_motorcycle(const std::size_t id, const FT time) const
{
  std::cout << " > Point " << this << " is visited by motorcycle #" << id << " at time: " << time << std::endl;
  CGAL_expensive_precondition(!has_motorcycle(id, time));
  return visiting_mcs_.insert(value_type(id, time));
}

template<typename MotorcycleGraphTraits_>
void
Dictionary_entry_base<MotorcycleGraphTraits_>::
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

template<typename MotorcycleGraphTraits_>
typename Dictionary_entry_base<MotorcycleGraphTraits_>::VMC_right_cit
Dictionary_entry_base<MotorcycleGraphTraits_>::
earliest_motorcycle() const
{
  CGAL_precondition(!visiting_mcs_.empty());
  return visiting_motorcycles().right.begin();
}

template<typename MotorcycleGraphTraits_>
typename Dictionary_entry_base<MotorcycleGraphTraits_>::VMC_left_it
Dictionary_entry_base<MotorcycleGraphTraits_>::
find_motorcycle(const std::size_t id) const
{
  // Since 'lower_bound' is used, it returns the first motorcycle fitting this
  VMC_left_it it = visiting_mcs_.left.lower_bound(id);
  if(it == visiting_mcs_.left.end() || it->first != id)
    return visiting_mcs_.left.end();

  return it;
}

template<typename MotorcycleGraphTraits_>
typename Dictionary_entry_base<MotorcycleGraphTraits_>::VMC_left_it
Dictionary_entry_base<MotorcycleGraphTraits_>::
find_motorcycle(const std::size_t id, const FT visiting_time) const
{
  VMC_left_it mit = find_motorcycle(id);
  if(mit == visiting_mcs_.left.end())
    return mit;

  bool is_valid_iterator = true;
  while(is_valid_iterator)
  {
    CGAL_assertion(mit->first == id);
    if(mit->second == visiting_time)
      return mit;

    ++mit;
    is_valid_iterator = (mit != visiting_mcs_.left.end() && mit->first == id);
  }

  return visiting_mcs_.left.end();
}

template<typename MotorcycleGraphTraits_>
typename Dictionary_entry_base<MotorcycleGraphTraits_>::VMC_left_it
Dictionary_entry_base<MotorcycleGraphTraits_>::
find_motorcycle(const std::size_t id,
                const FT min_visiting_time, const FT max_visiting_time,
                FT& visiting_time,
                const bool strictly_at_min, const bool strictly_at_max) const
{
  CGAL_precondition(min_visiting_time <= max_visiting_time);

  VMC_left_it mit  = find_motorcycle(id);
  if(mit == visiting_mcs_.left.end())
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
    is_valid_iterator = (mit != visiting_mcs_.left.end() && mit->first == id);
  }

  return visiting_mcs_.left.end();
}

template<typename MotorcycleGraphTraits_>
typename Dictionary_entry_base<MotorcycleGraphTraits_>::VMC_left_it
Dictionary_entry_base<MotorcycleGraphTraits_>::
find_motorcycle(const std::size_t id,
                const FT min_visiting_time, const FT max_visiting_time,
                const bool strictly_at_min, const bool strictly_at_max) const
{
  FT useless;
  return find_motorcycle(id, min_visiting_time, max_visiting_time, useless,
                         strictly_at_min, strictly_at_max);
}

template<typename MotorcycleGraphTraits_>
bool
Dictionary_entry_base<MotorcycleGraphTraits_>::
has_motorcycle(const std::size_t id) const
{
  return (find_motorcycle(id) != visiting_mcs_.left.end());
}

template<typename MotorcycleGraphTraits_>
bool
Dictionary_entry_base<MotorcycleGraphTraits_>::
has_motorcycle(const std::size_t id, const FT visiting_time) const
{
  return (find_motorcycle(id, visiting_time) != visiting_mcs_.left.end());
}

template<typename MotorcycleGraphTraits_>
bool
Dictionary_entry_base<MotorcycleGraphTraits_>::
has_motorcycle(const std::size_t id,
               const FT min_visiting_time, const FT max_visiting_time,
               FT& visiting_time,
               const bool strictly_at_min, const bool strictly_at_max) const
{
  CGAL_precondition(min_visiting_time <= max_visiting_time);
  return (find_motorcycle(id, min_visiting_time, max_visiting_time, visiting_time,
                         strictly_at_min, strictly_at_max) != visiting_mcs_.left.end());
}

template<typename MotorcycleGraphTraits_>
bool
Dictionary_entry_base<MotorcycleGraphTraits_>::
has_motorcycle(const std::size_t id,
               const FT min_visiting_time, const FT max_visiting_time,
               const bool strictly_at_min, const bool strictly_at_max) const
{
  CGAL_precondition(min_visiting_time <= max_visiting_time);

  FT useless;
  return (find_motorcycle(id, min_visiting_time, max_visiting_time, useless,
                         strictly_at_min, strictly_at_max) != visiting_mcs_.left.end());
}

template<typename MotorcycleGraphTraits_>
bool
Dictionary_entry_base<MotorcycleGraphTraits_>::
has_motorcycles() const
{
  return !visiting_mcs_.empty();
}

template<typename MotorcycleGraphTraits_>
bool
Dictionary_entry_base<MotorcycleGraphTraits_>::
has_simultaneous_collision() const
{
  CGAL_precondition(!visiting_motorcycles().empty());
  if(visiting_motorcycles().size() == 1)
    return false;

  // accessing 'right' of the bimap gives ordered time values
  // the first and second entries are thus the times for the two closest
  // points. If the times are equal (or almost), there is a simultaneous collision.
  VMC_right_cit first_mc_it = visiting_motorcycles().right.begin();
  VMC_right_cit second_mc_it = ++(visiting_motorcycles().right.begin());
  const FT first_time = first_mc_it->first;
  const FT second_time = second_mc_it->first;
  CGAL_assertion(first_time <= second_time);

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "closest motorcycles at times: " << first_time << " and " << second_time << std::endl;
#endif

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  // Add a little bit of tolerance
  return (CGAL::abs(second_time - first_time) < (std::numeric_limits<FT>::epsilon() *
                                                 CGAL::abs(first_time + second_time)));
#else
  return (second_time == first_time);
#endif
}

template<typename MotorcycleGraphTraits_>
typename Dictionary_entry_base<MotorcycleGraphTraits_>::size_type
Dictionary_entry_base<MotorcycleGraphTraits_>::
remove_motorcycle(const std::size_t id) const
{
  // Note that this will remove all instances of motorcycle 'id' in the multiset
  return visiting_mcs_.left.erase(id);
}

// -----------------------------------------------------------------------------
//                           Dictionary entry class
// This class represents a point that is involved in the motorcycle graph algorithm
//
// It is useful for robustness to regroup them in a single dictionary-like structure,
// sorted by location because intersections are computed on locations (instead
// of points)
// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits_,
         typename DictionaryEntryBaseContainer_ =
           std::list<Dictionary_entry_base<MotorcycleGraphTraits_> > >
class Dictionary_entry
{
  typedef Dictionary_entry<MotorcycleGraphTraits_, DictionaryEntryBaseContainer_>  Self;
  typedef Dictionary_entry_base<MotorcycleGraphTraits_>                            Base;
  typedef typename DictionaryEntryBaseContainer_::iterator                         DEBC_it;

public:
  typedef MotorcycleGraphTraits_                                         Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                            Triangle_mesh;

  typedef typename Geom_traits::FT                                       FT;
  typedef typename Geom_traits::Point_d                                  Point;

  typedef typename Geom_traits::Face_location                            Face_location;
  typedef typename Geom_traits::Barycentric_coordinates                  Barycentric_coordinates;

  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

  typedef typename Base::size_type                                       size_type;
  typedef typename Base::Visiting_motorcycles_container                  Visiting_motorcycles_container;
  typedef typename Base::VMC_it                                          VMC_it;
  typedef typename Base::VMC_left_it                                     VMC_left_it;
  typedef typename Base::VMC_right_cit                                   VMC_right_cit;
  typedef typename Base::Siblings_container                              Siblings_container;

  // Access
  const Face_location& location() const { return location_; }

    // the "base" is a pointer to an object of type 'Base' that is stored in a container
    // (because the base is common to multiple locations)
  DEBC_it& base() const { CGAL_precondition(base_ != DEBC_it()); return base_; }
  void set_base(DEBC_it new_base) const { base_ = new_base; }

    // return the point's location in another face
    // \pre the location is on the border of 'fd'
  const Face_location& location(face_descriptor fd) const;

  // Constructor
  Dictionary_entry(const Face_location& location) : location_(location), base_() { }

  // ---------------------------------------------------------------------------
  // Simple wrappers to artifically create a base
  const Point& point() const { return base()->point(); }
  bool is_blocked() const { return base()->is_blocked(); }
  void block() const { base()->block(); }
  Visiting_motorcycles_container& visiting_motorcycles() { return base()->visiting_motorcycles(); }
  const Visiting_motorcycles_container& visiting_motorcycles() const { return base()->visiting_motorcycles(); }
  Siblings_container& siblings() const { return base()->siblings(); }

  std::pair<VMC_it, bool> add_motorcycle(const std::size_t id, const FT time) const {
    return base()->add_motorcycle(id, time);
  }
  void add_motorcycles(const Visiting_motorcycles_container& foreign_visiting_mcs) const {
    return base()->add_motorcycles(foreign_visiting_mcs);
  }

  VMC_right_cit earliest_motorcycle() const { return base()->earliest_motorcycle(); }
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
  bool has_simultaneous_collision() const { return base()->has_simultaneous_collision(); }
  size_type remove_motorcycle(const std::size_t id) const { return base()->remove_motorcycle(id); }
  // ---------------------------------------------------------------------------

  // Functions
  const Face_location& sibling(face_descriptor fd) const
  {
    CGAL_precondition(!siblings().empty());
    if(location().first == fd)
      return location();

    Face_location dummy(std::make_pair(fd, Barycentric_coordinates()));
    typename Siblings_container::iterator it = siblings().find(dummy);
    CGAL_postcondition(it != siblings().end()); // couldn't find the sibling, is 'fd' correct ?
    return *it;
  }

  bool is_sibling(const Face_location& location) const
  {
    CGAL_precondition(!siblings().empty());
    return (siblings().find(location) != siblings().end());
  }

  // To build a set<Dictionary_entry> (now obsolete since an unordered_set is now used,
  // but leaving it available).
  friend bool operator<(const Self& lhs, const Self& rhs)
  {
    if(lhs.location().first == rhs.location().first)
    {
      // If the faces are equal, lexicographically compare the barycentric coordinates
      return std::lexicographical_compare(lhs.location().second.begin(), lhs.location().second.end(),
                                          rhs.location().second.begin(), rhs.location().second.end());
    }

    return lhs.location().first < rhs.location().first;
  }

  friend bool operator==(const Self& lhs, const Self& rhs)
  {
    return (lhs.location().first == rhs.location().first &&
            lhs.location().second == rhs.location().second);
  }

  friend std::size_t hash_value(const Self& dec)
  {
    boost::hash<face_descriptor> face_hasher;
    std::size_t seed = 0;
    boost::hash_combine(seed, face_hasher(dec.location().first));
    boost::hash_combine(seed, boost::hash_range(dec.location().second.begin(),
                                                dec.location().second.end()));
    return seed;
  }

  // Output
  friend std::ostream& operator<<(std::ostream& out, const Self& dec)
  {
    out << "  Location: " << dec.location().first << " barycentric coordinates: { "
        << dec.location().second[0] << "; " << dec.location().second[1]
        << "; " << dec.location().second[2] << "}" << std::endl;
    out << *(dec.base());
    return out;
  }

private:
  const Face_location location_; // location in the mesh
  mutable DEBC_it base_; // everything else (mutable because this class is the value type of a set)
};

// -----------------------------------------------------------------------------
//                           Dictionary class
// The points involved in the motorcycle graph.
// -----------------------------------------------------------------------------
template<typename MotorcycleGraphTraits_>
class Dictionary
{
  typedef Dictionary<MotorcycleGraphTraits_>                        Self;

public:
  typedef MotorcycleGraphTraits_                                    Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                       Triangle_mesh;

  typedef typename Geom_traits::FT                                  FT;
  typedef typename Geom_traits::Point_d                             Point;

  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;
  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                   descriptor_variant;

  typedef Dictionary_entry_base<MotorcycleGraphTraits_>             Dictionary_entry_base;
  typedef std::list<Dictionary_entry_base>                          DEB_container;
  typedef typename DEB_container::iterator                          DEBC_it;
  typedef typename DEB_container::const_iterator                    DEBC_cit;
  typedef Dictionary_entry<MotorcycleGraphTraits_, DEB_container>   Dictionary_entry;

  typedef typename Dictionary_entry::Face_location                  Face_location;

  typedef boost::unordered_set<Dictionary_entry>                    Dictionary_entry_container;
  typedef typename Dictionary_entry_container::iterator             DEC_it;

  // Access
  Dictionary_entry_container& entries() { return entries_; }
  const Dictionary_entry_container& entries() const { return entries_; }

  // Constructor
  Dictionary() : entries_(), entries_bases_() { }

  // Functions
  void erase(DEC_it pos, const bool erase_siblings = true);

  std::pair<DEC_it, bool> find(const Dictionary_entry& e) const;
  std::pair<DEC_it, bool> find(const Face_location& loc) const;

  std::pair<DEC_it, bool> get_sibling(const DEC_it e, const face_descriptor fd) const;

  std::pair<DEC_it, bool> insert(const Face_location& loc, const Point& p, const Triangle_mesh& mesh);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const Point& p, const std::size_t i, const FT time, const Triangle_mesh& mesh);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const std::size_t i, const FT time, const Triangle_mesh& mesh);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const Triangle_mesh& mesh);

  // Ouput
  friend std::ostream& operator<<(std::ostream& out, const Self& d)
  {
    DEC_it it = d.entries().begin(), end = d.entries().end();
    for(; it!=end; ++it)
      out << &*it << " " << *it << std::endl;

    return out;
  }

private:
  Dictionary_entry_container entries_;
  DEB_container entries_bases_;
};

template<typename MotorcycleGraphTraits_>
void
Dictionary<MotorcycleGraphTraits_>::
erase(DEC_it pos, const bool erase_siblings)
{
  // WARNING: this function also erase (by default) the potential entry's siblings !!

  face_descriptor fd = pos->location().first;

  DEBC_it common_base = pos->base();
  CGAL_assertion(common_base != DEBC_it());

  if(erase_siblings)
  {
    // remove the other siblings, if any
    typename Dictionary_entry_base::Siblings_container::iterator smit = common_base->siblings().begin();
    for(; smit!=common_base->siblings().end(); ++smit)
    {
      // Location for that face
      face_descriptor adj_fd = smit->first;
      CGAL_assertion(adj_fd != boost::graph_traits<Triangle_mesh>::null_face());

      if(adj_fd == fd)
        continue;

      // Find it in the entries (failure is not an option), and erase it
      std::pair<DEC_it, bool> entry = find(*smit);
      CGAL_assertion(entry.second);
      entries_.erase(entry.first);
    }
  }

  // erase the base
  entries_bases_.erase(common_base);

  // erase the main one
  entries_.erase(pos);
}

template<typename MotorcycleGraphTraits_>
std::pair<typename Dictionary<MotorcycleGraphTraits_>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits_>::
find(const Dictionary_entry& e) const
{
  DEC_it res = entries().find(e);
  return std::make_pair(res, (res != entries().end()));
}

template<typename MotorcycleGraphTraits_>
std::pair<typename Dictionary<MotorcycleGraphTraits_>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits_>::
find(const Face_location& loc) const
{
  Dictionary_entry dummy(loc);
  return find(dummy);
}

template<typename MotorcycleGraphTraits_>
std::pair<typename Dictionary<MotorcycleGraphTraits_>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits_>::
get_sibling(const DEC_it e, const face_descriptor fd) const
{
  if(e->location().first == fd)
    return std::make_pair(e, true);

  // @todo: would be better if the siblings gave a 'DEC_it' immediately
  // instead of having to call a 'find' here.
  const Face_location& location = e->sibling(fd);
  return find(location);
}

template<typename MotorcycleGraphTraits_>
std::pair<typename Dictionary<MotorcycleGraphTraits_>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits_>::
insert(const Face_location& location, const Point& p, const Triangle_mesh& mesh)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  Dictionary_entry entry(location); // only the location is set up so far
  std::pair<DEC_it, bool> is_insert_successful = entries().insert(entry);

  if(!is_insert_successful.second)
  {
    std::cerr << "Warning: point already exists in the dictionary: "
              << &*(is_insert_successful.first) << std::endl
              << *(is_insert_successful.first) << std::endl;
  }
  else // Insertion is successful
  {
    // Initialize the base with all motorcycle information
    Dictionary_entry_base deb(p);
    entries_bases_.push_back(deb);
    DEBC_it common_base = boost::prior(entries_bases_.end());
    is_insert_successful.first->set_base(common_base);

    // If the point is on border: initialize and insert its siblings
    descriptor_variant dv = PMP::get_descriptor_from_location(location, mesh);
    face_descriptor initial_fd = location.first;

    CGAL_precondition(common_base->siblings().empty());

    if(dv.which() != 2) // point is not strictly on a face -> it's on the border
    {
      // fill for the initial face
      common_base->siblings().insert(location);

      if(dv.which() == 0) // on a vertex
      {
        vertex_descriptor vd = boost::get<vertex_descriptor>(dv);
        halfedge_descriptor hd = halfedge(vd, mesh);

        CGAL::Face_around_target_iterator<Triangle_mesh> fit, fend;
        boost::tie(fit, fend) = CGAL::faces_around_target(hd, mesh);
        std::size_t number_of_incident_faces = std::distance(fit, fend);
        for(; fit!=fend; ++fit)
        {
          face_descriptor fd = *fit;
          if(fd == initial_fd || fd == boost::graph_traits<Triangle_mesh>::null_face())
            continue;

          const Face_location location_in_fd = PMP::locate_in_adjacent_face(location, fd, mesh);

          // Set up the sibling information
          common_base->siblings().insert(location_in_fd);

          std::pair<DEC_it, bool> is_insert_successful_in_fd = entries().insert(location_in_fd);

          // Insertion must be successful: we insert all the possible location
          // when we insert a location on a border, so we can't have had successful
          // insertion for the main entry but not for equivalent locations.
          CGAL_assertion(is_insert_successful_in_fd.second);

          // Set the base
          is_insert_successful_in_fd.first->set_base(common_base);
        }
        CGAL_postcondition_code(if(PMP::is_on_mesh_border(location, mesh)))
        CGAL_postcondition(common_base->siblings().size() == number_of_incident_faces - 1 /*the null face*/);
        CGAL_postcondition_code(else)
        CGAL_postcondition(common_base->siblings().size() == number_of_incident_faces);
      }
      else if(dv.which() == 1) // on a halfedge
      {
        halfedge_descriptor hd = boost::get<halfedge_descriptor>(dv);
        if(!is_border(edge(hd, mesh), mesh))
        {
          face_descriptor fd = face(opposite(hd, mesh), mesh);
          const Face_location location_in_fd = PMP::locate_in_adjacent_face(location, fd, mesh);
          std::pair<DEC_it, bool> is_insert_successful_in_fd = entries().insert(location_in_fd);
          CGAL_assertion(is_insert_successful_in_fd.second);
          is_insert_successful_in_fd.first->set_base(common_base);
          common_base->siblings().insert(location_in_fd);
          CGAL_postcondition(common_base->siblings().size() == 2);
        }
      }
    }

    std::cout << "New point in the dictionary: "
              << &*(is_insert_successful.first) << std::endl
              << *(is_insert_successful.first) << std::endl;
  }

  return is_insert_successful;
}

template<typename MotorcycleGraphTraits_>
std::pair<typename Dictionary<MotorcycleGraphTraits_>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits_>::
insert(const Face_location& loc, const Point& p, const std::size_t i, const FT time, const Triangle_mesh& mesh)
{
  std::pair<DEC_it, bool> entry = insert(loc, p, mesh);
  entry.first->add_motorcycle(i, time);

return entry;
}

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary<MotorcycleGraphTraits>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const std::size_t i, const FT time, const Triangle_mesh& mesh)
{
  Point p = CGAL::Polygon_mesh_processing::location_to_point(loc, mesh);
  return insert(loc, p, i, time, mesh);
}

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary<MotorcycleGraphTraits>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const Triangle_mesh& mesh)
{
  Point p = CGAL::Polygon_mesh_processing::location_to_point(loc, mesh);
  return insert(loc, p, mesh);
}


} // namespace Polyline_tracing
} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_DICTIONARY_H
