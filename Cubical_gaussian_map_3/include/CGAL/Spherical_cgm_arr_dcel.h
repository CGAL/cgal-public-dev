// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Kapelushnik Lior <liorkape@post.tau.ac.il>

/*! \file
 * spherical arrangements of none intersecting arcs of great circles on a sphere
 */

#ifndef CGAL_SPHERICAL_CGM_ARR_DCEL_H
#define CGAL_SPHERICAL_CGM_ARR_DCEL_H

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>

#include <list>

#include <CGAL/Arrangement_2.h>
#include "CGAL/Cgm_arr_dcel.h"

namespace CGAL {

/*
 Extend the planar-map vertex to hold information used by the spherical map

 Point_2 - a type for the cubical map vertex point
 Direction_3 - a direction type for the spherical map
*/
template <class T_Point_2, 
          class T_Direction_3 = Cartesian<Gmpq>::Direction_3 >
class Spherical_cgm_arr_vertex : public Cgm_arr_vertex<T_Point_2> { 

public:
  /*
   empty constructor
  */
  Spherical_cgm_arr_vertex(); 
  /*
   set the direction represented by a cubical map vertex

   direction - the spherical direction that the cubical vertex represents
  */
  inline void set_direction(const T_Direction_3 &direction);
  /*
   get the direction which a cubical vertex is representing

   return value - the direction represented by the vertex
  */
  inline const T_Direction_3 &get_direction() const;
  /*
   set as a representative vertex, in case the same direction has
   more than one cubical vertex, only one vertex will be representative
  */
  inline void set_rep();
  /*
   set as an arc vertex, the vertex is generated by an arc endpoint
  */
  inline void setReal();
  /*
   find if a vertex is generated by an arc endpoint direction

   return value - true if an arc vertex (represents an endpoint of an arc)
  */
  inline bool getReal() const;
  /*
   find if vertex is unique representative

   return value - true if vertex is a direction representative
  */
  inline bool is_rep() const;
  /*
   get the sphere vertex pointing to the cgm vertex

   return value - a pointer to Vertex_handle as defined in the Spherical_map
     which points to the cubical vertex
  */
  inline void *getSphereVertex() const {
    return m_sphereVertex;
  }
  /*
   set the sphere vertex pointing to the cgm vertex

   itPtr - a pointer to Vertex_handle as defined in the Spherical_map
     which points to the cubical vertex (will be created in spherical map update)
  */
  inline void setSphereVertex(void *itPtr) {
    m_sphereVertex = itPtr;
  }
  /*
   copy extended features from another vertex

   v - the vertex to copy extended features from
  */
  virtual void assign(const Spherical_cgm_arr_vertex & v); 
private:
  T_Direction_3 m_dir; // the direction of a vertex
  bool   m_isDist; // is vertex distinct representative in degtenerate case
  bool  m_isReal; // is the vertex part of an arc ending points

  // pointer to the sphere vertex handle that points to this cgm vertex
  void *m_sphereVertex;
};

/*
 empty constructor, initializes flags and data
*/
template <class T_Point_2, class T_Direction_3>
  Spherical_cgm_arr_vertex<T_Point_2, T_Direction_3>::
Spherical_cgm_arr_vertex():
  Cgm_arr_vertex<T_Point_2>(), m_dir(T_Direction_3(0,0,0)), m_isDist(false),
  m_isReal(false), m_sphereVertex(0) {};

/*
 set the direction represented by a cubical map vertex

 direction - the spherical direction that the cubical vertex represents
*/
template <class T_Point_2, class T_Direction_3>
inline void Spherical_cgm_arr_vertex<T_Point_2, T_Direction_3>::
set_direction(const T_Direction_3 &direction) {
  m_dir = direction;
}

/*
 get the direction which a cubical vertex is representing

 return value - the direction represented by the vertex
*/
template <class T_Point_2, class T_Direction_3>
inline const T_Direction_3 &Spherical_cgm_arr_vertex<T_Point_2, T_Direction_3>::
get_direction() const {
  return m_dir;
}

/*
 set as a representative vertex, in case the same direction has
 more than one cubical vertex, only one vertex will be representative
*/
template <class T_Point_2, class T_Direction_3>
inline void Spherical_cgm_arr_vertex<T_Point_2, T_Direction_3>::
set_rep() {
  m_isDist = true;
}

/*
 set as an arc vertex, the vertex is generated by an arc endpoint
*/
template <class T_Point_2, class T_Direction_3>
inline void Spherical_cgm_arr_vertex<T_Point_2, T_Direction_3>::
setReal() {
  m_isReal = true;
}

/*
 find if a vertex is generated by an arc endpoint direction

 return value - true if an arc vertex (represents an endpoint of an arc)
*/
template <class T_Point_2, class T_Direction_3>
inline bool Spherical_cgm_arr_vertex<T_Point_2, T_Direction_3>::
getReal() const {
  return m_isReal;
}

/*
 find if vertex is unique representative

 return value - true if vertex is a direction representative
*/
template <class T_Point_2, class T_Direction_3>
inline bool Spherical_cgm_arr_vertex<T_Point_2, T_Direction_3>::
is_rep() const {
  return m_isDist;
}

/*
 copy extended features from another vertex

 v - the vertex to copy extended features from
*/
template <class T_Point_2, class T_Direction_3>
  void Spherical_cgm_arr_vertex<T_Point_2, T_Direction_3>::
assign(const Spherical_cgm_arr_vertex & v)
{
  Cgm_arr_vertex<T_Point_2>::assign(v);
  m_dir = v.m_dir;
  m_isDist = v.m_isDist;
  m_isReal = v.m_isReal;
  set_is_real(v.is_real());
}

/*
 Extend the planar-map halfedge to hold information used by the spherical map

 X_monotone_curve_2 - the type of curve used by the cubical gaussian map
*/
template <class X_monotone_curve_2>
class Spherical_cgm_arr_halfedge : public Cgm_arr_halfedge<X_monotone_curve_2> 
{
public:
  /*
   empty constructor
  */
  Spherical_cgm_arr_halfedge();
  /*
   set mark state of the halfedge

   val - the new halfedge mark state
  */
  inline void setMark(bool val=true);
  /*
   get current mark status

   return value - the mark status
  */
  inline bool getMark();
  /*
   get the sphere halfedge pointing to the cgm halfedge

   return value - a pointer to Halfedge_handle as defined in the Spherical_map
     which points to the cubical halfedge
  */
  inline void *getSphereEdge() const {
      return m_sphereHalfedge;
  }

  /*
   set the sphere halfedge pointing to the cgm halfedge

   itPtr - a pointer to Halfedge_handle as defined in the Spherical_map
     which points to the cubical halfedge (will be created in spherical map update)
  */
  inline void setSphereEdge(void *itPtr) {
    m_sphereHalfedge = itPtr;
  }

  /*
   copy extended features from another halfedge

   e - the halfedge to copy extended features from
  */
  virtual void assign(const Spherical_cgm_arr_halfedge & e);
private:
  bool m_isMarked; // is the halfedge currently marked
  // pointer to the sphere halfedge that points to this cgm halfedge
  void *m_sphereHalfedge;
};

/*
 empty constructor, initializes flags and data
*/
template <class X_monotone_curve_2>
Spherical_cgm_arr_halfedge<X_monotone_curve_2>::
Spherical_cgm_arr_halfedge():
  Cgm_arr_halfedge<X_monotone_curve_2>(), m_isMarked(false),
  m_sphereHalfedge(0) {};

/*
 set mark state of the halfedge

 val - the new halfedge mark state
*/
template <class X_monotone_curve_2>
inline void Spherical_cgm_arr_halfedge<X_monotone_curve_2>::
setMark(bool val) {
  m_isMarked = val;
}

/*
 get current mark status

 return value - the mark status
*/
template <class X_monotone_curve_2>
inline bool Spherical_cgm_arr_halfedge<X_monotone_curve_2>::
getMark() {
  return m_isMarked;
}

/*
 copy extended features from another halfedge

 e - the halfedge to copy extended features from
*/
template <class X_monotone_curve_2>
void Spherical_cgm_arr_halfedge<X_monotone_curve_2>::
assign(const Spherical_cgm_arr_halfedge & e) {
  Cgm_arr_halfedge<X_monotone_curve_2>::assign(e);
  m_isMarked = e.m_isMarked;
}

/*
 Extend the planar-map halfedge to hold information used by the spherical map
*/
class Spherical_cgm_arr_face: public Cgm_arr_face
{
public:
  /*
   empty constructor
  */
  Spherical_cgm_arr_face(): 
    Cgm_arr_face(), m_isMarked(false), m_sphereFace(0) {};

  /*
   set the mark state of the face
    
     val - the new mark state
  */
  inline void setMark(bool val=true);

  /*
   get the mark value
    
   return value - the face mark state
  */
  inline bool getMark();
  
  /*
   get the sphere face pointing to the cgm face

   return value - a pointer to Face_handle as defined in the Spherical_map
     which points to the cubical face
  */
  inline void *getSphereFace() const {
      return m_sphereFace;
  }

  /*
   set the sphere face pointing to the cgm face

   itPtr - a pointer to Face_handle as defined in the Spherical_map
     which points to the cubical face (will be created in spherical map update)
  */
  inline void setSphereFace(void *itPtr) {
    m_sphereFace = itPtr;
  }

  /*
   copy extended features from another face

   f - the face to copy extended features from
  */
  virtual void assign(const Spherical_cgm_arr_face & f) {
    Cgm_arr_face::assign(f);
    m_isMarked = f.m_isMarked;
  }

 private:
  bool m_isMarked; // is the face currently marked
  void *m_sphereFace; // the sphere face containing the cgm face
};

//Spherical_cgm_arr_face::Spherical_cgm_arr_face():
//  m_isMarked(false), m_sphereFace(0) {};

//  Cgm_arr_face(), m_isMarked(false), m_sphereFace(0) {};

/*
 set the mark state of the face
    
 val - the new mark state
*/
inline void Spherical_cgm_arr_face::
setMark(bool val) {
  m_isMarked = val;
}

/*
get the mark value

return value - the face mark state
*/
inline bool Spherical_cgm_arr_face::
getMark() {
  return m_isMarked;
}

/*
 a cubical arrangement dcel that contains extended features to be used
 by the spherical arrangement

 Traits - a traits class for the cubical gaussian map
   and for the direction type of spherical arcs endpoints
*/
template <class Traits>
class Spherical_cgm_arr_dcel: public Arr_dcel_base<
  Spherical_cgm_arr_vertex<typename Traits::Point_2,
                           typename Traits::Kernel::Direction_3>,
  Spherical_cgm_arr_halfedge<typename Traits::X_monotone_curve_2>,
  Spherical_cgm_arr_face>
{
public:
  /*
   empty constructor
  */
  Spherical_cgm_arr_dcel();
};

/*
 empty constructor
*/
template <class Traits>
Spherical_cgm_arr_dcel<Traits>::Spherical_cgm_arr_dcel() {};

} //namespace CGAL

#endif
