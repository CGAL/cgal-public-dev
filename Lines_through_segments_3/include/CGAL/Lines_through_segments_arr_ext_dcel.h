// Copyright (c) 2010  Tel-Aviv University (Israel).
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
// $URL: $
// $Id: $
// 
//
// Author(s)     : Asaf Porat          <asafpor1@post.tau.ac.il>

#ifndef LINE_THROUGH_SEGMENTS_ARR_EXT_DCEL_H
#define LINE_THROUGH_SEGMENTS_ARR_EXT_DCEL_H

#include <CGAL/basic.h>
#include <CGAL/Arr_dcel_base.h>

/*************************************************************
 * Each vertex is extended with Boolean flag that determines whether it 
 * was added to output.
 * This flag is initialized to false. When an edge or a face that contains 
 * the vertex is added to the output list, this flag is set to true.
 *************************************************************/

namespace CGAL {

template <typename Point_2>
class Arrangement_vertex : public CGAL::Arr_vertex_base<Point_2> {
private:
      
  bool m_added_to_output;/* If true, this vertex is part of edge or face
                            that were already added to the output. */
public:
  /*! Constructor */
  Arrangement_vertex() {m_added_to_output = false;}

  /*! Set added to output flag. */
  void set_added_to_output(bool _added_to_output) 
  { 
    m_added_to_output = _added_to_output;
  }

  bool get_added_to_output() const
  {
    return m_added_to_output;
  }
      
  /*! Assign from another vertex.
   * \param v the other vertex
   */
  virtual void assign(Arrangement_vertex& v)
  {
    CGAL::Arr_vertex_base<Point_2>::assign(v);
  }
};

/*************************************************************
 * Each edge is extended with the following:
 *    1. Boolean flag that determines whether it was added to output.
 *       This flag is initialized to false. When a face that contains the
 *       edge is
 *       added to the output list, this flag is set to true.
 *    2. Ext_obj - pointer to the original rational segment that created the
 *       edge at the arrangement.
 *************************************************************/

template <typename X_monotone_curve_2,typename Ext_obj>
class Arrangement_segment_halfedge :
    public CGAL::Arr_halfedge_base<X_monotone_curve_2> {
private:
  /* A list of pointers to the original rational 
     segment that created the edge at the arrangement.*/
  std::list<const Ext_obj*> m_segment_list; 
  bool m_added_to_output;/* Each half edge is added twice to the output.
                            This Boolean indicates if this edge already
                            added to the output. */
      
public:
  typedef typename std::list<const Ext_obj* >::const_iterator
  const_iterator;
  const_iterator segs_begin() const 
  {
    return m_segment_list.begin(); 
  }
      
  const_iterator segs_end() const 
  {
    return m_segment_list.end(); 
  }

  /*! Constructor */
  Arrangement_segment_halfedge() { m_added_to_output = false;}
      
  /*! Obtain the segment */
  unsigned int num_of_segments() const { return m_segment_list.size(); }

  /*! Set the segment */
  void add_segment(const Ext_obj * _segment) 
  { 
    m_segment_list.push_back(_segment);
  }

  std::list<const Ext_obj*> get_segments_list()
  {
    return m_segment_list;
  }
      
  void set_segments_list(std::list<const Ext_obj*>& seg_list)
  {
    m_segment_list = seg_list;
  }
      
  /*! Assign from another halfedge.
   * \param he the other halfedge
   */
  virtual void assign(const Arrangement_segment_halfedge& he)
  {
    CGAL::Arr_halfedge_base<X_monotone_curve_2>::assign(he);
    m_segment_list = he.m_segment_list;
    m_added_to_output = he.m_added_to_output;
  }
      
  bool get_added_to_output()
  {
    return m_added_to_output;
  }

  void set_added_to_output(bool b)
  {
    m_added_to_output = b;
  }

  struct rebind
  {
    typedef Arrangement_segment_halfedge<X_monotone_curve_2,Ext_obj> other;
  };

};
  
/*************************************************************
 * Each edge is extended with:
 *    Counter to the number of plane faces.
 *    A plane face is a face that represents a plane that contains three 
 *    segments. 
 *    When two of these faces overlaps the new face represents a
 *    plane the contains 4 segments.
 *    This plane will be added to the output.
 *************************************************************/
template <typename Ext_obj>
class Arrangement_ext_face : public CGAL::Arr_face_base {
      
private:
  /* A list of pointers to the original rational 
     segment that created the face at the arrangement.*/
  std::list<const Ext_obj*> m_segment_list;

public:
  typedef typename std::list<const Ext_obj* >::const_iterator
  const_iterator;
  const_iterator segs_begin() const 
  {
    return m_segment_list.begin(); 
  }
      
  const_iterator segs_end() const 
  {
    return m_segment_list.end(); 
  }
  /*! Constructor */
  Arrangement_ext_face() 
  {
  }
      
  unsigned int num_of_overlap_plane_faces() const 
  { 
    return m_segment_list.size();
  }

  void clear()
  { 
    m_segment_list.clear();
  }

  void add_segment(const Ext_obj * _segment) 
  { 
    m_segment_list.push_back(_segment);
  }
      
  virtual void assign(const Arrangement_ext_face& f)
  {
    CGAL::Arr_face_base::assign(f);
    m_segment_list = f.m_segment_list;
  }
};

/*! A new dcel builder with extended features */
template <typename Traits,typename Ext_obj_>
class Lines_through_segments_arr_ext_dcel :
    public CGAL::Arr_dcel_base<Arrangement_vertex<typename Traits::Point_2>,
                               Arrangement_segment_halfedge<typename Traits::
                                                            X_monotone_curve_2,
                                                            Ext_obj_>,
                               Arrangement_ext_face<Ext_obj_> >
{
public:
  typedef Ext_obj_ Ext_obj;
  typedef typename Arrangement_segment_halfedge<
    typename Traits::X_monotone_curve_2,
    Ext_obj_>::const_iterator const_iterator;
  Lines_through_segments_arr_ext_dcel() {}
      
  template < typename Data_Traits>
  struct rebind
  {
    typedef Lines_through_segments_arr_ext_dcel<Data_Traits, Ext_obj> other;
  };
};

} //namespace CGAL

#endif
