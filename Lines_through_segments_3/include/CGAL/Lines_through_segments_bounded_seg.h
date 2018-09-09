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

#ifndef LINE_THROUGH_SEGMENTS_BOUNDED_SEG_H
#define LINE_THROUGH_SEGMENTS_BOUNDED_SEG_H

/*************************************************************
 * This class represents a bounded segment on the x axis.
 * The bounded segment is defined over rational numbers.
 *
 *************************************************************/

namespace CGAL {
  
template <typename NT>
class Lines_through_segments_rbound_unbound_union {

public:

private:
  NT m_bound;
  unsigned short m_is_bounded;

public:
   enum {
      LTS_BS_BOUNDED = 0,
      LTS_BS_UNBOUNDED_PLUS_INFINITY = 1,
      LTS_BS_UNBOUNDED_MINUS_INFINITY = 2,
      
   };

  Lines_through_segments_rbound_unbound_union()
  {
  }
      
      
  Lines_through_segments_rbound_unbound_union(const NT& bound,
                                              const unsigned short is_bounded)
  {
    m_bound = bound;
    m_is_bounded = is_bounded;
  }

  Lines_through_segments_rbound_unbound_union(const NT& bound)
  {
    m_bound = bound;
    m_is_bounded = LTS_BS_BOUNDED;
  }

  Lines_through_segments_rbound_unbound_union(const unsigned short is_bounded)
  {
    m_is_bounded = is_bounded;
  }

  Lines_through_segments_rbound_unbound_union(const Lines_through_segments_rbound_unbound_union& bound)
  {
    m_bound = bound.m_bound;
    m_is_bounded = bound.m_is_bounded;
  }
      
  bool is_bound()
  {
    return (m_is_bounded == LTS_BS_BOUNDED);
  }

  NT bound()
  {
    CGAL_assertion(m_is_bounded == LTS_BS_BOUNDED);
    return m_bound;
  }


  bool operator==(const Lines_through_segments_rbound_unbound_union& bound)
  {
    return ((m_is_bounded == LTS_BS_BOUNDED && 
             bound.m_is_bounded == LTS_BS_BOUNDED) && 
            (this->m_bound == bound.m_bound));
  }

  bool operator!=(const Lines_through_segments_rbound_unbound_union& bound)
  {
    return (!(*this == bound));
  }

  bool operator==(const NT& bound)
  {
    return (m_is_bounded == LTS_BS_BOUNDED && (this->m_bound == bound));
  }

  bool operator>=(const NT& bound)
  {
    return (*this >= Lines_through_segments_rbound_unbound_union(bound));
  }
      
  bool operator<=(const NT& bound)
  {
    return (*this <= Lines_through_segments_rbound_unbound_union(bound));
  }

  bool operator==(const unsigned short is_bound)
  {
    return (m_is_bounded == is_bound);
  }

  bool operator<(const Lines_through_segments_rbound_unbound_union& bound)
  {
    if ((m_is_bounded == LTS_BS_BOUNDED && 
         bound.m_is_bounded == LTS_BS_BOUNDED && 
         this->m_bound < bound.m_bound))
      return true;
         
    if (m_is_bounded == LTS_BS_UNBOUNDED_MINUS_INFINITY && 
        bound.m_is_bounded == LTS_BS_BOUNDED)
      return true;

    if (m_is_bounded == LTS_BS_BOUNDED && 
        bound.m_is_bounded == LTS_BS_UNBOUNDED_PLUS_INFINITY)
      return true;

    if (m_is_bounded == LTS_BS_UNBOUNDED_MINUS_INFINITY && 
        bound.m_is_bounded == LTS_BS_UNBOUNDED_PLUS_INFINITY)
      return true;

    return false;
  }

  bool operator>(const Lines_through_segments_rbound_unbound_union& bound)
  {
    return (!(*this <= bound));
  }

  bool operator<=(const Lines_through_segments_rbound_unbound_union& bound)
  {
    return ((*this < bound) || (*this == bound));
  }

  bool operator>=(const Lines_through_segments_rbound_unbound_union& bound)
  {
    return (!(*this < bound));
  }

  void operator=(const NT& bound)
  {
    m_is_bounded = LTS_BS_BOUNDED;
    m_bound = bound;
  }

  void operator=(const unsigned short is_bounded)
  {
    m_is_bounded = is_bounded;
  }
      
  std::string to_string() const
  {
    std::ostringstream o;
    if (m_is_bounded == LTS_BS_BOUNDED)
      o << m_bound;
    else if (m_is_bounded == LTS_BS_UNBOUNDED_PLUS_INFINITY)
      o << "LTS_BS_UNBOUNDED_PLUS_INFINITY";
    else
      o << "LTS_BS_UNBOUNDED_MINUS_INFINITY";
    return o.str();
  }
};
   
template <typename Rational>
class Lines_through_segments_bounded_seg {

private:      
  typedef Lines_through_segments_rbound_unbound_union<Rational>  LTS_rbound;
   
  LTS_rbound m_max_x; /* Holds the max value of x.*/
  bool m_segment_include_max_x; /* true if the segment contains max_x. */

  /* Holds the max value of x.*/
  LTS_rbound m_min_x; 
  bool m_segment_include_min_x; /* true if the segment contains min_x. */

  bool m_valid; /* True if the segment is valid, i.e, if min_x < max_x*/
public:
  Lines_through_segments_bounded_seg() 
  {
    m_valid = false;
    m_min_x = LTS_rbound::LTS_BS_UNBOUNDED_PLUS_INFINITY;
    m_max_x = LTS_rbound::LTS_BS_UNBOUNDED_MINUS_INFINITY;
    m_segment_include_min_x = false;
    m_segment_include_max_x = false;
  }

  Lines_through_segments_bounded_seg(
     const LTS_rbound& _max_x, 
     const LTS_rbound& _min_x,
     bool _segment_include_max_x,
     bool _segment_include_min_x)
  { 
    m_max_x = _max_x;
    m_min_x = _min_x;

    m_segment_include_min_x = _segment_include_min_x;
    m_segment_include_max_x = _segment_include_max_x;
         
    m_valid = (m_min_x < m_max_x || 
               (m_min_x == m_max_x && m_segment_include_min_x && 
                m_segment_include_max_x));
  }

  Lines_through_segments_bounded_seg(const Lines_through_segments_bounded_seg&
                                     _segment)
  { 
    m_max_x = _segment.m_max_x;
    m_min_x = _segment.m_min_x;

    m_segment_include_min_x = _segment.m_segment_include_min_x;
    m_segment_include_max_x = _segment.m_segment_include_max_x;

    m_valid = _segment.m_valid;
  }
      
  bool is_valid()
  {
    return m_valid;
  }

  LTS_rbound get_min()
  { 
    return m_min_x; 
  }
  LTS_rbound get_max()
  {
    return m_max_x; 
  }

  void set_min(const Rational& _min_x,
               bool _segment_include_min_x)
  { 
    m_min_x = _min_x;
        
    m_segment_include_min_x = _segment_include_min_x;
    m_valid = (m_min_x < m_max_x ||
               (m_min_x == m_max_x && m_segment_include_min_x &&
                m_segment_include_max_x));
  }
  void set_max(const Rational& _max_x,
               bool _segment_include_max_x)
  {
    m_max_x = _max_x; 
    m_segment_include_max_x = _segment_include_max_x;
    m_valid = (m_min_x < m_max_x ||
               (m_min_x == m_max_x && m_segment_include_min_x &&
                m_segment_include_max_x));
  }

  bool segment_contain_max(){ return this->m_segment_include_max_x; }
  bool segment_contain_min(){ return this->m_segment_include_min_x; }

  std::string to_string() const
  {
    std::ostringstream o;
        
    if (m_segment_include_max_x && m_segment_include_min_x)
      o << "Bounded segment  = [" << m_min_x <<"," << m_max_x
        << "] valid = " << m_valid << std::endl;
    else if (m_segment_include_max_x)
      o << "Bounded segment  = (" << m_min_x <<"," << m_max_x 
        << "] valid = " << m_valid << std::endl;
    else if (m_segment_include_min_x)
      o << "Bounded segment  = [" << m_min_x <<"," << m_max_x 
        << ") valid = " << m_valid << std::endl;
    else
      o << "Bounded segment  = (" << m_min_x <<"," << m_max_x 
        << ") valid = " << m_valid << std::endl;
    return o.str();
  }

  ~Lines_through_segments_bounded_seg(){}
};

} //namespace CGAL

#endif //LINE_THROUGH_SEGMENTS_BOUNDED_SEG_H
