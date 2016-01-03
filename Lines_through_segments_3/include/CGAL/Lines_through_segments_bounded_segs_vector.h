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

#ifndef LINE_THROUGH_SEGMENTS_BOUNDED_SEGS_VECTOR_H
#define LINE_THROUGH_SEGMENTS_BOUNDED_SEGS_VECTOR_H

/*************************************************************
 * This class represent a vector of 2 dimensional Bounded segments at the
 * segment [0,1].
 * The bounded segments are sorted over the segment [0,1] from left to right.
 * The bounded segments are defined with rational numbers.
 *
 *************************************************************/

#include <list>
#include <string>

#include <CGAL/basic.h>
#include <CGAL/Lines_through_segments_bounded_seg.h>

namespace CGAL {

template <typename Rational>
class Lines_through_segments_bounded_segs_vector {
  typedef Lines_through_segments_bounded_seg<Rational> Bounded_seg;
  typedef Lines_through_segments_rbound_unbound_union<Rational> Rbound;

private:
  std::list<Bounded_seg> m_bounded_segments;

public:
  typedef typename std::list<Bounded_seg>::iterator       iterator;
  typedef typename std::list<Bounded_seg>::const_iterator const_iterator;
      
  // iterator support
  iterator begin() { return m_bounded_segments.begin(); }
  const_iterator begin() const { return m_bounded_segments.begin(); }
  iterator end() { return m_bounded_segments.end(); }
  const_iterator end() const { return m_bounded_segments.end(); }
  typedef Bounded_seg&       reference;
      
public:
      
  Lines_through_segments_bounded_segs_vector() { }

  Lines_through_segments_bounded_segs_vector(const Lines_through_segments_bounded_segs_vector& segments)
  { 
    const_iterator it;
    for (it = segments.begin(); it != segments.end(); ++it)
    {
      add_bounded_seg(*it);
    }
  }
      
  unsigned int size(){ return m_bounded_segments.size(); }

  void add_bounded_seg(Bounded_seg segment)
  {
#if BOUNDED_SEGMENTS_VECTOR_DEBUG
    std::cout << "Add bounded segment " << segment << std::endl;
#endif

    if (segment.is_valid())
    {
      if (m_bounded_segments.size() == 0)
      {
        m_bounded_segments.push_front(segment);
      }
      else
      {
        iterator it;

        /* Enter the segments sorted from left to right */
        for (it = m_bounded_segments.begin();
             it != m_bounded_segments.end(); ++it)
        {
          if (segment.get_max() <= (*it).get_min())
          {
            m_bounded_segments.insert(it, segment);
            return;
          }
        }
        m_bounded_segments.insert(it, segment);
      }
    }
  }

  /*************************************************************
   * Function description:
   * -------------------- 
   * The following constructor gets 2 vectors of bounded segments, 
   * and union them into one vector.
   *
   * Input:
   *
   *  first_bounded_segments - vector of bounded segments.
   *  second_bounded_segments - vector of bounded segments.
   *
   * Output:
   * 
   *
   *************************************************************/

  void merge_bounded_segs_vectors(Lines_through_segments_bounded_segs_vector& first_bounded_segments,
                                  Lines_through_segments_bounded_segs_vector& second_bounded_segments)
  {
#if BOUNDED_SEGMENTS_VECTOR_DEBUG         
    std::cout << "first_bounded_segments = " << std::endl;
    std::cout << first_bounded_segments << std::endl;
    std::cout << "second_bounded_segments = " << std::endl;
    std::cout << second_bounded_segments << std::endl;
#endif
    iterator first_it;
    iterator second_it;

    /* Clear all old bounded segments. */
    m_bounded_segments.clear();
    
    for (first_it = first_bounded_segments.begin();
         first_it != first_bounded_segments.end(); ++first_it)
    {
      for (second_it = second_bounded_segments.begin();
           second_it != second_bounded_segments.end(); ++second_it)
      {
        Rbound min_x;
        Rbound max_x;
        bool segment_contain_min = false;
        bool segment_contain_max = false;
               
        if ((*first_it).get_max() < (*second_it).get_max())
        {
          max_x = (*first_it).get_max();
          segment_contain_max = (*first_it).segment_contain_max();
        }
        else if ((*first_it).get_max() > (*second_it).get_max())
        {
          max_x = (*second_it).get_max();
          segment_contain_max = (*second_it).segment_contain_max();
        }
        else if ((*first_it).get_max() == (*second_it).get_max())
        {
          max_x = (*second_it).get_max();
          segment_contain_max = (*second_it).segment_contain_max() &&
            (*first_it).segment_contain_max();
        }
               
        if ((*first_it).get_min() > (*second_it).get_min())
        {
          min_x = (*first_it).get_min();
          segment_contain_min = (*first_it).segment_contain_min();
        }
        else if ((*first_it).get_min() < (*second_it).get_min())
        {
          min_x = (*second_it).get_min();
          segment_contain_min = (*second_it).segment_contain_min();
        }
        else if ((*first_it).get_min() == (*second_it).get_min())
        {
          min_x = (*second_it).get_min();
          segment_contain_min = (*second_it).segment_contain_min() &&
            (*first_it).segment_contain_min();
        }
               
        this->add_bounded_seg(Bounded_seg(max_x,min_x, segment_contain_max,
                                          segment_contain_min));
      }
    }
#if BOUNDED_SEGMENTS_VECTOR_DEBUG   
    std::cout << "res bounded_segments = " << std::endl;
    std::cout << *this << std::endl;
#endif
         
  }

  void clear()
  {
    this->m_bounded_segments.clear();
  }
      
  std::string to_string()
  {
    iterator it;
    std::ostringstream o;
    o << "size = " << this->m_bounded_segments.size() << std::endl;
    for (it = this->m_bounded_segments.begin();
         it != this->m_bounded_segments.end(); ++it)
    {
      o << (*it) << std::endl;
    }
         
    return o.str();
  }

  ~Lines_through_segments_bounded_segs_vector()
  {
    m_bounded_segments.clear();
  }
};

} //namespace CGAL

#endif //LINE_THROUGH_SEGMENTS_BOUNDED_SEGS_VECTOR_H
