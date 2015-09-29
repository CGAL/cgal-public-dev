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


#ifndef LINE_THROUGH_SEGMENTS_FIND_OVERLAP_LINES_H
#define LINE_THROUGH_SEGMENTS_FIND_OVERLAP_LINES_H

#include <CGAL/Lines_through_segments_3/internal.h>

#include <utility>
#include <vector>

#include <boost/array.hpp>

namespace CGAL { namespace LTS { namespace internal {

struct Compare_points_on_line {
  /* Find the maximal common segment, I.e. the segment the the maximal number of
     lines include. */
  template <typename T>
  bool operator() (const std::pair<bool, T>& x,
                   const std::pair<bool, T>& y) const
  { return this->operator()(x.second, y.second); }

  template <typename T>
  bool operator()(const T& x, const T& y) const {
    if (x.x() < y.x())
      return true;
    else if (x.x() == y.x() && 
             x.y() < y.y())
      return true;
    else if (x.x() == y.x() && 
             x.y() == y.y() && 
             x.z() <= y.z())
      return true;
    return false;
  }
};

} // internal

/// XXX try to get rid of the traits argument and externalize the
/// handling of overlap into a functor
template <typename Traits,
          typename ForwardIterator,
          typename OutputIterator,
          typename With_segments>
OutputIterator find_overlap(Traits /*t*/, ForwardIterator begin, ForwardIterator end, 
                            OutputIterator out_iter, With_segments with_segments) {
  typedef typename Traits::Rational_kernel::Point_3 Rational_point_3;
  typedef typename Traits::Rational_kernel::Segment_3 Rational_segment_3;
  
  typedef std::pair<bool, Rational_point_3> List_element;
  std::vector<List_element> points_vector;
    
  ForwardIterator save = begin;
  internal::Compare_points_on_line comp;
  while(begin != end) {
    if (comp((*begin)->source(), (*begin)->target()))
    {
      points_vector.push_back(std::make_pair(true,(*begin)->source()));
      points_vector.push_back(std::make_pair(false,(*begin)->target()));
    }
    else
    {
      points_vector.push_back(std::make_pair(true,(*begin)->target()));
      points_vector.push_back(std::make_pair(false,(*begin)->source()));
    }
    ++begin;
  }
    
  std::sort(points_vector.begin(), points_vector.end(), comp);
         
  int num_of_overlap_lines = 0;
  bool started = false;
  Rational_point_3 start_point;
  for (typename std::vector<List_element>::iterator
         it = points_vector.begin();
       it != points_vector.end(); ++it)
  {
    if (it->first)
    {
      num_of_overlap_lines++;
    }
    else
    {
      num_of_overlap_lines--;
    }
         
    if (!started && num_of_overlap_lines == 4)
    {
      started = true;
      start_point = it->second;
    }
    else if (started && num_of_overlap_lines == 3)
    {
      Rational_segment_3 output_seg(start_point, it->second);
      out_iter = insert_transversal(out_iter, output_seg, 
                                    *save, *(++save), *(++save), 
                                    *(++save), with_segments);
    }
  }
  return out_iter;
}

} // LTS
} // CGAL

#endif /* LINE_THROUGH_SEGMENTS_FIND_OVERLAP_LINES_H */
