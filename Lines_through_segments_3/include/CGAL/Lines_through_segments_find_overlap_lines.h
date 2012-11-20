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

#include <CGAL/Lines_through_segments_output_obj.h>

namespace CGAL {

template <typename Lines_through_segments_traits_3,
          typename With_segments>
class Lines_through_segments_find_overlap_lines {

public:
  typedef Lines_through_segments_traits_3       Traits_3;
private:
  typedef typename Traits_3::Rational_kernel    Rational_kernel;
  typedef typename Traits_3::Alg_kernel         Alg_kernel;
  typedef typename Rational_kernel::Point_3     Rational_point_3;
  typedef typename Rational_kernel::Segment_3   Rational_segment_3;
   
private:
  class List_element {
  public:
      
    bool m_is_start; /* TRUE start, FALSE end */
    Rational_point_3 m_point;
      
    List_element(bool is_start, /* TRUE start, FALSE end */
                 Rational_point_3 point)
    {
      m_is_start = is_start; /* TRUE start, FALSE end */
      m_point = point;
    }

    List_element(Rational_point_3 point)
    {
      m_is_start = false; /* TRUE start, FALSE end */
      m_point = point;
    }

    List_element(const List_element &list_el)
    {
      m_is_start = list_el.m_is_start; /* TRUE start, FALSE end */
      m_point = list_el.m_point;
    }

  };
   
  class Compare_points_on_line {
  public:
    bool operator() (const List_element& first,
                     const List_element& second)
    {
      if (first.m_point.x() < second.m_point.x())
        return true;
      else if (first.m_point.x() == second.m_point.x() && 
               first.m_point.y() < second.m_point.y())
        return true;
      else if (first.m_point.x() == second.m_point.x() && 
               first.m_point.y() == second.m_point.y() && 
               first.m_point.z() <= second.m_point.z())
        return true;
      return false;
    }
  };

public:
  /* Find the maximal common segment, I.e. the segment the the maximal number of
     lines include. */
  template <typename Insert_iterator, typename Segs_vector>
  void operator()(unsigned int S1,
                  Segs_vector& segs,
                  Insert_iterator* insert_it)
  {
    Compare_points_on_line comp;
    std::list<List_element> points_list;
    unsigned int saved_S1 = S1;
    for (; S1 < segs.size(); ++S1)
    {
      if (comp(List_element(segs[S1]->source()),
               List_element(segs[S1]->target())))
      {
        points_list.push_back(List_element(true,segs[S1]->source()));
        points_list.push_back(List_element(false,segs[S1]->target()));
      }
      else
      {
        points_list.push_back(List_element(true,segs[S1]->target()));
        points_list.push_back(List_element(false,segs[S1]->source()));
      }
    }
    points_list.sort(comp);
         
    typename std::list<List_element>::iterator it;
    int num_of_overlap_lines = 0;
    bool started = false;
    Rational_point_3 start_point;
    for (it = points_list.begin(); it != points_list.end(); ++it)
    {
      if (it->m_is_start)
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
        start_point = it->m_point;
      }
      else if (started && num_of_overlap_lines == 3)
      {
        typedef boost::false_type With_arrangement;

        Rational_segment_3 output_seg(start_point, it->m_point);

        Lines_through_segments_arr_gen_func<
          Traits_3, With_segments, With_arrangement
          > m_arr_g_func;

        started = false;
        
        m_arr_g_func.insert_transversal_to_output(insert_it,
                                                  output_seg,
                                                  With_segments(),
                                                  segs[saved_S1],
                                                  segs[saved_S1 + 1],
                                                  segs[saved_S1 + 2],
                                                  segs[saved_S1 + 3]);
      }
    }
  }
};

} //namespace CGAL

#endif /* LINE_THROUGH_SEGMENTS_FIND_OVERLAP_LINES_H */
