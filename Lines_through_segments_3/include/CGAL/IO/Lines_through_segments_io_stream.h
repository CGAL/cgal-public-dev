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

#include <sstream>
#include <CGAL/Lines_through_segments_3.h>
#ifndef LINE_THROUGH_SEGMENTS_IO_STREAM_H
#define LINE_THROUGH_SEGMENTS_IO_STREAM_H

namespace CGAL {

template <typename Lines_through_segments_traits_3, typename With_segments>
inline std::ostream&
operator<<(std::ostream& out,
           Lines_through_segments_arr_object<Lines_through_segments_traits_3,
                                             With_segments>& to_print)
{
  out << to_print.to_string();
  return out;
}


template <typename Rational>
inline std::ostream&
operator<<(std::ostream& out,
           const Lines_through_segments_bounded_seg<Rational>& to_print)
{
  out << to_print.to_string();
  return out;
}

template <typename Rational>
inline std::ostream&
operator<<(std::ostream& out,
           const Lines_through_segments_rbound_unbound_union<Rational>&
           to_print)
{
  out << to_print.to_string();
  return out;
}

template <typename Rational>
inline std::ostream&
operator<<(std::ostream & out,
           Lines_through_segments_bounded_segs_vector<Rational>& to_print)
{
  out << to_print.to_string();
  return out;
}

template <typename Lines_through_segments_traits_3, typename Insert_iterator,
          typename With_segments>
inline std::ostream&
operator<<(std::ostream& out,
           Lines_through_segments_impl<Lines_through_segments_traits_3,
           With_segments>& to_print)
{
  out << to_print.to_string();
  return out;
}

template <typename Point_and_segment_pair, typename Compare>
inline std::ostream&
operator<<(std::ostream& out,
           Lines_through_segments_isolated_points<Point_and_segment_pair,
                                                  Compare>& to_print)
{
  out << to_print.to_string();
  return out;
}

template <typename Traits_3, typename Point, typename Number_type>
inline std::ostream&
operator<<(std::ostream & out,
           const Lines_through_segments_point_adapt_3<Traits_3, Point,
                                                      Number_type>& to_print)
{
  out << to_print.to_string();
  return out;
}
   
template <typename Traits_3, typename Point_, typename Number_type>
inline std::ostream&
operator<<(std::ostream& out,
           const Lines_through_segments_point_adapt_2<Traits_3, Point_,
                                                      Number_type>& to_print)
{
  out << to_print.to_string();
  return out;
}

template <typename Traits_3>
inline std::ostream&
operator<<(std::ostream & out,
           const Lines_through_segments_through_3<Traits_3>& to_print)
{
  out << to_print.to_string();
  return out;
}

template <typename Traits_3>
inline std::ostream&
operator<<(std::ostream& out,
           const Lines_through_segments_mapped_2<Traits_3>& to_print)
{
  out << to_print.to_string();
  return out;
}

template <typename Transversal, typename Segments>
inline std::ostream&
operator<<(std::ostream& out,
           const std::pair<Transversal, Segments>& to_print)
{
  out << "S1 = " << *to_print.second[0] << std::endl;
  out << "S2 = " << *to_print.second[1] << std::endl;
  out << "S3 = " << *to_print.second[2] << std::endl;
  out << "S4 = " << *to_print.second[3] << std::endl;
  out << to_print.first << std::endl;
      
  return out;
}

} //namespace CGAL

#endif //LINE_THROUGH_SEGMENTS_IO_STREAM_H
