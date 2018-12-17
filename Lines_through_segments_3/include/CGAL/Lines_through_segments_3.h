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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s): Asaf Porat        <asafpor1@post.tau.ac.il>
//            Efi Fogel         <efif@post.tau.ac.il>

#ifndef LINE_THROUGH_SEGMENTS_3_H
#define LINE_THROUGH_SEGMENTS_3_H

/*! \file
*************************************************************
* The following class computes all of the lines that crosses 4 segments in 3D.
* The class is a functor.
*
* Input:
*
* Lines_through_segments_traits_3 - dedicated traits class.
*
*************************************************************
*/
#include <CGAL/Lines_through_segments_impl.h>
#include <CGAL/Lines_through_segments_3/internal.h>
#include <CGAL/Lines_through_segments_3/find_overlap.h>
#include <CGAL/Lines_through_segments_3/exceptions.h>
#include <CGAL/Lines_through_segments_output_obj.h>

#include <utility>
#include <memory>

#include <boost/shared_ptr.hpp>
#include <boost/variant.hpp>

namespace CGAL {

template <typename Traits_3_,
          typename With_segments_ = boost::false_type>
class Lines_through_segments_3 {

public:
  typedef Traits_3_                                     Traits_3;
  typedef typename Traits_3::Rational_kernel            Rational_kernel;
  typedef typename Rational_kernel::Segment_3           Segment_3;
  typedef typename Traits_3::Alg_kernel                 Alg_kernel;
  typedef With_segments_                                With_segments;

  typedef CGAL::Arrangement_2<
      typename Traits_3::Traits_arr_on_plane_2
    , Lines_through_segments_arr_ext_dcel<
          typename Traits_3::Traits_arr_on_plane_2
        , Segment_3>
    > Arr_on_plane;

  typedef CGAL::Arrangement_on_surface_2<
      typename Traits_3::Traits_arr_on_sphere_2
    , CGAL::Arr_spherical_topology_traits_2<
          typename Traits_3::Traits_arr_on_sphere_2
        , Lines_through_segments_arr_ext_dcel<
              typename Traits_3::Traits_arr_on_sphere_2
            , Segment_3> >
    > Arr_on_sphere;

  /* Output objects */
  typedef typename Rational_kernel::Line_3              Line_3;

  typedef Lines_through_segments_through_3<Traits_3>    Through_3;
  typedef typename Lines_through_segments_through_3<Traits_3>::Point_3
                                                        Through_point_3;
  typedef typename Lines_through_segments_through_3<Traits_3>::Segment_3
                                                        Through_segment_3;
  typedef typename Lines_through_segments_through_3<Traits_3>::Point_3_segment_3
    Through_point_3_segment_3;

  typedef Lines_through_segments_mapped_2<Traits_3>     Mapped_2;
  typedef typename Lines_through_segments_mapped_2<Traits_3>::X_monotone_curve_2
    Mapped_x_monotone_curve_2;
  typedef typename Lines_through_segments_mapped_2<Traits_3>::Point_2
                                                        Mapped_point_2;

  typedef typename Lines_through_segments_mapped_2<Traits_3>::General_polygon_2
                                                        Mapped_general_polygon_2;

  typedef typename Mapped_2::Mapped_transversal         Mapped_transversal;
  typedef typename Through_3::Through_transversal       Through_transversal;

   typedef typename Lines_through_segments_output_obj<
      Traits_3, Segment_3>::Transversal                 Transversal;

  typedef typename
  Lines_through_segments_output_obj<Traits_3, Segment_3>::Transversal_with_segments
  Transversal_with_segments;

  typedef typename Rational_kernel::Segment_3           Rational_segment_3;

  typedef typename std::vector< boost::shared_ptr<Arr_on_plane> >::iterator Arr_on_plane_iterator;
  typedef std::pair< Arr_on_plane_iterator, Arr_on_plane_iterator> Arr_on_plane_range;

  typedef typename std::vector< boost::shared_ptr<Arr_on_sphere> >::iterator Arr_on_sphere_iterator;
  typedef std::pair< Arr_on_sphere_iterator, Arr_on_sphere_iterator> Arr_on_sphere_range;

private:
  /// used planar arrangements
  std::vector< boost::shared_ptr<Arr_on_plane> > m_planar_arrangements;
  /// used spherical arrangements
  std::vector< boost::shared_ptr<Arr_on_sphere> > m_spherical_arrangements;

  /// generated segments
  std::vector<Rational_segment_3> m_segs;

  // kernels
  Rational_kernel                 m_rational_kernel;
  Alg_kernel                      m_alg_kernel;

public:
  /*! Constructor
   * \param alg_kernel algebraic kernel.
   * \param rational_kernel rational kernel.
   */
  Lines_through_segments_3(const Alg_kernel& alg_kernel = Alg_kernel(),
                           const Rational_kernel& rational_kernel = Rational_kernel()) :
    m_rational_kernel(rational_kernel), m_alg_kernel(alg_kernel) {}

  /*************************************************************
   * The following function gets input iterator of segments and output all
   * common lines to an output iterator.
   * A common line is a line that passes through 4 segments in 3 space.
   * In a general case only 0-2 lines will pass through 4 segments.
   *
   * A line equation of general line which intersects with 3 lines L1,L2,L3
   * is represented as an hyperbola at L1,L2 space. Intersection of 2
   * hyperbolas at this space represent a line which intersect with 4
   * segments.
   *
   * For each two segments the function gets over the n-2 segments and adds
   * for each, an object to the 2 dimensional arrangement derived by the two
   * segments.
   * Next, the function get over all the intersection points (vertices) at
   * the arrangement and output a common line for each vertex.
   *
   * The run time is O(N^3 * logn + k), where K is the output size.
   * At the best case the function will run only at O(n^3 * logn) - where
   * there aren't any common lines.
   *
   *  L2.t * L1.t * C402 - L2.t * C403 - L1.t * C401 - C400 = 0
   *
   *************************************************************/

  // Solve. The typename of the iterator is either Transversal or
  // Transversal_with_segments based on the compile-time flag
  // With_segments.
  template <typename Input_iterator,
            typename Output_iterator>
  void operator() (Input_iterator begin, Input_iterator end,
                   Output_iterator output_iterator,
                   bool rational_output = false)
  {
    typedef Lines_through_segments_impl<Traits_3, With_segments>
      Lines_through_segments_impl;

    if(begin == end) return;

    std::vector<const Rational_segment_3*> segs;
    std::transform(begin, end, std::back_inserter(segs),
                   &std::addressof<Rational_segment_3>);

    if (segs.empty()) return;

    /* Holds the number of overlapping segments.
     * If equal to the number of remaining segments, than all of the
     * remanining segments overlap,
     * and the common line for each four is returned.
     */
    unsigned int num_of_overlap_segs = 0;
    unsigned int S1, S2, S3;

    for (S1 = 0; (S1+3) != segs.size(); ++S1)
    {
      for (S2 = S1+1; (S2+2) != segs.size(); ++S2)
      {
#if LINES_DEBUG
        std::cout << change_color(CGAL_YELLOW, "NEW ARRANGEMENT S1 = ",
                                  *segs[S1], "\nS2 = ", *segs[S2])
                  << std::endl;
#endif
        try
        {
          Lines_through_segments_impl line_through_segs_obj(
            segs[S1], segs[S2], &m_alg_kernel, &m_rational_kernel);

          /* For each line add a new curve to the arrangement, the
           * intersection of the hyperbolas represent common lines. */
          for (S3 = S2+1; S3 != segs.size(); ++S3)
          {
#if LINES_DEBUG
            std::cout <<  change_color(CGAL_GREEN,
                                       "\nAdd Element to Arr = ",
                                       *segs[S3])
                      << std::endl;
#endif
            line_through_segs_obj.add_element_to_arrangement(*segs[S3]);
          }

          line_through_segs_obj.find_all_lines(rational_output, output_iterator);

          // extract the arrangments
          if(!line_through_segs_obj.arrangement_on_plane()->is_empty())
            m_planar_arrangements.push_back(line_through_segs_obj.arrangement_on_plane());
          if(!line_through_segs_obj.arrangement_on_sphere()->is_empty())
            m_spherical_arrangements.push_back(line_through_segs_obj.arrangement_on_sphere());

          num_of_overlap_segs = 0;
        }
        catch(const CGAL::Lines_through_segments_exp_2_lines_overlap& /* e */)
        {
          const Rational_segment_3* temp = segs[S2];
          segs[S2] = segs[segs.size() - 1];
          segs[segs.size() - 1] = temp;
#if LINES_DEBUG
          std::cout << e.what() << std::endl;
          std::cout << change_color(CGAL_YELLOW, "2 lines over = (",
                                    *segs[S1],",",
                                    *segs[S2],")")
                    << std::endl;
#endif
          --S2;
          ++num_of_overlap_segs;
          if (num_of_overlap_segs == segs.size() - S1)
          {
            // XXX we really need a Traits_3 value around here
            LTS::find_overlap(Traits_3(), segs.begin() + S1, segs.end(), output_iterator, With_segments());
            return;
          }
        }
      }
    }
  }

  Arr_on_plane_range
  planar_arrangements() {
    return std::make_pair(m_planar_arrangements.begin(),
                          m_planar_arrangements.end());
  }

  Arr_on_sphere_range
  spherical_arrangements() {
    return std::make_pair(m_spherical_arrangements.begin(),
                          m_spherical_arrangements.end());
  }

};
} //namespace CGAL

#endif /*LINE_THROUGH_SEGMENTS_3_H*/
