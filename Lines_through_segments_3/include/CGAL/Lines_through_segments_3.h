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
#include <boost/variant.hpp>
#include <CGAL/Lines_through_segments_impl.h>
#include <CGAL/Lines_through_segments_find_overlap_lines.h>
#include <CGAL/Lines_through_segments_exceptions.h>
#include <CGAL/Lines_through_segments_output_obj.h>

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
  typedef boost::false_type                             With_arrangement;
   
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
  // typedef typename Lines_through_segments_mapped_2<Traits_3>::Rational_point_2
  // Mapped_rat_point_2;
   
  typedef typename Lines_through_segments_mapped_2<Traits_3>::General_polygon_2  
                                                        Mapped_general_polygon_2;

  typedef typename Mapped_2::Mapped_transversal         Mapped_transversal;
  typedef typename Through_3::Through_transversal       Through_transversal;
      
   typedef typename Lines_through_segments_output_obj<
      Traits_3, 
      Segment_3>::Transversal   Transversal;

  typedef typename
  Lines_through_segments_output_obj<Traits_3, Segment_3>::Transversal_with_segments 
  Transversal_with_segments;

   typedef typename Rational_kernel::Segment_3           Rational_segment_3;
   typedef typename 
   Lines_through_segments_mapped_2_with_arrangement<Traits_3,
                                                    Rational_segment_3>::Arrangement_2
   Arr_on_plane;

   typedef typename 
   Lines_through_segments_through_3_with_arrangement<Traits_3,
                                                     Segment_3>::Arrangement_2
   Arr_on_sphere;
   
protected:
   std::vector<Rational_segment_3> m_segs;
   Rational_kernel                 m_rational_kernel;
   Alg_kernel                      m_alg_kernel;
   const Rational_kernel*          m_rational_kernel_ptr;
   const Alg_kernel*               m_alg_kernel_ptr;

   
public:
  /*! Empty constructor */
  Lines_through_segments_3() {}
      
  /*! Constructor
   * \param alg_kernel algebraic kernel.
   * \param rational_kernel rational kernel.
   */
  Lines_through_segments_3(const Alg_kernel& alg_kernel,
                           const Rational_kernel& rational_kernel) :
    m_rational_kernel_ptr(&rational_kernel),
    m_alg_kernel_ptr(&alg_kernel)
  {
  }

  /*! Destructor */
  ~Lines_through_segments_3()
  {
    m_segs.clear();
  }
  
  /*************************************************************
   * The following function gets input iterator of segments and output all
   * common lines to an output iterator.
   * A common line is a line that passes through 4 segments in 3 space.
   * In a general case only 0-2 lines will pass through 4 segments.
   * 
   * A line equation of general line which intersects with 3 lines L1,L2,L3
   * is represented as an hyperbola at L1,L2 space. Intersection of 2
   * hyperbolas at this space represent a line which intersect with 4 
   * segmetns.
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
                   bool copy_segs_container = true,
                   bool rational_output = false)
  {
     this->get_all_segments(begin, end, 
                            output_iterator,
                            With_arrangement(),
                            copy_segs_container,
                            rational_output);
  }
   
protected:   
  template <typename Input_iterator,
            typename Output_iterator, 
            typename With_arrangement>
  void get_all_segments(Input_iterator begin, Input_iterator end,
                        Output_iterator output_iterator,
                        With_arrangement with_arrangement,
                        bool copy_segs_container,
                        bool rational_output)
  {
    typedef Lines_through_segments_impl<Traits_3, Output_iterator,
                                        With_segments, With_arrangement>
      Lines_through_segments_impl;
      
     if (copy_segs_container)
     {
        m_segs.clear();
        for (Input_iterator it = begin; it != end; ++it) 
           m_segs.push_back(*it);
     }
  
    std::vector<const Rational_segment_3*> segs;

    if (copy_segs_container)
    {
       for (unsigned int index = 0; index < m_segs.size(); ++index)
       {
          segs.push_back(&m_segs[index]);
       }
    }
    else
    {
       for (Input_iterator it = begin; it != end; ++it) 
          segs.push_back(&(*it));
    }
    
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
          Lines_through_segments_impl
            line_through_segs_obj(segs[S1], segs[S2], &output_iterator,
                                  m_alg_kernel_ptr,m_rational_kernel_ptr);
                    
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

          line_through_segs_obj.find_all_lines(rational_output);

          if (With_arrangement())
          {
             this->set_arrangements(line_through_segs_obj.arrangement_on_plane(),
                                    line_through_segs_obj.arrangement_on_sphere());
          }
#if LTS_DRAW_ARR
          line_through_segs_obj.draw_arr();
#endif
          num_of_overlap_segs = 0;
        }
        catch(const CGAL::Lines_through_segments_exp_2_lines_overlap& e)
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
            Lines_through_segments_find_overlap_lines<Traits_3, With_segments>
              find_overlap;
            find_overlap(S1, segs, &output_iterator);
            return;
          }
        }
      }
    }
  }
  virtual void set_arrangements(Arr_on_plane* arr_on_plane,
                                 Arr_on_sphere* arr_on_sphere)
  {
  }
};

template <typename Lines_through_segments_traits_3,
          typename With_segments_ = boost::false_type>
class Lines_through_segments_3_with_arrangements :
      public Lines_through_segments_3<Lines_through_segments_traits_3,
                                      With_segments_>
{
private:
  typedef Lines_through_segments_3<Lines_through_segments_traits_3,
                                   With_segments_>              Base;
  typedef
  Lines_through_segments_3_with_arrangements<Lines_through_segments_traits_3,
                                             With_segments_>    Self;

public:
  typedef typename Base::Traits_3                       Traits_3;
  typedef With_segments_                                With_segments;
  typedef boost::true_type                              With_arrangement;
  typedef typename Base::Rational_kernel                Rational_kernel;
  typedef typename Rational_kernel::Segment_3           Segment_3;
  typedef typename Base::Alg_kernel                     Alg_kernel;

  typedef Lines_through_segments_output_obj<Traits_3, 
                                            Segment_3>   LTS_output_obj;
   
   
  /* Output objects */
  typedef typename Base::Line_3              Line_3;
  typedef typename LTS_output_obj::Mapped_2_with_arr    Mapped_2;
  typedef typename Mapped_2::Point_2                    Mapped_point_2;
//  typedef typename Mapped_2::Rational_point_2           Mapped_rat_point_2;
  typedef typename Mapped_2::X_monotone_curve_2
    Mapped_x_monotone_curve_2;
  typedef typename Mapped_2::General_polygon_2
    Mapped_general_polygon_2;

  typedef typename LTS_output_obj::Through_3_with_arr   Through_3;

  typedef typename Through_3::Point_3                   Through_point_3;
  typedef typename Through_3::Segment_3                 Through_segment_3;
  typedef typename Through_3::Point_3_segment_3
    Through_point_3_segment_3;
   
  typedef typename Mapped_2::Mapped_transversal         Mapped_transversal;
  typedef typename Through_3::Through_transversal       Through_transversal;
   
  typedef typename LTS_output_obj::Transversal_with_arr Transversal;
   
  typedef typename LTS_output_obj::Transversal_with_segments_with_arr 
    Transversal_with_segments;

private:

  std::list<typename Base::Arr_on_plane*> m_planar_arrangements;
  std::list<typename Base::Arr_on_sphere*> m_spherical_arrangements;

public:
  // Solve. The typename of the iterator is either Transversal or
  // Transversal_with_segments based on the compile-time flag
  // With_segments.
  template <typename Input_iterator,
            typename Output_iterator>
  void operator() (Input_iterator begin, Input_iterator end,
                   Output_iterator output_iterator,
                   bool copy_segs_container = true,
                   bool rational_output = false)
  {
     
    Base::get_all_segments(begin, end, 
                            output_iterator, 
                           With_arrangement(),
                           copy_segs_container,
                           rational_output);
  }
   
  /*! Empty constructor */
  Lines_through_segments_3_with_arrangements() 
    : Base() 
  {}
   
  ~Lines_through_segments_3_with_arrangements()
  {
    typename std::list<typename Base::Arr_on_plane* >::iterator it_planar;
    for (it_planar = m_planar_arrangements.begin();
         it_planar != m_planar_arrangements.end();
         ++it_planar)
    {
      delete *it_planar;
    }
      
    typename std::list<typename Base::Arr_on_sphere* >::iterator it_spherical;
    for (it_spherical = m_spherical_arrangements.begin();
         it_spherical != m_spherical_arrangements.end();
         ++it_spherical)
    {
      delete *it_spherical;
    }
      
  }
   
  Lines_through_segments_3_with_arrangements(const Alg_kernel& alg_kernel,      
                                             const Rational_kernel&
                                             rational_kernel) :
     Base(alg_kernel, rational_kernel)
  {
  }
   
private:
  void set_arrangements(typename Base::Arr_on_plane* arr_on_plane,
                        typename Base::Arr_on_sphere* arr_on_sphere)
  {
    m_planar_arrangements.push_back(arr_on_plane);
    m_spherical_arrangements.push_back(arr_on_sphere);
  }
};
   
} //namespace CGAL

#endif /*LINE_THROUGH_SEGMENTS_3_H*/
