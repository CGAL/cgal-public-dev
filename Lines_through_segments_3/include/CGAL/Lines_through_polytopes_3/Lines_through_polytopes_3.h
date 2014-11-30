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

#ifndef LINE_THROUGH_POLYTOPES_3_H
#define LINE_THROUGH_POLYTOPES_3_H

/*! \file
*************************************************************
* The following class computes all of the lines that crosses 4 polytopes in 3D.
* The class is a functor. 
*
* Input:
*
* Lines_through_segments_traits_3 - dedicated traits class.
*
*************************************************************
*/

#include <CGAL/Lines_through_polytopes_3/Lines_through_polytopes_impl.h>
namespace CGAL {

template <typename Lines_through_segments_traits_3>
class Lines_through_polytopes_3 {
      
public:
  typedef Lines_through_segments_traits_3      Traits_3;
  typedef typename Traits_3::Rational_kernel   Rational_kernel;
  typedef typename Traits_3::Alg_kernel        Alg_kernel;

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
  typedef CGAL::Polyhedron_3<Rational_kernel>            Rational_polyhedron_3;
  typedef typename Lines_through_segments_output_obj<Traits_3,
                                                     Rational_polyhedron_3>::Transversal
                                                     Transversal;

  typedef typename
  Lines_through_segments_output_obj<Traits_3, Rational_polyhedron_3>::Transversal_with_segments 
  Transversal_with_segments;
      
protected:
   Rational_kernel                 m_rational_kernel;
   Alg_kernel                      m_alg_kernel;
   const Rational_kernel*          m_rational_kernel_ptr;
   const Alg_kernel*               m_alg_kernel_ptr;

public:

  /*! Empty constructor */
  Lines_through_polytopes_3() {}
      
  /*! Constructor
   * \param alg_kernel algebraic kernel.
   * \param rational_kernel rational kernel.
   */
  Lines_through_polytopes_3(const Alg_kernel& alg_kernel,
                           const Rational_kernel& rational_kernel) :
    m_rational_kernel_ptr(&rational_kernel),
    m_alg_kernel_ptr(&alg_kernel)
  {
  }

  /*! Destructor */
  ~Lines_through_polytopes_3()
  {
  }

  /*************************************************************
   * The following function gets input iterator of polytopes and output all
   * common lines to an output iterator.
   * A common line is a line that passes through 4 polytopes in 3 space.
   * 
   *************************************************************/

  template <typename Input_iterator, typename Insert_iterator>
  void operator() (Input_iterator polytopes_start,
                   Input_iterator polytopes_end,
                   Insert_iterator insert_iterator)
  {
    typedef Lines_through_polytopes_impl<Lines_through_segments_traits_3,
      Insert_iterator>
      Lines_through_polytopes_impl;

    typedef typename Rational_kernel::Segment_3            Rational_segment_3;
    typedef typename Rational_kernel::Point_3              Rational_point_3;
    typedef CGAL::Polyhedron_3<Rational_kernel>            Rational_polyhedron_3;
    typedef typename Rational_polyhedron_3::Edge_iterator  Edge_iterator;
    int ii = 0;
    int jj = 0;
         
         
    for (; polytopes_start != polytopes_end; ++polytopes_start)
    {
      Input_iterator temp = polytopes_start;
      std::advance(temp,3);
      if (temp ==  polytopes_end)
        break;

      Input_iterator polytopes_start_second = polytopes_start;
      polytopes_start_second++;
            
      for (; polytopes_start_second != polytopes_end; polytopes_start_second++)
      {
        temp = polytopes_start_second;
        std::advance(temp,2);
        if (temp == polytopes_end)
          break;

        std::cout << change_color(CGAL_YELLOW,"NEW COUPLE OF POLYTOPES") << std::endl;
        std::cout << *polytopes_start << std::endl;
        std::cout << *polytopes_start_second << std::endl;
               
        for (Edge_iterator e_p1 = polytopes_start->edges_begin(); 
             e_p1 != polytopes_start->edges_end();
             ++e_p1)
        {
          for (Edge_iterator e_p2 = polytopes_start_second->edges_begin(); 
               e_p2 != polytopes_start_second->edges_end();
               ++e_p2)
          {
            // #if 0
            //                      std::cout << "jj = " << jj << std::endl;
            std::cout 
              << change_color(CGAL_GREEN,
                              "New Arrangement\ne_p1 = ",
                              e_p1->vertex()->point().x()," ",
                              e_p1->vertex()->point().y()," ",
                              e_p1->vertex()->point().z()," ",
                              e_p1->opposite()->vertex()->point().x()," ",
                              e_p1->opposite()->vertex()->point().y()," ",
                              e_p1->opposite()->vertex()->point().z()," ",
                              "\ne_p2 = ",
                              e_p2->vertex()->point().x()," ",
                              e_p2->vertex()->point().y()," ",
                              e_p2->vertex()->point().z()," ",
                              e_p2->opposite()->vertex()->point().x()," ",
                              e_p2->opposite()->vertex()->point().y()," ",
                              e_p2->opposite()->vertex()->point().z())
              << std::endl;

            //#endif                
            //                     if (jj == 0)//6 && jj != 9 && jj != 19 && jj != 24 && jj < 27)
            {
              Lines_through_polytopes_impl
                lines_through_poly_obj(*polytopes_start,e_p1,
                                       *polytopes_start_second, e_p2, &insert_iterator,
                                       m_alg_kernel_ptr,m_rational_kernel_ptr);

              if (lines_through_poly_obj.is_valid())
              {
                Input_iterator polytopes_start_third = polytopes_start_second;
                polytopes_start_third++;
                           
                if (polytopes_start_third == polytopes_end)
                  return;
                           
                for (polytopes_start_third;
                     polytopes_start_third != polytopes_end;
                     polytopes_start_third++)
                {
                  std::cout << change_color(CGAL_BLUE,"ADD POLYTOPE") << std::endl;
#if 0
                  std::cout << *polytopes_start_third << std::endl;
#endif                         
                  lines_through_poly_obj.add_polytope(*polytopes_start_third);
                }
                lines_through_poly_obj.find_all_lines();
#if LTS_DRAW_ARR
                lines_through_poly_obj.draw_arr();
                lines_through_poly_obj.draw_rest_arr();
#endif
              }
              else
              {
                std::cout << change_color(CGAL_RED,"INVALID") << std::endl;
              }
            }
            jj++;
                     
          }
                  
                  
          ii++;
                  
        }

      }
    }
         
  }
   template <typename Arr_on_plane,
             typename Arr_on_sphere>
   void set_arrangements(Arr_on_plane* arr_on_plane,
                         Arr_on_sphere* arr_on_sphere)
   {
   }
   
};
   
} //namespace CGAL

#endif /*LINE_THROUGH_POLYTOPES_3_H*/
