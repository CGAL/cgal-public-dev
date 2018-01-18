#ifndef ARRANGEMENT_GENERAL_FUNCTIONS_H
#define ARRANGEMENT_GENERAL_FUNCTIONS_H

#include "Graphics.h"

#include <stdio.h>
#include <stdlib.h>
#include <CGAL/basic.h>
#include <CGAL/Arr_enums.h>

const int SAMPLE_RATE = 500;

template <typename Rational_kernel,
          typename Alg_kernel,
          typename Arrangement_2,
          typename Traits_2_adapt,
          typename Point_2_adapt>
class Arrangement_general_functions
{
   typedef typename Arrangement_2::Traits_2 Traits_2;
   typedef typename Rational_kernel::Point_2 RPoint_2;
   typedef class Graphics<typename Rational_kernel::FT,
                          typename Traits_2::Curve_2,
                          RPoint_2> Graphics;
   typedef std::vector<RPoint_2> PointVector;			
   Graphics *m_graphics;
   char *argv[1];
   Traits_2_adapt m_traits_2_adapt;
   
public:   
   Arrangement_general_functions()
   {

      argv[0] = (char*)malloc(20);
      strcpy(argv[0],"ASAF DEBUG\n");
      
      m_graphics = new Graphics(1, argv);
   }

   ~Arrangement_general_functions()
   {
      free(argv[0]);
      delete m_graphics;
   }
   


private:

   void draw_arrangement_curve(const typename Traits_2::X_monotone_curve_2& curve, char* color)
   {
      PointVector	 p;
      sample_x_monotone_curve<PointVector,RPoint_2>(curve,p);

      /* Patch - the left and right are not working properly for vertical 
         lines. Remove points - blue rectangles.
      */
      if (p.size() != 2 || (p[0] != p[1]))
         draw_poly_segment(p.begin(), p.end(), color) ;
      return;
   }

public:
      void operator() (Arrangement_2& arr )
   {
      typedef typename Rational_kernel::FT Rational;
      typename Arrangement_2::Edge_iterator iter;
      for (iter = arr.edges_begin (); iter != arr.edges_end (); ++iter)
      {
         if ((iter->twin()->face()->num_of_overlap_plane_faces() == 1 &&
              iter->face()->num_of_overlap_plane_faces() == 1 &&
              *iter->face()->segs_begin() != 
              *iter->twin()->face()->segs_begin()) ||
             (iter->twin()->face()->num_of_overlap_plane_faces() == 1 &&
              iter->face()->num_of_overlap_plane_faces() == 1 &&
              *iter->face()->segs_begin() != 
              *iter->segs_begin()) ||
             (arr.number_of_originating_curves(iter) >= 2 &&
              iter->twin()->face()->num_of_overlap_plane_faces() < 2 &&
              iter->face()->num_of_overlap_plane_faces() < 2))
         {
            this->draw_arrangement_curve(iter->curve(), "black"); //blue
         }
         else
         {
            this->draw_arrangement_curve(iter->curve(), "black");
         }
                  
         

      }

      typename Arrangement_2::Face_iterator fiter;
      for (fiter = arr.faces_begin (); fiter != arr.faces_end (); ++fiter)
      {
         if (!fiter->is_unbounded())
         {
            Rational sum_x = 0;
            Rational sum_y = 0;
            int size = 0;
            typename Traits_2::X_monotone_curve_2 left_most;
            double left_x;
            typename Arrangement_2::Ccb_halfedge_circulator circ,curr;
            curr = circ = fiter->outer_ccb();
            
            left_most = circ->curve();
            left_x = CGAL::to_double(left_most.source().x()) + CGAL::to_double(left_most.target().x());
            
            do {
               PointVector	 p;
               sample_x_monotone_curve<PointVector,RPoint_2>(curr->curve(),p);
               typename PointVector::iterator pit;

               for (pit = p.begin();pit!= p.end(); ++pit)
               {
                  sum_x += pit->x();
                  sum_y += pit->y();
                  size++;
               }
               if (CGAL::to_double(curr->curve().source().x()) + 
                   CGAL::to_double(curr->curve().target().x()) < left_x)
               {
                  left_most = curr->curve();
                  left_x = CGAL::to_double(left_most.source().x()) + CGAL::to_double(left_most.target().x());

               }
               
            } while (++curr != circ);
            
            sum_x /= size;
            sum_y /= size;
            
            std::stringstream n_faces;
            n_faces << fiter->num_of_overlap_plane_faces();

            if (left_most.source().x() == left_most.target().x())
            {
               Point_2_adapt source = left_most.source();
               Point_2_adapt target = left_most.target();

               m_graphics->add_text(n_faces.str().c_str(),
               to_screen_double_x(source.x()),
               to_screen_double_y((
                  source.y() + target.y())/2));




               // m_graphics->draw_cross(to_screen_double_x(left_most.source().x()) + 2,
               //                        to_screen_double_y((left_most.source().y() + 
               //                                            left_most.target().y())/2),c);
            }
            else
            {
               m_graphics->add_text(n_faces.str().c_str(),
                                    CGAL::to_double(sum_x),CGAL::to_double(sum_y));
//               m_graphics->draw_cross(CGAL::to_double(sum_x),CGAL::to_double(sum_y),c);
            }
         }
      }

      m_graphics->Display();
      return;
   }

   
private:

   template <typename Point_vector, typename Point_2>
   void sample_non_vertical_curve(const typename Traits_2::X_monotone_curve_2& curve,Point_vector& p)
   {
      typedef typename Alg_kernel::FT Algebraic;
      typedef typename Rational_kernel::FT Rational;

      Point_2_adapt source = curve.source();
      Point_2_adapt target = curve.target();

      Algebraic sx = source.x();
      Algebraic tx = target.x();

      /* t x y + u x + v y + w = 0. */
      /* t x y + v y = - u *x - w. */
      /* y  = - (u * x + w)/(t * x + v). */
   
      //sample points along the curve

      Rational n ;
      Rational step (1,SAMPLE_RATE) ;
      bool first = true;
	
      for (n = 0; n < 1; n +=step)
      {
         Algebraic x;
         Algebraic y;
         if ( n + step < 1)
            x = sx + n * (tx - sx);
         else
            x = tx ;
          y = m_traits_2_adapt.get_y_val(curve,x);
         
// = - (curve.u() * x + curve.w()) / 
//             (curve.t() * x + curve.v());
         p.push_back(Point_2(to_screen_double_x(x),to_screen_double_y(y)));
      }	
      return;
   }
      template <typename Point_vector, typename Point_2, typename Algebraic>
      void sample_vertical_curve(const typename Traits_2::X_monotone_curve_2& curve,Point_vector& p,
      Algebraic& source_x, Algebraic& source_y, Algebraic& target_x, Algebraic& target_y)
   {
      Point_2 s(to_screen_double_x(source_x), to_screen_double_y(source_y));
      Point_2 t(to_screen_double_x(target_x), to_screen_double_y(target_y));

      p.push_back(s);
      p.push_back(t);
            
      return;
   }

   template <typename Point_vector, typename Point_2>
   void sample_x_monotone_curve(const typename Traits_2::X_monotone_curve_2& curve,Point_vector& p)
   {
      typedef typename Alg_kernel::FT Algebraic;

      Algebraic source_x;
      Algebraic source_y;
      Algebraic target_y;
#if USE_RATIONAL_ARC_TRAITS      
      if (curve.left_infinite_in_x() == CGAL::ARR_INTERIOR)
      {
         if (curve.left_infinite_in_y() == CGAL::ARR_INTERIOR)
#endif
         {
            Point_2_adapt source = curve.left();
            source_x = source.x();
            source_y = source.y();
         }
#if USE_RATIONAL_ARC_TRAITS
         else
         {
            source_y = - 100;
            source_x = m_traits_2_adapt.convert_real_to_algebraic(
               curve.left_x());
         }
      }
      else
      {
         m_traits_2_adapt.get_horizontal_asymptote_y_val(curve, source_y);
    
         source_x = -100;
      }
#endif
      
      Algebraic target_x;
#if USE_RATIONAL_ARC_TRAITS
      if (curve.right_infinite_in_x() == CGAL::ARR_INTERIOR)
      {
         if (curve.right_infinite_in_y() == CGAL::ARR_INTERIOR)
#endif
         {
            Point_2_adapt target = curve.right();
            target_x = target.x();
            target_y = target.y();
         }
#if USE_RATIONAL_ARC_TRAITS
         else
         {
            target_y = 100;
            target_x = m_traits_2_adapt.convert_real_to_algebraic(
               curve.right_x());
         }
      }
      else
      {
         target_x = 100;
         m_traits_2_adapt.get_horizontal_asymptote_y_val(curve, target_y);
      }
#endif
      if (source_x == target_x)
      {
         sample_vertical_curve<Point_vector,Point_2>(curve,p,
         source_x,source_y,target_x, target_y);
      }
      else 
      {
         if (target_y == source_y)
         {
            sample_vertical_curve<Point_vector,Point_2>(curve,p,
            source_x,source_y,target_x, target_y);
         }
         else  /*!Is_vertical_2 */
         {
            sample_non_vertical_curve<Point_vector,Point_2>(curve,p);
         }
      }
   }


   
   template <typename Input_iterator>
   void draw_poly_segment(Input_iterator begin, Input_iterator end, char* color)
   {
      Input_iterator curr_iterator = begin;
      RPoint_2 prev(*curr_iterator);
      ++curr_iterator;
      RPoint_2 curr(*curr_iterator);
      for ( ; curr_iterator != end ; ++curr_iterator)
      {
         curr = *curr_iterator;
         m_graphics->draw_edge(prev,curr,color);
                  
         prev = curr;
      }
      return;
   }
private:
   template <typename Number_type>
   double to_screen_double_y(const Number_type y)
   {
      return ((1 - CGAL::to_double(y)));
   }

   template <typename Number_type>
   double to_screen_double_x(const Number_type x)
   {
      return ((CGAL::to_double(x)));
   }

};
   
#endif	//ARRANGEMENT_GENERAL_FUNCTIONS_H
