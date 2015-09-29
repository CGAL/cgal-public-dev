// Copyright (c) 2006 Tel-Aviv University (Israel).
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
// $URL: svn+ssh://ophirset@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_2/include/CGAL/IO/Fig_stream_Conic_arc_2.h $
// $Id: Fig_stream_Conic_arc_2.h 5459 2007-11-29 10:24:13Z ophirset $
// 
// Author(s)     : Ron Wein           <wein@post.tau.ac.il>

#ifndef CGAL_FIG_STREAM_CONIC_ARC_2_H
#define CGAL_FIG_STREAM_CONIC_ARC_2_H

#include <CGAL/IO/Fig_stream.h>
#include <CGAL/Arr_enums.h>
#include <list>

namespace CGAL {

/*!
 * Write an x-monotone conic arc to a FIG stream.
 */
template <class Conic_traits>
void write_x_monotone_conic_arc 
(CGAL::Fig_stream<typename Conic_traits::Rat_kernel>& fs,
 const typename Conic_traits::X_monotone_curve_2&     cv)
{
  typedef typename Conic_traits::Rat_kernel  Rat_kernel;
  typedef typename Rat_kernel::Segment_2     Segment_2;
  typedef typename Rat_kernel::Point_2       Point_2;

  Bbox_2 bbox (CGAL::to_double(fs.bounding_rect().xmin()), 
               CGAL::to_double(fs.bounding_rect().ymin()), 
               CGAL::to_double(fs.bounding_rect().xmax()), 
               CGAL::to_double(fs.bounding_rect().ymax()));
  
  // Draw a curves conic arc: As the arc is x-monotone, its source and its 
  // target has the extreme x-coordinates.
  int     n = min(fs.width(), 100 );
  if (cv.orientation() == COLLINEAR)
  {
    // In case of a line segment we need 1 sample points.
    n = 1;
  }
  
  if (n <= 0)
    return;
  
  typedef std::pair<double, double>    App_point_2;
  int                                  i;    
  
  App_point_2  *pts = new App_point_2 [n + 1];
  cv.polyline_approximation (n, bbox, pts);

/*=======
    Conic_traits traits;
    Arr_parameter_space left_boundary = traits.parameter_space_in_x_2_object() (cv, MIN_END);
    Arr_parameter_space right_boundary = traits.parameter_space_in_x_2_object() (cv, MAX_END);
    
    if (cv.orientation() == CGAL::COLLINEAR)
    {
        if (left_boundary==ARR_INTERIOR && right_boundary==ARR_INTERIOR)
        {
            // In case of a linear segment:
            Alg_segment_2   seg = Alg_segment_2 (cv.source(), cv.target());
            fs << seg;
        }
        else if (left_boundary!=ARR_INTERIOR && right_boundary!=ARR_INTERIOR)
        {
            // In case of a linear line:
            Alg_line_2   line (cv.source(), cv.target());
            fs << line;
        }
        else
        {
            // In case of a ray.
            Alg_ray_2   ray;
            if (left_boundary!=ARR_INTERIOR)
            {
                ray = Alg_ray_2 (cv.right(), cv.left());
            }
            else
            {
                CGAL_assertion (right_boundary!=ARR_INTERIOR);
                ray = Alg_ray_2 (cv.left(), cv.right());
            }
>>>>>>> .r41166
*/

  for (i = 1; i <= n; i++)
  {
    const double r = 1000.0;
    double s_x = long (pts[i-1].first * r) / r;
    double s_y = long (pts[i-1].second * r) / r;
    double t_x = long (pts[i].first * r) / r;
    double t_y = long (pts[i].second * r) / r;
    if ((s_x == t_x) and (s_y == t_y))
      continue;

/*=======
            fs << ray;
        }
    } 
    else if (CGAL::compare (cv.r(), cv.s()) == CGAL::EQUAL &&
             CGAL::sign (cv.t()) == CGAL::ZERO)
    {
        // In case of a circular arc:
        Algebraic    x_mid = (cv.source().x() + cv.target().x()) / 2;
        Alg_point_2  q = Alg_point_2(x_mid, 0);
        Alg_point_2  p = cv.get_point_at_x (q);
        
        fs.write_circular_arc (cv.source(), p, cv.target());
    }
    else
    {
        Alg_kernel ker;
        Algebraic  x_start = cv.left().x();
        if (cv.is_left_unbounded())
        {
            x_start = fs.bounding_rect().xmin();
            
            // vertical asymptote
            if (left_boundary==ARR_INTERIOR)
            {
                Alg_point_2 p = ker.construct_point_on_2_object() 
                    (cv.get_vertical_asymptote(), 0);
                Algebraic x_asymptote = ker.compute_x_2_object() (p);
                x_asymptote += 1e-6;
                x_start = max(x_start, x_asymptote);
            }
        }
>>>>>>> .r41166
*/

    Segment_2 seg (Point_2(s_x, s_y), Point_2(t_x, t_y));
    fs << seg;
  }
  delete[] pts;

/*=======
        Algebraic  x_end = cv.right().x();
        if (cv.is_right_unbounded())
        {
            x_end = fs.bounding_rect().xmax();
            
            // vertical asymptote
            if (right_boundary==ARR_INTERIOR)
            {
                Alg_point_2 p = ker.construct_point_on_2_object() 
                    (cv.get_vertical_asymptote(), 0);
                Algebraic x_asymptote = ker.compute_x_2_object() (p);
                x_asymptote -= 1e-6;
                x_end = min(x_end, x_asymptote);
            }
        }
 
        
        // Represent the arc as a spline with 5 control points.
        Algebraic   x;
        Alg_point_2 q;
        Alg_point_2 cps[5];
        int         i;
        
        for (i = 0; i < 5; i++)
        {
            x = (x_start*(4 - i) + x_end*i) / 4;
            
            q = Alg_point_2(x, 0);
            cps[i] = cv.get_point_at_x (q);      
        }
        
        fs.write_spline (cps + 0, cps + 5, 1.0);
    }
    
    return;
>>>>>>> .r41166
*/
}

/*!
 * Write a conic arc to a FIG stream.
 */
template <class Conic_traits>
void write_conic_arc
(CGAL::Fig_stream<typename Conic_traits::Rat_kernel>& fs,
 const typename Conic_traits::Curve_2&                cv)
{
    typedef typename Conic_traits::X_monotone_curve_2   Conic_arc_2;
    
    // Subdivide the arc into x-monotone sub-arcs.
    Conic_traits           traits;
    std::list<Object> xcvs;
    
    traits.make_x_monotone_2_object() (cv,
                                       std::back_inserter(xcvs));
    
    // Write the x-monotone sub-arcs.
    typename std::list<Object>::iterator xit;
    
    for (xit = xcvs.begin(); xit != xcvs.end(); ++xit)
    {
        Conic_arc_2 arc;
        if (assign(arc, *xit))
            write_x_monotone_conic_arc<Conic_traits> (fs, arc);
    }
    
    return;
}

} //namespace CGAL

#endif
