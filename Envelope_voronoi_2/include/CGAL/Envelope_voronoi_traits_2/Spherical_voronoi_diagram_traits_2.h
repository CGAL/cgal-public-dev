// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
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
// $Id:  $
// 
//
// Author(s): Ophir Setter          <ophir.setter@post.tau.ac.il>
//

/*!
  \file   Spherical_voronoi_diagram_traits_2.h
  \brief  The file contains a traits class for the computation of spherical
  Voronoi diagram on the sphere (VD of points). The points have to reside 
  exactly on the unit sphere. This traits class should be united with the 
  traits class for the power diagram on the sphere.
*/


#ifndef CGAL_SPHERICAL_VORONOI_DIAGRAM_TRAITS_2_H
#define CGAL_SPHERICAL_VORONOI_DIAGRAM_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/functions_on_enums.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Envelope_3/Envelope_base.h>

#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>


namespace CGAL {

template <typename T_Kernel>
class Spherical_voronoi_diagram_traits_2 
: public Arr_geodesic_arc_on_sphere_traits_2<T_Kernel>
{
 public:
  typedef T_Kernel                              Kernel;
  typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel>    Base;
  
  // Category tags:
  typedef typename Base::Has_left_category      Has_left_category;
  typedef typename Base::Has_merge_category     Has_merge_category;

  typedef typename Base::Left_side_category     Left_side_category;
  typedef typename Base::Bottom_side_category   Bottom_side_category;
  typedef typename Base::Top_side_category      Top_side_category;
  typedef typename Base::Right_side_category    Right_side_category;
  

  // Traits objects
  typedef typename Base::Point_2                Point_2;
  typedef typename Base::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Base::Curve_2                Curve_2;
  typedef typename Base::Multiplicity           Multiplicity;
    
  typedef Point_2                               Site_2;

  typedef typename Base::Plane_3                Plane_3;

  typedef typename Kernel::Ray_3                Ray_3;
  typedef typename Kernel::Point_3              Point_3;
  typedef typename Kernel::Direction_3          Direction_3;
  typedef typename Kernel::Vector_3             Vector_3;
  
  static inline Plane_3 bisector(const Point_2 &p1, const Point_2 &p2)
  {
    // The function currently works for nomalized DIFFERENT directions only and 
    // not general directions to avoid the need for sqrt.
    Kernel kernel;

    CGAL_envelope_voronoi_precondition (kernel.equal_3_object()(p1, p2) != true);
    CGAL_envelope_voronoi_precondition_msg (p1.dx()*p1.dx() + p1.dy()*p1.dy() + p1.dz()*p1.dz() == 1, "The directions must be normalized.");
    CGAL_envelope_voronoi_precondition_msg (p2.dx()*p2.dx() + p2.dy()*p2.dy() + p2.dz()*p2.dz() == 1, "The directions must be normalized.");

//    typedef typename Kernel::Vector_3     Vector_3;
//    typedef typename Kernel::FT           FT;

//     Vector_3 v1 = p1.vector();
//     Vector_3 v2 = p2.vector();

//     // The resulting vector is the normal to the plane that is the bisector
//     // of these two normalized directions (The plane should be directed 
//     // towards p1).
//     Vector_3 diff = v1 - v2;
//     FT x = kernel.compute_x_3_object() (diff);
//     FT y = kernel.compute_y_3_object() (diff);
//     FT z = kernel.compute_z_3_object() (diff);
    
//     typename Kernel::Plane_3 p = 
//       kernel.construct_plane_3_object() (x, y, z, 0);

    typename Kernel::Point_3 q1 (p1.dx(), p1.dy(), p1.dz());
    typename Kernel::Point_3 q2 (p2.dx(), p2.dy(), p2.dz());
    
    return kernel.construct_bisector_3_object()(q1, q2);
  }

public:
  class Compare_distance_at_point_2
  {
  public:
    
    Comparison_result operator()(const Site_2& h1,
                                 const Site_2& h2,
                                 const Point_2& p) const
    {
      Plane_3 bis = bisector(h1, h2);
      
      Kernel kernel;
      typename Kernel::Ray_3 r = 
        kernel.construct_ray_3_object()(ORIGIN, p);
      typename Kernel::Vector_3 v = 
        kernel.construct_vector_3_object()(r);
      typename Kernel::Point_3 po = 
        kernel.construct_translated_point_3_object()(ORIGIN, v);
      return CGAL::opposite(bis.oriented_side(po));
    }
  };
  
  Compare_distance_at_point_2 compare_distance_at_point_2_object() const
  {
    return Compare_distance_at_point_2();
  }
  
  class Construct_point_on_x_monotone_2
  {
  public:
    Point_2 operator() (const X_monotone_curve_2& xcurve) const
    {
      Kernel kernel;

      typename Kernel::Construct_point_3 cons_point_3 =
        kernel.construct_point_3_object();
      typename Kernel::Construct_ray_3 cons_ray_3 =
        kernel.construct_ray_3_object();
      typename Kernel::Construct_vector_3 cons_vector_3 =
        kernel.construct_vector_3_object();
      typename Kernel::Construct_direction_3 cons_dir_3 =
        kernel.construct_direction_3_object();

      Point_2 s = xcurve.source();
      Point_2 t = xcurve.target();
      
      Point_3 o = cons_point_3(ORIGIN);

      Ray_3 ray_s = cons_ray_3(o, s);
      Vector_3 vec_s = cons_vector_3(ray_s);
      Direction_3 dir_s = cons_dir_3(vec_s);

      Ray_3 ray_t = cons_ray_3(o, t);
      Vector_3 vec_t = cons_vector_3(ray_t);
      Direction_3 dir_t = cons_dir_3(vec_t);

      Ray_3 ray_norm = cons_ray_3(o, xcurve.normal());
      Vector_3 vec_norm = cons_vector_3(ray_norm);


      if (kernel.equal_3_object()(dir_s,
                                  kernel.construct_opposite_direction_3_object()
                                  (dir_t)))
      {
        Vector_3 vec_res = 
          kernel.construct_cross_product_vector_3_object()
          (vec_t, vec_norm);
        return kernel.construct_direction_3_object()(vec_res);
      }

      Vector_3 vec_res = vec_s + vec_t;

      Vector_3 cross = 
        kernel.construct_cross_product_vector_3_object()(vec_s, vec_t);
      Direction_3 dir_cross = cons_dir_3(cross);
      
      if (kernel.equal_3_object()(dir_cross, xcurve.normal()) == false)
        vec_res = kernel.construct_opposite_vector_3_object()(vec_res);

      return kernel.construct_direction_3_object()(vec_res);
    }
  };

  Construct_point_on_x_monotone_2 construct_point_on_x_monotone_2_object() const
  {
    return Construct_point_on_x_monotone_2();
  }

  class Compare_dominance_2
  {
  public:
    Comparison_result operator()(const Site_2& h1,
                                 const Site_2& h2) const
      
    {
      CGAL_envelope_voronoi_assertion(h1 == h2);
      return EQUAL;
    }
  };
  
  Compare_dominance_2 compare_dominance_2_object() const
  {
    return Compare_dominance_2();
  }
 
  class Compare_distance_above_2
  {
  private:
    typedef Spherical_voronoi_diagram_traits_2<Kernel>     Traits;
    
    /*! The base traits instance (in case it has state) */
    const Traits * m_traits;
    
  public:

    /*! Constructor
     * \param traits the traits instance
     */
    Compare_distance_above_2(const Traits * traits) : m_traits(traits) {}
    
    Comparison_result operator()(const Site_2& h1,
                                 const Site_2& h2,
                                 const X_monotone_curve_2& cv) const
    {
      typedef typename Kernel::Direction_2  Direction_2;

      Kernel kernel;

      Comparison_result res = m_traits->compare_y(h1, h2);
      if (res != EQUAL)
        return CGAL::opposite(res);

//       // Continue here & check the same place on the power diagram
//       // curve is not vertical - d1 != d2
//       const Direction_2 d1 = m_traits->project_xy(h1);
//       const Direction_2 d2 = m_traits->project_xy(h2);
//       CGAL_envelope_voronoi_assertion(kernel.equal_2_object()(d1, d2) != true);

//       typename Base::Direction_3 normal = cv.normal();

//       Direction_2 p = (cv.is_directed_right()) ?
//         Direction_2(-(normal.dy()), normal.dx()) :
//         Direction_2(normal.dy(), -(normal.dx()));
  
//       return (kernel.counterclockwise_in_between_2_object()(p, d1, d2)) ?
//         SMALLER : LARGER;

      // Each point is contained inside its area of dominance.
      // (Do not forget to check if the curve is directed right.)

      CGAL_envelope_voronoi_assertion(cv.is_vertical() == true);

      Point_3 o = kernel.construct_point_3_object()(ORIGIN);
      Plane_3 cv_plane = kernel.construct_plane_3_object()(o, cv.normal());
      
      typename Kernel::Construct_translated_point_3   translate = 
        kernel.construct_translated_point_3_object();
      res = kernel.oriented_side_3_object()(cv_plane, 
                                            translate(o, h1.vector()));
      if (cv.is_directed_right())
        return CGAL::opposite(res);
      return res;
    }
  };

  Compare_distance_above_2 compare_distance_above_2_object() const
  {
    return Compare_distance_above_2(this);
  }

  class Construct_bisector_2
  {
  public:

    template <class OutputIterator>
      OutputIterator operator()(const Site_2& s1,
                                const Site_2& s2,
                                OutputIterator o) const
    {
      // this is a temporary implementation that will change when the
      // traits class would enable to build Curve_2 from a plane.
      typename Kernel::Plane_3 bis ((typename Kernel::Plane_3)
                                    bisector(s1, s2));
      
#if defined(CGAL_FULL_X_MONOTONE_GEODESIC_ARC_ON_SPHERE_IS_SUPPORTED)
      Curve_2 c(bis);
#else
      // we use 4 different directions to make sure that there is no 
      // arc bigger then 180 degrees.
      typedef typename Kernel::Vector_3                Vector_3;
      typedef typename Kernel::Direction_3             Direction_3;

      Kernel kernel;
      Vector_3 v1 = kernel.construct_base_vector_3_object() (bis, 1);
      Vector_3 v2 = kernel.construct_base_vector_3_object() (bis, 2);

      Direction_3 d1 = kernel.construct_direction_3_object() (v1);
      Direction_3 d2 = kernel.construct_direction_3_object() (v2);
      Direction_3 d3 = kernel.construct_opposite_direction_3_object() (d1);
      Direction_3 d4 = kernel.construct_opposite_direction_3_object() (d2);

      Curve_2 c1 (d1, d2);
      Curve_2 c2 (d2, d3);
      Curve_2 c3 (d3, d4);
      Curve_2 c4 (d4, d1);
#endif

      // Temporary!!! Initiating new traits locally.
      typedef std::list<CGAL::Object>                  Object_list;
      Object_list x_monotones;

      Base base;

#if defined(CGAL_FULL_X_MONOTONE_GEODESIC_ARC_ON_SPHERE_IS_SUPPORTED)
      base.make_x_monotone_2_object()(c, std::back_inserter(x_monotones));
#else
      base.make_x_monotone_2_object()(c1, std::back_inserter(x_monotones));
      base.make_x_monotone_2_object()(c2, std::back_inserter(x_monotones));
      base.make_x_monotone_2_object()(c3, std::back_inserter(x_monotones));
      base.make_x_monotone_2_object()(c4, std::back_inserter(x_monotones));
#endif

      for (Object_list::iterator it = x_monotones.begin(); 
           it != x_monotones.end();
           ++it)
      {
        X_monotone_curve_2 curve;
        if (assign(curve, *it))
        {
          *o++ = curve;
          continue;
        }
        // all have to be x-mono curves.
        CGAL_envelope_voronoi_assertion(false);
      }
      
      return o;
    }
  };

  Construct_bisector_2 construct_bisector_2_object() const
  {
    return Construct_bisector_2();
  }

};

} //namespace CGAL

#endif