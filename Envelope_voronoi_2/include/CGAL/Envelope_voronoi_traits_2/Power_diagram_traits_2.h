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
// Author(s): Ophir Setter          <ophir.setter@post.tau.ac.il>
//

/*! \file   Power_diagram_traits_2.h
 * \brief  A traits class for computing the power diagram of rational circles
 * in the plane using the divide-and-conquer algorithm.
 * The coordinates of the centers of the circles must be rational and the
 * square of the radii must also be rational (in order to use rational
 * arithmetic only).
 * The traits class supports also circles with radius 0 (points) which can
 * result in the computation of the regular Voronoi diagram of points.
 */

#ifndef CGAL_POWER_DIAGRAM_TRAITS_2_H
#define CGAL_POWER_DIAGRAM_TRAITS_2_H

#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>

namespace CGAL {

template <typename Kernel_>
class Power_diagram_traits_2 : public Arr_linear_traits_2<Kernel_> {
  using Kernel = Kernel_;
  using Self = Power_diagram_traits_2<Kernel>;
  using Base = Arr_linear_traits_2<Kernel>;
  using FT = typename Kernel::FT;
  using Plane_3 = typename Kernel::Plane_3;
  using Line_2 = typename Kernel::Line_2;
  using Kernel_circle_2 = typename Kernel::Circle_2;
  using Kernel_point_2 = typename Kernel::Point_2;

public:
  using Point_2 = typename Base::Point_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
  using Site_2 = Kernel_circle_2;

  class Construct_bisector_2 {
  public:

    /*! Construct the bisector of two sites.
     * The function constructs the bisector of two circles as defined by the
     * power metric.
     * \return An output-iterator with "value-type" of X_monotone_curve_2.
     */
    template <typename OutputIterator>
      OutputIterator operator()(const Site_2& s1, const Site_2& s2,
                                OutputIterator o) const {
      // Given two circles c_1 and c_2 with (x_1, y_1) and (x_2, y_2)
      // as the respective centers and r_1 and r_2 as the respective radii,
      // the bisector of the circles is:
      // -2(x_1 - x_2)x -2(y_1 - y_2)y +
      // (x_1^2 - x_2^2) + (y_1^2 - y_2^2) - (r_1^2 - r_2^2) = 0

      //! \todo See below...
      // I am too lazy to support all kernels. There may be a way to do
      // this by constructing the line perpendicular to the segment between
      // the centers that goes throught a point $p$.
      // The point $p$ should be the translated vector between the two center
      // with some value that depends on the radii of the circles.
      // I am using the "illegal" constructor of CGAL::Line_2 object...

      Kernel kernel;
      auto cons_cen_2 = kernel.construct_center_2_object();
      auto comp_x_2 = kernel.compute_x_2_object();
      auto comp_y_2 = kernel.compute_y_2_object();
      auto comp_squa_rad_2 = kernel.compute_squared_radius_2_object();

      Kernel_point_2 c_1 = cons_cen_2(s1);
      Kernel_point_2 c_2 = cons_cen_2(s2);

      // If the centers are equal then there is no bisector curve.
      if (kernel.equal_2_object()(c_1, c_2) == true) return o;

      FT x_1 = comp_x_2(c_1);
      FT y_1 = comp_y_2(c_1);

      FT x_2 = comp_x_2(c_2);
      FT y_2 = comp_y_2(c_2);

      FT sqr_r_1 = comp_squa_rad_2(s1);
      FT sqr_r_2 = comp_squa_rad_2(s2);

      const FT minus_two(-2);
      FT a = minus_two*(x_1 - x_2);
      FT b = minus_two*(y_1 - y_2);
      FT c = square(x_1) - square(x_2) + square(y_1) - square(y_2) -
        sqr_r_1 + sqr_r_2;

      typename Kernel::Line_2 l(a, b, c);

      CGAL_envelope_voronoi_assertion_code
        (typename Kernel::Line_2 bis =
         kernel.construct_bisector_2_object()(c_1, c_2));
      CGAL_envelope_voronoi_assertion_code
        (typename Kernel::Line_2 opp_l =
         kernel.construct_opposite_line_2_object()(l));
      CGAL_envelope_voronoi_assertion((compare(sqr_r_1, sqr_r_2) != EQUAL) ||
                                      kernel.equal_2_object()(l, bis) ||
                                      kernel.equal_2_object()(opp_l, bis));

      *o++ = X_monotone_curve_2(l);
      return o;
    }
  };

  Construct_bisector_2 construct_bisector_2_object() const
  { return Construct_bisector_2(); }

  class Compare_distance_above_2 {
  protected:
    const Self* m_traits;

  public:
    Compare_distance_above_2(const Self * traits) : m_traits(traits) {}

    Comparison_result operator()( const Site_2& h1, const Site_2& h2,
                                  const X_monotone_curve_2& cv) const {
      // The radii of the circles does not affect the slope of the bisector,
      // therefore, it also does not affect which of the sites dominates which
      // side of the bisector. We use a geometric test on the centers of the
      // circles to determine that.

      Kernel ker;
      typename Kernel::Construct_center_2 cons_cen_2 =
        ker.construct_center_2_object();
      Kernel_point_2 c_1 = cons_cen_2(h1);
      Kernel_point_2 c_2 = cons_cen_2(h2);
      Comparison_result res = ker.compare_y_2_object()(c_2, c_1);
      if (res == EQUAL) res = ker.compare_x_2_object()(c_1, c_2);
      CGAL_envelope_voronoi_assertion_code
        (Kernel_point_2 p = construct_above_point(cv));
      CGAL_envelope_voronoi_assertion
        (compare(m_traits->distance(p, h1), m_traits->distance(p, h2)) == res);
      return res;
    }

  protected:
    /*! The function constructs a point in the region above the given
     * curve. "Above" means to the left of the curve when going from the
     * lexicographically smaller endpoint to the lexicographically larger
     * endpoint. This function is used only for validation.
     * \param cv The input curve.
     * \return A point in the region above the given curve.
    */
    Kernel_point_2 construct_above_point(const X_monotone_curve_2& cv) {
      Kernel ker;
      Kernel_point_2 p1, p2;
      if (cv.is_segment()) {
        p1 = ker.construct_point_on_2_object()(cv.segment(), 0);
        p2 = ker.construct_point_on_2_object()(cv.segment(), 1);
      }
      else if(cv.is_ray()) {
        p1 = ker.construct_point_on_2_object()(cv.ray(), 0);
        p2 = ker.construct_point_on_2_object()(cv.ray(), 1);
      }
      else {
        p1 = ker.construct_point_on_2_object()(cv.line(), 0);
        p2 = ker.construct_point_on_2_object()(cv.line(), 1);
      }

      typename Kernel::Vector_2 vec;
      Kernel_point_2 q;
      if (ker.compare_xy_2_object()(p1, p2) == SMALLER) {
        vec = ker.construct_vector_2_object()(p1, p2);
        q = p1;
      }
      else {
        vec = ker.construct_vector_2_object()(p2, p1);
        q = p2;
      }

      typename Kernel::Vector_2 res_vec =
        ker.construct_perpendicular_vector_2_object()(vec, COUNTERCLOCKWISE);

      return ker.construct_translated_point_2_object()(q, res_vec);
    }
  };

  Compare_distance_above_2 compare_distance_above_2_object() const
  { return Compare_distance_above_2(this); }

  class Compare_distance_at_point_2 {
  public:
    Compare_distance_at_point_2(const Self * traits) : m_traits(traits) {}

    Comparison_result operator()(const Site_2& h1, const Site_2& h2,
                                 const Point_2& p) const {
      return CGAL::compare(m_traits->distance(p, h1), m_traits->distance(p, h2));
    }

  protected:
    const Self* m_traits;
  };

  Compare_distance_at_point_2 compare_distance_at_point_2_object() const {
    return Compare_distance_at_point_2(this);
  }

  class Construct_point_on_x_monotone_2 {
  public:

    /*! Constructs an interior point on an x-monotone curve.
     * \param xcurve The x-monotone curve.
     * \return An interior point of the x-monotone curve.
     */
    Point_2 operator()(const X_monotone_curve_2& xcurve) const {
      Kernel ker;
      if(xcurve.is_segment())
        return ker.construct_midpoint_2_object()(xcurve.left(), xcurve.right());
      else if(xcurve.is_ray())
        return ker.construct_point_on_2_object()(xcurve.ray(), 1);

      CGAL_envelope_voronoi_assertion(xcurve.is_line());
      return ker.construct_point_on_2_object()(xcurve.line(), 0);
    }
  };

  Construct_point_on_x_monotone_2 construct_point_on_x_monotone_2_object() const
  { return Construct_point_on_x_monotone_2(); }

  class Compare_dominance_2 {
  public:

    Compare_dominance_2(const Self* traits) : m_traits(traits) {}

    /*! Compares the dominance relation between two sites with no bisector.
     * Compares which of the two sites dominates the whole plane in case both
     * sites do not have a bisector. Two sites have no bisector only if they
     * have the same center.
     * \param h1 The first site (circle).
     * \param h2 The second site (circle).
     * \pre Both circles must have the same center.
     * \return
     */
    Comparison_result operator()(const Site_2& h1, const Site_2& h2) const {
      Kernel ker;
      CGAL_envelope_voronoi_precondition_code
        (Kernel_point_2 c_1 = ker.construct_center_2_object()(h1););
      CGAL_envelope_voronoi_precondition_code
        (Kernel_point_2 c_2 = ker.construct_center_2_object()(h2););

      CGAL_envelope_voronoi_precondition(ker.equal_2_object()(c_1, c_2));

      // The circle with the larger radii dominates.
      auto comp_sqr_rad_2 = ker.compute_squared_radius_2_object();
      FT sqr_r_1 = comp_sqr_rad_2(h1);
      FT sqr_r_2 = comp_sqr_rad_2(h2);

      // The circle that dominates the plane is the circle with the bigger
      // radii. This means that points on the plane are CLOSER to it (this is
      // the reason we return the opposite result).
      Comparison_result res = compare(sqr_r_2, sqr_r_1);
      CGAL_envelope_voronoi_assertion_code
        (typename Kernel::Point_2 org = ker.construct_point_2_object()(ORIGIN));

      CGAL_envelope_voronoi_assertion
        (compare(m_traits->distance(org, h1), m_traits->distance(org, h2)) ==
         res);

      return res;
    }

  protected:
    const Self* m_traits;
  };

  Compare_dominance_2 compare_dominance_2_object() const
  { return Compare_dominance_2(this); }

protected:
  /*! Computes the "distance" from a point p to a circle s.
   * Computes the "power distance" from a point p to a circle s.
   * The power distance from the point (x, y) to the circle (c, r) is
   * defined to be:
   * (x - c_x)^2 + (y - c_y)^2 - r^2
   *
   * we return the value without the x^2 and y^2 because it will be the same
   * for all circle sites of the diagram (this function is used only in the
   * functor above).
   * \param p The point in $R^2$
   * \param s The circle.
   * \return The "power distance".
   */
  FT distance(const Point_2& p, const Site_2& s) const {
    Kernel ker;
    typename Kernel::Compute_x_2 comp_x_2 = ker.compute_x_2_object();
    typename Kernel::Compute_y_2 comp_y_2 = ker.compute_y_2_object();

    Kernel_point_2 c = ker.construct_center_2_object()(s);
    FT sqr_r = ker.compute_squared_radius_2_object()(s);
    FT c_x = comp_x_2(c);
    FT c_y = comp_y_2(c);
    FT x = comp_x_2(p);
    FT y = comp_y_2(p);
    const FT minus_two(-2);
    FT res =
      minus_two*c_x*x + minus_two*c_y*y + square(c_x) + square(c_y) - sqr_r;
    return res;
  }
};

} //namespace CGAL

#endif // CGAL_POWER_DIAGRAM_TRAITS_2_H
