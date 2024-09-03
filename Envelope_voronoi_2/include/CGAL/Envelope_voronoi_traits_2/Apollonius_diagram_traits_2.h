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
// $Id:
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>

/*! \file Apollonius_diagram_traits_2.h
 * The file contains an implementation of the Apollonius diagram based on our
 * extension for the conic traits. THIS IS AN EARLY VERSION. Currenbly the
 * traits class based on the CKvA is more stable and should be used.
 */

#ifndef CGAL_APOLLON_DIAGRAM_TRAITS_2_H
#define CGAL_APOLLON_DIAGRAM_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/functions_on_enums.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Arr_tracing_traits_2.h>
#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>

namespace CGAL {

#ifdef DEBUG_CONIC_TRAITS
template <typename Rat_kernel_, typename Alg_kernel_, typename Nt_traits_>
class Apollonius_diagram_traits_2 :
  public CGAL::Arr_tracing_traits_2<Arr_conic_traits_2<Rat_kernel_, Alg_kernel_,
                                                       Nt_traits_>>
#else
template <typename Rat_kernel_, typename Alg_kernel_, typename Nt_traits_>
class Apollonius_diagram_traits_2 :
    public Arr_conic_traits_2<Rat_kernel_, Alg_kernel_, Nt_traits_>
#endif
{
public:
  using Rat_kernel = Rat_kernel_;
  using Alg_kernel = Alg_kernel_;
  using Nt_traits = Nt_traits_;
  using Rational = typename Rat_kernel::FT;
  using Algebraic = typename Alg_kernel::FT;
  using Base = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
  using Self = Apollonius_diagram_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
  using Multiplicity = unsigned int;
  using Rat_point_2 = typename Rat_kernel::Point_2;
  using Rat_line_2 = typename Rat_kernel::Line_2;
  using Rat_vector_2 = typename Rat_kernel::Vector_2;
  using Point_2 = typename Base::Point_2;
  using Curve_2 = typename Base::Curve_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
  using Intersection_curve = std::pair<X_monotone_curve_2,    Multiplicity>;
  using Has_boundary_category = typename Base::Has_boundary_category;

  class _Apollonius_site {
    friend std::ostream & operator<<(std::ostream & os,
                                     const _Apollonius_site& site)
    { return os << site.center() << " " << site.r(); }

    friend std::istream & operator>>(std::istream & is,
                                     _Apollonius_site& site) {
      Rational x, y;
      is >> x >> y >> site._r;
      site._center = Rat_point_2(x, y);
      return is;
    }

  protected:
    Rat_point_2 _center;
    Rational _r;

  public:
    _Apollonius_site() {}

  _Apollonius_site(const Rat_point_2& center, const Rational& r) :
    _center(center), _r(r)
    {}

    const Rat_point_2& center() const { return _center; }

    const Rational& r() const { return _r; }

    const Algebraic distance(const Point_2& p) const {
      Nt_traits nt_traits;
      Algebraic x = p.x() - this->center().x();
      Algebraic y = p.y() - this->center().y();
      return nt_traits.sqrt(x*x + y*y) - this->r();
    }
  };

  using Site_2 = _Apollonius_site;

  class Compare_distance_at_point_2 {
  public:
    Comparison_result operator()(const Site_2& h1, const Site_2& h2,
                                 const Point_2& p) const
    { return CGAL::compare(h1.distance(p), h2.distance(p)); }
  };

  Compare_distance_at_point_2 compare_distance_at_point_2_object() const
  { return Compare_distance_at_point_2(); }

  class Construct_point_on_x_monotone_2 {
  public:
    Point_2 operator()(const X_monotone_curve_2& xcurve) const {
      // this is not a good implementation. We need a point in the internal
      // of the x-monotone curve.
      CGAL_warning("Need a better implementation.");
      return xcurve.source();
    }
  };

  Construct_point_on_x_monotone_2 construct_point_on_x_monotone_2_object() const
  { return Construct_point_on_x_monotone_2(); }

  class Compare_dominance_2 {
  public:
    Comparison_result operator()(const Site_2& h1, const Site_2& h2) const {
      Rational r1 = h1.r();
      Rational r2 = h2.r();
      return CGAL::compare(r2, r1);
    }
  };

  Compare_dominance_2 compare_dominance_2_object() const
  { return Compare_dominance_2(); }

  class Compare_distance_above_2 {
  public:
    Comparison_result operator()(const Site_2& h1, const Site_2& h2,
                                 const X_monotone_curve_2& cv) const {
      // We take the smaller circle (which its center is in the x-range
      // of cv) and check if it is below or above cv.
      Rat_kernel ker;

      CGAL::Sign s = CGAL::compare(h1.r(), h2.r());
      if (s == ZERO) {
        Comparison_result res =
          ker.compare_y_2_object()(h2.center(), h1.center());
        if (res == EQUAL)
          return ker.compare_x_2_object()(h2.center(), h1.center());
        return res;
      }

      Self self;
      Comparison_result res = LARGER;
      if (cv.is_facing_up()) res = SMALLER;
      CGAL_envelope_voronoi_assertion_msg(res != EQUAL,
                                          "The center can be on the bisector.");
      return res * s;
    }
  };

  Compare_distance_above_2 compare_distance_above_2_object() const
  { return Compare_distance_above_2(); }

  class Construct_bisector_2 {
  public:
    template <typename OutputIterator>
      OutputIterator operator()(const Site_2& s1, const Site_2& s2,
                                OutputIterator o) const {
      Rat_kernel ker;
      const Rational xi = ker.compute_x_2_object()(s1.center());
      const Rational yi = ker.compute_y_2_object()(s1.center());
      const Rational ri = s1.r();

      const Rational xj = ker.compute_x_2_object()(s2.center());
      const Rational yj = ker.compute_y_2_object()(s2.center());
      const Rational rj = s2.r();

      Nt_traits nt_traits;

      // if one site is completly inside the other site return empty
      // object.
      const Rational xmx = xi - xj;
      const Rational ymy = yi - yj;
      const Rational rmr = ri - rj;

      const Rational xx_2 = xmx * xmx;
      const Rational yy_2 = ymy * ymy;
      const Rational RS = rmr * rmr;

      if (nt_traits.convert(xx_2 + yy_2) <= RS) return o;

      // The bisector is a hyperbolic arc with:
      //  r = 4 [(xi-xj)^2 - RS]
      //  s = same as r (exchange x and y)
      //  t = 8 (xi - xj) (yi - yj)
      //  u = 4 [-(yi + yj)(yi - yj)(xi - xj) - (xi + xj)(xi - xj)^2 + RS(xi + xj)]
      //  v = same as u (exchange x and y)
      //  w := (xi^2 - xj^2)^2 + (yi^2 - yj^2)^2 + (xi^2 - yj^2)^2 +
      //	     (yi^2 - xj^2)^2 - (xi^2 - yi^2)^2 - (xj^2 - yj^2)^2
      // 	      -2RS (xi^2 + yi^2 + xj^2 + yj^2) + RS^2
      // where:
      //  RS = (ri-rj)^2
      const Rational r = 4 * (xx_2 - RS);
      const Rational s = 4 * (yy_2 - RS);

      const Rational t = 8 * xmx * ymy;

      const Rational xpx = xi + xj;
      const Rational ypy = yi + yj;

      const Rational u = 4 * (RS*xpx - ypy*ymy*xmx - xpx*xmx*xmx);
      const Rational v = 4 * (RS*ypy - xpx*xmx*ymy - ypy*ymy*ymy);

      const Rational xi_2 = xi * xi;
      const Rational xj_2 = xj * xj;
      const Rational yi_2 = yi * yi;
      const Rational yj_2 = yj * yj;

      const Rational xi2xj2 = (xi_2 - xj_2);
      const Rational yi2yj2 = (yi_2 - yj_2);
      const Rational xi2yj2 = (xi_2 - yj_2);
      const Rational yi2xj2 = (yi_2 - xj_2);
      const Rational xi2yi2 = (xi_2 - yi_2);
      const Rational xj2yj2 = (xj_2 - yj_2);
      const Rational w = xi2xj2 * xi2xj2 + yi2yj2 * yi2yj2
        + xi2yj2 * xi2yj2 + yi2xj2 * yi2xj2
        - xi2yi2 * xi2yi2 - xj2yj2 * xj2yj2
        - 2 * RS * (xi_2 + yi_2 + xj_2 + yj_2)
        + RS*RS;

      // In the mean time we have to have two (different) points on the curve.
      // So we intersect the hyperbola / line with a line parallel to the major
      // axis of the hyperbola / perpendicular to the line
      // \todo: change this to work with normal constructor.
      Rat_point_2 p_small = s1.center();
      Rat_point_2 p_big = s2.center();
      if (ri > rj) {
        p_small = s2.center();
        p_big = s1.center();
      }

      // We use to approximated points. one is above the center of the circle,
      // and one is below.
      Rat_line_2 axis = ker.construct_line_2_object()(p_small, p_big);
      Rat_line_2 perpend =
        ker.construct_perpendicular_line_2_object()(axis, p_small);
      Rat_vector_2 vec = ker.construct_vector_2_object () (perpend);
      Rat_vector_2 opp_vec = ker.construct_opposite_vector_2_object()(vec);

      // these are the approximation points.
      Rat_point_2 app_point1 =
        ker.construct_translated_point_2_object()(p_small, vec);
      Rat_point_2 app_point2 =
        ker.construct_translated_point_2_object()(p_small, opp_vec);

      // and these are their lines (which will intersect the hyperbola).
      Rat_line_2 inter_line1 =
        ker.construct_perpendicular_line_2_object()(perpend, app_point1);
      Rat_line_2 inter_line2 =
        ker.construct_perpendicular_line_2_object()(perpend, app_point2);

      // The coefficients of the lines:
      Rational a1 = ker.compute_a_2_object()(inter_line1);
      Rational b1 = ker.compute_b_2_object()(inter_line1);
      Rational c1 = ker.compute_c_2_object()(inter_line1);

      Rational a2 = ker.compute_a_2_object()(inter_line2);
      Rational b2 = ker.compute_b_2_object()(inter_line2);
      Rational c2 = ker.compute_c_2_object()(inter_line2);

      // compute orientation in case that this is an hyperbola.
      Orientation orient = CGAL::orientation(app_point1, p_big, app_point2);

      Algebraic x1 = nt_traits.convert(ker.compute_x_2_object()(app_point1));
      Algebraic y1 = nt_traits.convert(ker.compute_y_2_object()(app_point1));
      Algebraic x2 = nt_traits.convert(ker.compute_x_2_object()(app_point2));
      Algebraic y2 = nt_traits.convert(ker.compute_y_2_object()(app_point2));

      // if this is a line, we need it COLLINEAR
      Curve_2 c;
      if (ri != rj)
        c = Curve_2(r, s, t, u, v, w, orient,
                    Point_2(x1, y1), 0, 0, 0, a1, b1, c1,
                    Point_2(x2, y2), 0, 0, 0, a2, b2, c2, true, true);
      else
        c = Curve_2(r, s, t, u, v, w, CGAL::COLLINEAR,
                     Point_2(x1, y1), 0, 0, 0, a1, b1, c1,
                     Point_2(x2, y2), 0, 0, 0, a2, b2, c2, true, true);

      Self self;
      std::list<Object> x_mono_list;
      self.make_x_monotone_2_object()(c, std::back_inserter(x_mono_list));

      for (auto it = x_mono_list.begin(); it != x_mono_list.end(); ++it) {
        X_monotone_curve_2 x;
        if (assign(x, *it)) *o++ = x;
      }
      return o;
    }
  };

  Construct_bisector_2 construct_bisector_2_object() const
  { return Construct_bisector_2(); }
};

} //namespace CGAL

#endif
