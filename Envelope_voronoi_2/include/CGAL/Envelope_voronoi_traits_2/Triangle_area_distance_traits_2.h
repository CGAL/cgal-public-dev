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

/*! \file Triangle_area_distance_traits_2.h
 * \brief This file contains a traits class for computing the VD defined by the
 * triangle area 2-point distance function
 * See: G. Barequet et al. - 2-Point site Voronoi diagrams,
 * Disc. Appl. Math. 122 (2002)
 * \todo Change the sites to be Kernel::Segment_2 and not
 * std::pair<Point_2, Point_2>.
 */

#ifndef CGAL_TRIANGLE_AREA_DISTANCE_TRAITS_2_H
#define CGAL_TRIANGLE_AREA_DISTANCE_TRAITS_2_H

#include <CGAL/number_utils.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>

namespace CGAL {

template <typename Kernel_>
class Triangle_area_distance_traits_2 : public Arr_linear_traits_2<Kernel_> {
  using Kernel = Kernel_;
  using Base = Arr_linear_traits_2<Kernel>;
  using Self = Triangle_area_distance_traits_2<Kernel>;
  using FT = typename Kernel::FT;
  using Plane_3 = typename Kernel::Plane_3;
  using Line_2 = typename Kernel::Line_2;
  using Ray_2 = typename Kernel::Ray_2;

public:
  using Point_2 = typename Kernel::Point_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
  using Site_2 = std::pair<Point_2, Point_2>;

  Triangle_area_distance_traits_2() {}

  Triangle_area_distance_traits_2(const Kernel& ker) : Base(ker) {}

  class Construct_bisector_2 {
  public:
    Construct_bisector_2(const Self* traits) : m_traits(traits) {}

    /*! Constructs the bisector between two sites.
     * The function constructs returns a range of x-monotone curve that
     * together represets the bisector of two site, pq and rs.
     * The site pq consists of the two points $p$ and $q$ and the site
     * rs consits of the two points $r$ and $s$.
     * In the general case the bisector of two sites consists of two
     * rational lines.
     * This function splits the two lines into 4 arcs (if needed) and
     * makes them x-monotone curves.
     * The geometric function where the lines are constructed is the
     * bisector funtion in Triangle_area_distance_traits_2.
     */
    template <typename OutputIterator>
      OutputIterator operator()(const Site_2& pq, const Site_2& rs,
                                OutputIterator o) const {
      typename Kernel::Line_2 l1, l2;
      std::size_t n = m_traits->bisector(pq, rs, l1, l2);
      if (n == 0) return o;
      if (n == 1) {
        *o++ = X_monotone_curve_2(l1);
        return o;
      }

      CGAL_envelope_voronoi_assertion(n == 2);

      // If the two lines intersect, we have to make them 4 rays instead
      // of 2 lines (the envelope code requires it).
      CGAL::Object obj = m_traits->Kernel::intersect_2_object()(l1, l2);
      if (obj.is_empty()) {
        // parallel lines.
        *o++ = X_monotone_curve_2(l1);
        *o++ = X_monotone_curve_2(l2);
      }
      else {
        auto const_ray = m_traits->construct_ray_2_object();
        const Point_2* inter = object_cast<Point_2>(&obj);
        *o++ = X_monotone_curve_2(const_ray(*inter, l1));
        *o++ = X_monotone_curve_2(const_ray(*inter, l1.opposite()));
        *o++ = X_monotone_curve_2(const_ray(*inter, l2));
        *o++ = X_monotone_curve_2(const_ray(*inter, l2.opposite()));
      }

      return o;
    }

  protected:
    const Self* m_traits;
  };

  Construct_bisector_2 construct_bisector_2_object() const
  { return Construct_bisector_2(this); }

  class Compare_distance_above_2 {
  public:
    Compare_distance_above_2(const Self* traits) : m_traits(traits) {}

    /*! Checks to which of the sites dominates the area above the curve.
     * The function determines which of the sites dominates the area above
     * The x-monotone curve.
     * The bisector lines partition the plane into a maximum of 4~regions.
     * If one of the points of a Voronoi site is inside a specific region,
     * then the region is dominated by the Voronoi site of the point.
     * When crossing a bisector we move from a region that is dominated by
     * one site to a region that is dominated by the other site
     * (except for the case where all four points are collinear).
     * We compute the number of curves we need to cross to get from one of
     * the points of the first site and from the region above the input curve
     * to the upper-most face. (If we have a vertical line then we take the
     * left one.)
     * If both numbers are even or odd together we return that the region
     * above the input curve is dominated by the first site. If one is odd
     * and the other is even then the region is dominated by the second site.
     */
    Comparison_result operator()(const Site_2& h1, const Site_2& h2,
                                 const X_monotone_curve_2& cv) const {
      CGAL_envelope_voronoi_assertion(h1.first != h2.first ||
                                      h1.second != h2.second);
      CGAL_envelope_voronoi_assertion(h1.first != h2.second ||
                                      h1.second != h2.first);

      // If the 4 points are collinear, then the area is dominated by the
      // pair of points that creates the segment with the smaller length.
      //! \todo Maybe there is a faster way to check this...
      if (m_traits->collinear_2_object()(h1.first, h1.second, h2.first) &&
          m_traits->collinear_2_object()(h1.first, h1.second, h2.second)) {
        // The smaller segment dominates.
        auto comp_sq = m_traits->compute_squared_distance_2_object();
        return CGAL::compare(comp_sq(h1.first, h1.second),
                             comp_sq(h2.first, h2.second));
      }

      // We need to take one of the points of the first site which is not
      // a point of the second site.
      Point_2 p = h1.first;

      // We construct the two bisector line, and determin on where cv
      // is. Then we consider all possible cases.
      Line_2 l1, l2;
      CGAL_envelope_voronoi_assertion_code(size_t n = )
        m_traits->bisector(h1, h2, l1, l2);
      CGAL_envelope_voronoi_assertion(n != 0);

      // we need a point on the x-curve that is not the intersection point
      // of l1 and l2.
      typename Kernel::Construct_point_on_2 cons_on =
        m_traits->construct_point_on_2_object();
      Point_2 p_on_xcv = m_traits->construct_point_on_x_monotone_2_object()(cv);

      // now we insepect the location of p_on_xcv relative to the lines,
      // and the location of p relative to the line. When passing a bisector
      // the sign that this function returns should change.
      CGAL::Comparison_result p_l1 = compare_point_and_line(p, l1);
      CGAL::Comparison_result p_l2;
      if (p_l1 == EQUAL) {
        const Point_2& q = h1.second;

        p_l1 = compare_point_and_line(q, l1);
        p_l2 = compare_point_and_line(q, l2);
      }
      else
        p_l2 = compare_point_and_line(p, l2);

      CGAL::Comparison_result p_xcv_l1 = compare_point_and_line(p_on_xcv, l1);
      CGAL::Comparison_result p_xcv_l2 = compare_point_and_line(p_on_xcv, l2);

      CGAL_envelope_voronoi_assertion(p_l1 != EQUAL && p_l2 != EQUAL);
      CGAL_envelope_voronoi_assertion(p_xcv_l1 == EQUAL || p_xcv_l2 == EQUAL);
      CGAL_envelope_voronoi_assertion(p_xcv_l1 != EQUAL || p_xcv_l2 != EQUAL);

      if (p_xcv_l1 == EQUAL) return CGAL::opposite(p_xcv_l2 * p_l1 * p_l2);

      // else
      return CGAL::opposite(p_xcv_l1 * p_l1 * p_l2);

      // Line_2 sup_line;
      // if (cv.is_segment())
      //   sup_line = cv.segment().supporting_line();
      // if (cv.is_ray())
      //   sup_line = cv.ray().supporting_line();
      // else
      //   sup_line = cv.line();

      // CGAL::Comparison_result on_l1 = POSITIVE;
      // if (sup_line == l2 || sup_line == l2.opposite()) on_l1 = NEGATIVE;

      // // test for intersection.
      // CGAL::Object obj = CGAL::intersection(l1, l2);
      // if (obj.is_empty()) { // parallel lines.
      //   CGAL::Comparison_result l1_above_l2 = POSITIVE;
      //   CGAL::Comparison_result pq_is_top = POSITIVE;

      //   // vertical lines
      //   if (l1.is_vertical()) {
      //     l1_above_l2 =
      //       CGAL::opposite(CGAL::compare(l1.x_at_y(0), l2.x_at_y(0)));
      //     if (p.x() < q.x()) {
      //       if (p.x() > r.x() || p.x() > s.x())
      //         pq_is_top = NEGATIVE;
      //     }
      //     else if (q.x() > r.x() || q.x() > s.x()) pq_is_top = NEGATIVE;
      //   }
      //   else {
      //     l1_above_l2 = CGAL::compare(l1.y_at_x(0), l2.y_at_x(0));
      //     if (p.y() > q.y()) {
      //       if (p.y() < r.y() || p.y() < s.y()) pq_is_top = NEGATIVE;
      //     }
      //     else if (q.y() < r.y() || q.y() < s.y()) pq_is_top = NEGATIVE;
      //   }

      //   if (l1_above_l2 != EQUAL)
      //     return CGAL::opposite(on_l1 * l1_above_l2 * pq_is_top);

      //   // else
      //   return pq_is_top;
      // }
      // else {
      //   CGAL_envelope_voronoi_assertion(cv.is_line() == false);

      //   // not parallel lines.
      //   const Point_2 *inter = object_cast<Point_2> (&obj);
      //   Point_2 p_on_cv;
      //   if (cv.is_segment())
      //     p_on_cv = CGAL::midpoint(cv.segment().vertex(0),
      //                              cv.segment().vertex(1));
      //   else p_on_cv = cv.ray().point(1);

      //   CGAL::Comparison_result cv_on_top_of_line = POSITIVE;
      //   if (on_l1 == POSITIVE) {
      //     if (l2.is_vertical()) {
      //     }
      //   }
      // }
    }

  protected:
    const Self* m_traits;

    CGAL::Comparison_result
    compare_point_and_line(const Point_2& p, const Line_2& l) const {
      if (l.is_vertical()) {
        return CGAL::opposite(m_traits->compare_x_at_y_2_object()(p, l));
        // Original code used:
        // CGAL::compare(p.x(), l.x_at_y(p.y()));
      }
      else {
        return m_traits->compare_y_at_x_2_object()(p, l);
        // Original code used:
        // CGAL::compare(p.y(), l.y_at_x(p.x()));
        }
      }
  };

  Compare_distance_above_2 compare_distance_above_2_object() const
  { return Compare_distance_above_2(this); }

  class Compare_distance_at_point_2 {
  public:
    Compare_distance_at_point_2(const Self* traits) : m_traits(traits) {}

    Comparison_result operator()(const Site_2& h1, const Site_2& h2,
                                 const Point_2& p) const {
      return CGAL::compare(_distance(p, h1), _distance(p, h2));
    }

  private:
    /*! The function computes the distance between a site and a planar point.
     * In this traits class the distance is the area of the triangle created
     * by the given point and the points of the site.
     * We actually return the double area of the triangle (and not the area
     * itself) because there is no point in dividing the number.
     * \param x The point in the plane.
     * \param s A Voronoi site composed of 2 planar ponits.
     *
     * \return
     */
    FT _distance(const Point_2& x, const Site_2& s) const {
      return CGAL::abs(
        m_traits->compute_area_2_object()(s.first, s.second, x));

      //// Maybe the implementation below can be faster because it does
      //// not divide the result by two.
      //       const Point_2 &p = s.first;
      //       const Point_2 &q = s.second;

      //       FT A = p.y() - q.y();
      //       FT B = q.x() - p.x();
      //       FT C = p.x()*q.y() - p.y()*q.x();

      //       return CGAL::abs(A*x.x() + B*x.y() + C);
    }

  protected:
    const Self* m_traits;
  };

  Compare_distance_at_point_2 compare_distance_at_point_2_object() const
  { return Compare_distance_at_point_2(this); }

  class Construct_point_on_x_monotone_2 {
  public:
    //! The function constructs a point in the interior of an x-mono curve
    Point_2 operator()(const X_monotone_curve_2& xcurve) const {
      Kernel ker;
      if(xcurve.is_segment())
        return ker.construct_midpoint_2_object()(xcurve.left(), xcurve.right());
      else if(xcurve.is_ray())
        return ker.construct_point_on_2_object()(xcurve.ray(), 1);

      CGAL_envelope_voronoi_assertion(xcurve.is_line());
      return ker.construct_point_on_2_object()(xcurve.line(), 1);
    }
  };

  Construct_point_on_x_monotone_2 construct_point_on_x_monotone_2_object() const
  { return Construct_point_on_x_monotone_2(); }

  class Compare_dominance_2 {
  public:
    Comparison_result operator()(const Site_2& h1, const Site_2& h2) const {
      // only if the 4 points are collinear and they create two equal segment
      // there is no bisector.
      // We can always return equal.
      CGAL_envelope_voronoi_assertion_code(Kernel k;);
      CGAL_envelope_voronoi_assertion
        (k.compute_squared_distance_2_object()(h1.first, h1.second) ==
         k.compute_squared_distance_2_object()(h2.first, h2.second));
      CGAL_envelope_voronoi_assertion
        (k.collinear_2_object()(h1.first, h1.second, h2.first) &&
        k.collinear_2_object()(h1.first, h1.second, h2.second));

      return EQUAL;
    }
  };

  Compare_dominance_2 compare_dominance_2_object() const
  { return Compare_dominance_2(); }

public:
  /*! Returns the bisector of two Voronoi site.
   * The function returns the bisector of two Voronoi sites under the
   * triangle-area distance function. In the general case, the bisector
   * is composed of two rational lines.
   * \param pq The first site, containing the two points p and q.
   * \param rs The second site, containing the two points r and s.
   * \param l1 The first line of the bisector.
   * \param l2 The second line of the bisector.
   *
   * \return The number of lines consisting the bisector. If this is a double
   *         line, then l1 and l2 are equal.
   */
  std::size_t bisector(const Site_2& pq, const Site_2& rs,
                       Line_2& l1, Line_2 &l2) const {
    /* Triangle area is defined by a determinent of a matrix:
     * http://mathworld.wolfram.com/TriangleArea.html
     * which results in the following distance function:
     *
     *   2d(x, {p, q}) = \abs{(p_y - q_y)x + (q_x - p_x)y + p_x*q_y - p_y*q_x}
     *
     *  which, in turn results with the following two lines as bisectors.
     */
    const Point_2& p = pq.first;
    const Point_2& q = pq.second;

    const Point_2& r = rs.first;
    const Point_2& s = rs.second;

    FT A_pq = p.y() - q.y();
    FT B_pq = q.x() - p.x();
    FT C_pq = p.x()*q.y() - p.y()*q.x();

    FT A_rs = r.y() - s.y();
    FT B_rs = s.x() - r.x();
    FT C_rs = r.x()*s.y() - r.y()*s.x();

    FT A1 = A_pq - A_rs;
    FT B1 = B_pq - B_rs;
    FT C1 = C_pq - C_rs;

    if (A1 == 0 && B1 == 0) return 0;

    l1 = Line_2(A1, B1, C1);
    l2 = Line_2(A_pq + A_rs, B_pq + B_rs, C_pq + C_rs);

    CGAL_envelope_voronoi_assertion(l1.is_degenerate() == false);
    CGAL_envelope_voronoi_assertion(l2.is_degenerate() == false);

      //! \todo Is there another way to know this?
    if ((l1 == l2) || (l1 == l2.opposite())) return 1;
    return 2;
  }
};

} //namespace CGAL

#endif // CGAL_TRIANGLE_AREA_DISTANCE_TRAITS_2_H
