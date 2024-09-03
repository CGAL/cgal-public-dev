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

/*! \file Algebraic_curve_envelope_traits_2.h
 */

#ifndef CGAL_ENVVOR_ALGEBRAIC_CURVE_TRAITS_2_H
#define CGAL_ENVVOR_ALGEBRAIC_CURVE_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/algebraic_kernel_1_tools.h>
#include <CGAL/Envelope_voronoi_2/envelope_voronoi_assertions.h>

namespace CGAL {

/*! A first adaptable traits for points diagram.
 * \todo Change name + move to another file.
 */
// template <typename T_CurvedKernel>
// class Adaptable_points_traits {
// public:
//   using Curved_kernel_2 = T_CurvedKernel;
//   using Self = Adaptable_points_traits<Curved_kernel_2>;
//   using Algebraic_kernel_2 = typename Curved_kernel_2::Curve_kernel_2;
//   using Rational = typename Algebraic_kernel_2::Boundary;
//   using Surface_3 = std::pair<Rational, Rational>;
//   using Poly_rat_1 = CGAL::Polynomial<Rational>;
//   using Poly_rat_2 = CGAL::Polynomial<Poly_rat_1>;
//   using Poly_rat_3 = CGAL::Polynomial<Poly_rat_2>;
//   using PT_rat_3 = CGAL::Polynomial_traits_d<Poly_rat_3>;
//   using Monomial = std::pair<CGAL::Exponent_vector, Rational>;

//   struct Construct_distance_function_2 {
//     Poly_rat_3 operator() (const Surface_3& s) const {
//       std::vector<Monomial> monomials;
//       CGAL::Exponent_vector ev(3);

//       ev[0]=2; ev[1]=0; ev[2]=0; monomials.push_back(Monomial(ev,
//                                                               Rational(1)));
//       ev[0]=0; ev[1]=2; ev[2]=0; monomials.push_back(Monomial(ev,
//                                                               Rational(1)));
// //    ev[0]=0; ev[1]=0; ev[2]=2; monomials.push_back(Monomial(ev, c));

// //    ev[0]=1; ev[1]=1; ev[2]=0; monomials.push_back(Monomial(ev,d));
// //    ev[0]=1; ev[1]=0; ev[2]=1; monomials.push_back(Monomial(ev,e));
// //    ev[0]=0; ev[1]=1; ev[2]=1; monomials.push_back(Monomial(ev,f));

//       ev[0]=1; ev[1]=0; ev[2]=0; monomials.push_back(Monomial(ev,
//                                                               -2 * s.first));
//       ev[0]=0; ev[1]=1; ev[2]=0; monomials.push_back(Monomial(ev,
//                                                               -2 * s.second));
//       ev[0]=0; ev[1]=0; ev[2]=1; monomials.push_back(Monomial(ev,
//                                                               Rational(-1)));

//       Rational t = s.first*s.first + s.second*s.second;
//       ev[0]=0; ev[1]=0; ev[2]=0; monomials.push_back(Monomial(ev, t));

//       typename PT_rat_3::Construct_polynomial construct_polynomial;
//       return construct_polynomial(monomials.begin(), monomials.end());
//     }
//   };

//   Construct_distance_function_2 construct_distance_function_2_object() const
//   { return Construct_distance_function_2(); }
// };

template <typename AdaptableTraits>
class Algebraic_curve_envelope_traits_2 :
    public AdaptableTraits::Curved_kernel_2 {
public:
  using Adaptable_traits_2 = AdaptableTraits;
  using Curved_kernel_2 = typename Adaptable_traits_2::Curved_kernel_2;
  using Algebraic_kernel_2 = typename Curved_kernel_2::Curve_kernel_2;

  using Self = Algebraic_curve_envelope_traits_2< Adaptable_traits_2 >;

  using Rational = typename Algebraic_kernel_2::Boundary;
  using Point_2 = typename Curved_kernel_2::Point_2;
  using Curve_2 = typename Curved_kernel_2::Curve_2;
  using X_monotone_curve_2 = typename Curved_kernel_2::X_monotone_curve_2;
  using Multiplicity = typename Curved_kernel_2::Multiplicity;

  using Poly_rat_1 = CGAL::Polynomial<Rational>;
  using Poly_rat_2 = CGAL::Polynomial<Poly_rat_1>;
  using Poly_rat_3 = CGAL::Polynomial<Poly_rat_2>;
  using PT_rat_1 = CGAL::Polynomial_traits_d<Poly_rat_1>;
  using PT_rat_2 = CGAL::Polynomial_traits_d<Poly_rat_2>;
  using PT_rat_3 = CGAL::Polynomial_traits_d<Poly_rat_3>;

  using Integer = typename Fraction_traits<Rational>::Numerator_type;
  using Poly_int_1 = CGAL::Polynomial<Integer>;
  using Poly_int_2 = CGAL::Polynomial<Poly_int_1>;
  using Poly_int_3 = CGAL::Polynomial<Poly_int_2>;
  using PT_int_1 = CGAL::Polynomial_traits_d<Poly_int_1>;
  using PT_int_2 = CGAL::Polynomial_traits_d<Poly_int_2>;
  using PT_int_3 = CGAL::Polynomial_traits_d<Poly_int_3>;
  using Monomial = std::pair<CGAL::Exponent_vector, Rational>;

  using Obj_list = std::list<CGAL::Object>;

  using Surface_3 = typename Adaptable_traits_2::Site_2;

  // I think that I meant that the surfaces are of the form z = f(x, y).
  using Xy_monotone_surface_3 = Poly_int_3;
  using Intersetion_cache_key =
    std::pair<Xy_monotone_surface_3, Xy_monotone_surface_3>;
  using Intersetion_cache_value = std::list<CGAL::Object>;
  using Intersetion_cache =
    std::map<Intersetion_cache_key, Intersetion_cache_value >;

  /*! Cache for trisector in order not to create trisectors once. */
  mutable Intersetion_cache m_inter_cache;

  //! \todo Should take the state from the base class.
  Adaptable_traits_2*       m_base_traits;

public:
  Algebraic_curve_envelope_traits_2< Adaptable_traits_2 >
  (Adaptable_traits_2 *base_traits)
  : m_base_traits(base_traits) {}

  Intersetion_cache& intersetion_cache() const { return m_inter_cache; }

  /*! Creates xy-monotone surfaces from a general surface
   *  return a past-the-end iterator.
   */
  class Make_xy_monotone_3 {
  protected:
    const Self* m_traits;

  public:
    Make_xy_monotone_3(const Self* traits) : m_traits(traits) {}

    // create xy-monotone surfaces from a general surface
    // return a past-the-end iterator
    template <typename OutputIterator>
    OutputIterator operator()(const Surface_3& s, bool is_lower,
                              OutputIterator o) const {
      Poly_rat_3 rat_surf = m_traits->m_base_traits->
        construct_distance_function_2_object()(s);

      Poly_int_3 surf;
      Integer dummy_int;
      typename Fraction_traits<Poly_rat_3>::Decompose()
        (rat_surf, surf, dummy_int);

      *o++ = Xy_monotone_surface_3(surf);
      return o;
    }
  };

  /*! Get a Make_xy_monotone_3 functor object. */
  Make_xy_monotone_3
  make_xy_monotone_3_object() const
  { return Make_xy_monotone_3(this); }

  /*! Constructs the projected boundary of the surface. In our case
    We deal with unbounded surfaces. */
  class Construct_projected_boundary_2 {
  public:
    // insert into the OutputIterator all the (2d) curves of the boundary of
    // the vertical projection of the surface on the xy-plane
    // the OutputIterator value type is X_monotone_curve_2
    template <typename OutputIterator>
    OutputIterator
    operator()(const Xy_monotone_surface_3& s, OutputIterator o) const
    { return o; }
  };

  /*! Get a Construct_projected_boundary_2 functor object. */
  Construct_projected_boundary_2
  construct_projected_boundary_2_object() const
  { return Construct_projected_boundary_2(); }

  class Construct_projected_intersections_2 {
  protected:
    const Self* m_self;

  public:
    Construct_projected_intersections_2(const Self* self) : m_self(self) {}

    // insert into OutputIterator all the (2d) projections on the xy plane of
    // the intersection objects between the 2 surfaces
    // the data type of OutputIterator is Object
    template <typename OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const {
      CGAL_envelope_voronoi_precondition(s1 != s2);

      // We want to find the algebraic curve which is the intersection
      // between s1 and s2 (which are both polynomials with 3 variables.
      return o;
    }
  };

  /*! Get a Construct_projected_intersections_2 functor object. */
  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const
  { return Construct_projected_intersections_2(this); }

  class Compare_z_at_xy_3 {

  protected:
    const Self* m_self;

  public:
    Compare_z_at_xy_3(const Self* self) : m_self(self) {}

  public:
    // check which of the surfaces is closer to the envelope at the xy
    // coordinates of point (i.e. lower if computing the lower envelope, or
    // upper if computing the upper envelope)
    // precondition: the surfaces are defined in point
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      CGAL_error();
      return EQUAL;
    }

    // check which of the surfaces is closer to the envelope at the xy
    // coordinates of cv (i.e. lower if computing the lower envelope, or upper
    // if computing the upper envelope)
    // precondition: the surfaces are defined in all points of cv, and the
    //               answer is the same for each of these points
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const
    {
      CGAL_error();
      return EQUAL;
    }

    // check which of the surfaces is closer to the envelope.
    // (i.e. lower if computing the lower envelope, or upper
    // if computing the upper envelope)
    // precondition: there is no intersections between the surfaces.
    Comparison_result operator()(const Xy_monotone_surface_3& s1,
                                 const Xy_monotone_surface_3& s2) const {
      CGAL_error();
      return EQUAL;
    }
  };

  /*! Get a Compare_z_at_xy_3 functor object. */
  Compare_z_at_xy_3
  compare_z_at_xy_3_object() const { return Compare_z_at_xy_3(this); }

  class Compare_z_at_xy_above_3 {
  protected:
    const Self* m_traits;

  public:
    Compare_z_at_xy_above_3(const Self* traits) : m_traits(traits) {}

    // check which of the surfaces is closer to the envelope on the points above
    // the curve cv (i.e. lower if computing the lower envelope, or upper if
    // computing the upper envelope)
    // precondition: the surfaces are defined above cv
    //               the choise between s1 and s2 for the envelope is the same
    //               for every point in the infinitesimal region above cv
    //               the surfaces are EQUAL over the curve cv
    Comparison_result
    operator()(const X_monotone_curve_2& cv,
               const Xy_monotone_surface_3& s1,
               const Xy_monotone_surface_3& s2) const {
      std::pair<Rational, Rational> p_above =
        rational_point_above_or_below_arc(cv, true);
      Comparison_result res = compare_surfaces_at(s1, s2, p_above);
      return res;
    }
  };

  /*! Get a Compare_z_at_xy_above_3 functor object. */
  Compare_z_at_xy_above_3
  compare_z_at_xy_above_3_object() const
  { return Compare_z_at_xy_above_3(this); }

  class Compare_z_at_xy_below_3 {
  protected:
    const Self* m_traits;

  public:
    Compare_z_at_xy_below_3(const Self* traits) : m_traits(traits) {}

    Comparison_result
    operator()(const X_monotone_curve_2& cv,
               const Xy_monotone_surface_3& s1,
               const Xy_monotone_surface_3& s2) const
    {
      std::pair<Rational, Rational> p_below =
        rational_point_above_or_below_arc(cv, false);
      Comparison_result res = compare_surfaces_at(s1, s2, p_below);
      return res;
    }
  };

  /*! Get a Compare_z_at_xy_below_3 functor object. */
  Compare_z_at_xy_below_3
  compare_z_at_xy_below_3_object() const
  { return Compare_z_at_xy_below_3(this); }

protected:
  static Comparison_result compare_surfaces_at(const Poly_int_3& s1,
                                               const Poly_int_3& s2,
                                               const std::pair<Rational,
                                               Rational>& point) {
    std::vector<Poly_rat_1> rat_polys;
    typename PT_rat_1::Construct_polynomial construct_poly;
    rat_polys.push_back(construct_poly(point.first));
    rat_polys.push_back(construct_poly(point.second));
    rat_polys.push_back(construct_poly(Rational(0), Rational(1)));

    typename PT_int_3::Substitute substitue;
    Poly_rat_1 rat_poly1 = substitue(s1, rat_polys.begin(), rat_polys.end());
    Poly_rat_1 rat_poly2 = substitue(s2, rat_polys.begin(), rat_polys.end());

    Poly_int_1 poly1, poly2;
    Integer dummy_int;
    typename Fraction_traits<Poly_rat_1>::Decompose()
      (rat_poly1, poly1, dummy_int);
    typename Fraction_traits<Poly_rat_1>::Decompose()
      (rat_poly2, poly2, dummy_int);

    //! \todo ker should be taken from the algebraic kernel.
    typename Curved_kernel_2::Curve_kernel_2::Algebraic_kernel_1 ker;
    return CGAL::compare_smallest_nonnegative_roots(ker, poly1, poly2);
  }

//! Construct a rational point above or below an algebraic curve x-monotone arc.
/*! Construct a rational point above or below an algebraic curve x-monotone arc.
  "Above" is defined to be the area to the area to the left of the curve when
  traveling on it from a lexicographic smaller point to a lexicographic larger
  point. The function handles vertical arcs.
  \param arc X-monotone arc
  \param above True if the function should construct an above point. False
  if the function should construct a below point.
  \return A rational pair (representing x and y respectively) which forms a
  rational point above/below the given x-monotone arc.
*/
  static std::pair<Rational, Rational>
    rational_point_above_or_below_arc(
      const X_monotone_curve_2& arc,
      bool above)
  {
    using Curve_analysis_2 = typename Curved_kernel_2::Curve_2;
    using size_type = typename Curve_analysis_2::size_type;
    using Status_line_1 = typename Curve_analysis_2::Status_line_1;

    //! \todo Again, take state!!!
    Algebraic_kernel_2 ker;
    auto upper_boundary_y_2 = ker.upper_boundary_y_2_object();
    autolower_boundary_y_2 = ker.lower_boundary_y_2_object();

    if (arc.is_vertical()) {
      const Curve_analysis_2& ca = arc.curve();

      // for counting arc numbers. a vertical line should be an event (at least
      // I think so).
      Status_line_1 supporting_status_line =
        ca.status_line_for_x(arc.x());
      CGAL_envelope_voronoi_assertion(supporting_status_line.is_event());

      // We use the supporting status line of the curve to count the number of
      // arcs on the interval to the left which are below the lower endpoint of
      // the arc.
      // (Using first for the arcs to the left of each event and
      // second for the arcs to the right of each event.)
      const std::pair<size_type, size_type>& minus_inf_branches =
        supporting_status_line.number_of_branches_approaching_minus_infinity();
      size_type num_of_arcs = above ? minus_inf_branches.first :
        minus_inf_branches.second;
      if (arc.is_finite(ARR_MIN_END)) {
        const typename Point_2::Coordinate_2 min_coord =
          arc.curve_end(ARR_MIN_END).xy();
        auto comp_xy_2 = ker.compare_xy_2_object();
        for (size_type i = 0; i < supporting_status_line.number_of_events(); ++i)
        {
          if (comp_xy_2(min_coord, supporting_status_line.xy_coordinate_2(i))
              == SMALLER)
            break;

          // Using first for the arcs to the left of each event and
          // second for the arcs to the right of each event.
          num_of_arcs += above ?
            supporting_status_line.number_of_incident_branches(i).first :
            supporting_status_line.number_of_incident_branches(i).second;
        }
      }

      // going to the interval left of the vertical arc in case of above and
      // to the interval on the right in case of below.
      size_type interval_id = above ?
        supporting_status_line.index() :
        supporting_status_line.index() + 1;

      Rational x_rat = ca.boundary_value_in_interval(interval_id);
      Status_line_1 interval_status_line =
        arc.curve().status_line_at_exact_x(x_rat);

      CGAL_envelope_voronoi_assertion (num_of_arcs >= 0 && num_of_arcs <   \
                      interval_status_line.number_of_events());
      CGAL_envelope_voronoi_assertion(interval_status_line.is_event() == false);

      Rational y_rat(0);
      if (interval_status_line.number_of_events() == 0)
        return make_pair(x_rat, y_rat);

      if (num_of_arcs == 0)
        y_rat =
          lower_boundary_y_2(nterval_status_line.xy_coordinate_2(num_of_arcs));
      else
        y_rat = above ?
          upper_boundary_y_2(
            interval_status_line.xy_coordinate_2(num_of_arcs-1)) :
          lower_boundary_y_2(
            interval_status_line.xy_coordinate_2(num_of_arcs-1)) ;

      return make_pair(x_rat, y_rat);
    }
    else {
      // this is the x value of the point
      Rational x_rat = arc.boundary_in_x_range_interior();

      // this is the status line at the x value of the point
      Status_line_1 status_line = arc.curve().status_line_at_exact_x(x_rat);

      // this is the upper boundary of the y value of the arc at x
      Rational y_rat = above ?
        upper_boundary_y_2(status_line.xy_coordinate_2(arc.arcno())) :
        lower_boundary_y_2(status_line.xy_coordinate_2(arc.arcno()));

      return make_pair(x_rat, y_rat);
    }
  }
};

} //namespace CGAL

#endif // CGAL_VDL3_SINGLE_CELL_ENVELOPE_TRAITS_3_H
