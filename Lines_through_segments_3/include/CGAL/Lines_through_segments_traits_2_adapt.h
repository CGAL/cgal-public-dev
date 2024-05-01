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

#ifndef LINES_THROUGH_SEGMENTS_TRAITS_2_ADAPT_H
#define LINES_THROUGH_SEGMENTS_TRAITS_2_ADAPT_H

/*! \file
 *
 * The file contains different classes to adapt between the different
 * arrangement traits classes.
 * Currently, Conic and algebraic traits classes are supported.
 */
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Lines_through_segments_traits_3.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Lines_through_segments_general_functions.h>

#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#if USE_SQRT_TRAITS
#include <CGAL/Algebraic_kernel_2_1.h>
#include <CGAL/Lazy_exact_nt.h>
#endif
#include <CGAL/Root_of_traits.h>


namespace CGAL {

template <typename Traits_3_>
class Lines_through_segments_traits_on_plane_adapt {
public:
  using Traits_3 = Traits_3_;

private:
  using Rational_kernel = typename Traits_3::Rational_kernel;
  using Rational = typename Rational_kernel::FT;
  using Alg_kernel = typename Traits_3::Alg_kernel;
  using Integer = CORE::BigInt;

  using Traits_arr_on_plane_2 = typename Traits_3::Traits_arr_on_plane_2;
  using Rational_segment_2 = typename Rational_kernel::Segment_2;
  using Rational_segment_3 = typename Rational_kernel::Segment_3;
  using Rational_point_2 = typename Rational_kernel::Point_2;

  // typedef typename Traits_arr_on_plane_2::Point_2 Point_2;

  /* Specific typedefs for conic arc traits. */
  using Nt_traits = CGAL::CORE_algebraic_number_traits;
  using Conic_traits_arr_on_plane_2_ =
    CGAL::Arr_conic_traits_2<Rational_kernel, Alg_kernel, Nt_traits>;
  using Conic_traits_arr_on_plane_2 =
    CGAL::Arr_consolidated_curve_data_traits_2<Conic_traits_arr_on_plane_2_,
                                               const Rational_segment_3*>;

  using Conic_point_2 = typename Conic_traits_arr_on_plane_2::Point_2;
  using Conic_curve_2_wd = typename Conic_traits_arr_on_plane_2::Curve_2;
  using Conic_curve_2 = typename Conic_traits_arr_on_plane_2_::Curve_2;
  using Conic_x_monotone_curve_2 =
    typename Conic_traits_arr_on_plane_2::X_monotone_curve_2;

  /* Specific typedefs for Rational arc traits . */
#if USE_SQRT_TRAITS
#if USE_LAZY
  using AK1 = CGAL::Algebraic_kernel_2_1<Rational >;
#else
  using AK1 = CGAL::Algebraic_kernel_2_1<Rational >;
#endif
  using Sqrt_ext = typename  AK1::Algebraic_real_1;
#else
  using AK1 = CGAL::Algebraic_kernel_d_1<Integer>;
#endif

  using Rational_arc_traits_arr_on_plane_2_ =
    CGAL::Arr_rational_function_traits_2<AK1>;
  using Rational_arc_traits_arr_on_plane_2 =
    CGAL::Arr_consolidated_curve_data_traits_2
    <Rational_arc_traits_arr_on_plane_2_,const Rational_segment_3*>;

  using Rational_arc_point_2 =
    typename Rational_arc_traits_arr_on_plane_2::Point_2;
  using Rational_arc_x_monotone_curve_2 =
    typename Rational_arc_traits_arr_on_plane_2::X_monotone_curve_2;
  using Rational_arc_curve_2 =
    typename Rational_arc_traits_arr_on_plane_2_::Curve_2;
  using Rational_arc_curve_2_wd =
    typename Rational_arc_traits_arr_on_plane_2::Curve_2;

  using Polynomial = CORE::Polynomial<Rational>;
  using CORE_big_float = typename CORE::BigFloat;
  using CORE_big_float_interval = typename CORE::BFInterval;

  using Traits_algebraic_kernel =
    typename Rational_arc_traits_arr_on_plane_2::Algebraic_kernel_d_1;
  using Polynomial_1 =
    typename Rational_arc_traits_arr_on_plane_2::Polynomial_1;
  using Bound = typename Traits_algebraic_kernel::Bound;
  using Bound_pair = std::pair<Bound,Bound>;

public:
  using Algebraic_real_1 =
    typename Rational_arc_traits_arr_on_plane_2::Algebraic_real_1;
  using Algebraic = typename Alg_kernel::FT;

  /**************************************************************
   * The following function return the orientation of a curve.
   *
   * Output:
   *      CGAL::COLLINEAR or CGAL::COUNTERCLOCKWISE or CGAL::CLOCKWISE
   ***************************************************************/

  //! \brief
  CGAL::Orientation orientation(const Conic_curve_2_wd& curve,
                                const Alg_kernel* ker)
  { return this->orientation_conic_traits(curve,ker); }

  //! \brief
  CGAL::Orientation orientation(const Conic_x_monotone_curve_2& curve,
                                const Alg_kernel* ker)
  { return this->orientation_conic_traits(curve,ker); }


  //! \brief
  template <typename Curve_2>
  CGAL::Orientation orientation_conic_traits(const Curve_2& curve,
                                             const Alg_kernel* ker)
  { return curve.orientation(); }

  //! \brief
  CGAL::Orientation orientation(const Rational_arc_x_monotone_curve_2& curve,
                                const Alg_kernel* ker)
  { return this->orientation_rational_arc_traits(curve,ker); }

  //! \brief
  CGAL::Orientation orientation(const Rational_arc_curve_2_wd& curve,
                                const Alg_kernel* ker)
  { return this->orientation_rational_arc_traits(curve,ker); }

  //! \brief
  template <typename Curve_2>
  CGAL::Orientation orientation_rational_arc_traits(const Curve_2& curve,
                                                    const Alg_kernel* ker) {
    Algebraic left_x;
    Algebraic right_x;
    Algebraic left_y;
    Algebraic right_y;

    convert_point (curve.left(),left_x, left_y);
    convert_point (curve.right(),right_x, right_y);

    Algebraic mid_x = left_x + (right_x - left_x)/2;
    Algebraic mid_y = left_y + (right_y - left_y)/2;

    typename Alg_kernel::Point_2 mid_p(mid_x, mid_y);
    typename Alg_kernel::Point_2 right_p(right_x, right_y);
    typename Alg_kernel::Point_2 left_p(left_x, left_y);
    CGAL::Orientation orient = ker->orientation_2_object()(left_p, mid_p,
                                                           right_p);
    if (orient == LEFT_TURN) return CGAL::COUNTERCLOCKWISE;
    else if (orient == RIGHT_TURN) return CGAL::CLOCKWISE;
    return CGAL::COLLINEAR;
  }

  /**************************************************************
   * The following functions returns the middle point of a curve.
   *
   * Output:
   *      Point_2
   ***************************************************************/

  //! \brief
  template <typename Point_2>
  void get_mid_point(const Conic_curve_2_wd& curve, Point_2& output_p)
  { return get_mid_point_conic_traits(curve,output_p); }

  //! \brief
  template <typename Point_2>
  void get_mid_point(const Conic_x_monotone_curve_2& curve,
                     Point_2& output_p)
  { return get_mid_point_conic_traits(curve,output_p); }

  //! \brief
  template <typename Curve_2,typename Point_2>
  void get_mid_point_conic_traits(const Curve_2& curve, Point_2& output_p) {
    Point_2 source = curve.source();
    Point_2 target = curve.target();

    Algebraic mid_x = source.x() + (target.x()- source.x())/2;

    Rational t = curve.t();
    Rational u = curve.u();
    Rational v = curve.v();
    Rational w = curve.w();

    if (t == 0 && v == 0)
      output_p = Point_2(mid_x, ((-u) * mid_x - w));
    else
      output_p = Point_2(mid_x, ((-u) * mid_x - w)/(t * mid_x + v));
  }

  //! \brief
  template <typename Point_2>
  void get_mid_point(const Rational_arc_curve_2_wd& curve, Point_2& output_p)
  { return get_mid_point_rational_arc_traits(curve,output_p); }

  //! \brief
  template <typename Point_2>
  void get_mid_point(const Rational_arc_x_monotone_curve_2& curve,
                     Point_2& output_p)
  { return get_mid_point_rational_arc_traits(curve,output_p); }

  //! \brief
  template <typename Curve_2,typename Point_2>
  void get_mid_point_rational_arc_traits(const Curve_2& curve,
                                         Point_2& output_p) {
    Point_2 source = curve.source();
    Point_2 target = curve.target();

    Algebraic mid_x = source.x() + (target.x()- source.x())/2;

#if 1
    Algebraic y_numer = curve.numerator().evaluate(mid_x);
    Algebraic y_denom = curve.denominator().evaluate(mid_x);
#else
    Polynomial core_numer = convert_polynomial(curve.numerator());
    Polynomial core_denom = convert_polynomial(curve.denominator());

    Algebraic y_numer = core_numer.eval(mid_x);
    Algebraic y_denom = core_denom.eval(mid_x);
#endif

    Algebraic mid_y(y_numer / y_denom);

    output_p = Point_2(mid_x,mid_y);
  }

  /**************************************************************
   * The following functions returns the y val, given a x val on curve.
   *
   * Output:
   *      Point_2
   ***************************************************************/

  //! \brief
  template <typename Number_type>
  Algebraic get_y_val(const Conic_curve_2_wd& curve, const Number_type& x)
  { return get_y_val_conic_traits(curve,x); }

  //! \brief
  template <typename Number_type>
  Algebraic get_y_val(const Conic_x_monotone_curve_2& curve,
                      const Number_type& x)
  { return get_y_val_conic_traits(curve,x); }

  //! \brief
  Rational get_y_val_ratioal(const Conic_x_monotone_curve_2& curve,
                             const Rational& x)
  { return (((-curve.u()) * x - curve.w())/(curve.t() * x + curve.v())); }

  //! \brief
  Rational get_y_val_ratioal(const Rational_arc_x_monotone_curve_2& curve,
                             const Rational& x) {
    CGAL_error_msg("Not implemented.");
    return 0;
  }

  template <typename Curve_2,typename Number_type>
  Algebraic get_y_val_conic_traits(const Curve_2& curve,
                                   const Number_type& x) {
    Rational t = curve.t();
    Rational u = curve.u();
    Rational v = curve.v();
    Rational w = curve.w();

    /* Get the y value of curve at x */
    Algebraic temp_y;
    if (t == 0 && v == 0) temp_y = ((-u) * x - w);
    else temp_y = (((-u) * x - w)/(t * x + v));
    return temp_y;
  }

  //! \brief
  template <typename Number_type>
  Algebraic get_y_val(const Rational_arc_curve_2_wd& curve,
                      const Number_type& x)
  { return get_y_val_rational_arc_traits(curve,x); }

  Rational get_x_val(const Conic_curve_2_wd& curve, const Rational& y)
  { return (-curve.w() - curve.v() * y) / (curve.t() * y + curve.u()); }

  //! \brief
  Rational get_x_val_of_vertical_cruve(
     const Conic_x_monotone_curve_2& curve) {
     CGAL_assertion(is_vertical(curve));
     return (-curve.w()/curve.u());
  }

  //! \brief
  Rational get_x_val_of_vertical_cruve(
     const Rational_arc_x_monotone_curve_2& curve) {
     CGAL_error_msg("Not implemented.");
     return 0;
  }

  //! \brief
  Rational get_x_val(const Conic_x_monotone_curve_2& curve, const Rational& y)
  { return (-curve.w() - curve.v() * y) / (curve.t() * y + curve.u()); }

  //! \brief
  Rational get_x_val(const Rational_arc_x_monotone_curve_2& curve,
                     const Rational& y)
  { return get_x_val_rational_arc_traits(curve,y); }

  //! \brief
  Rational get_x_val(const Rational_arc_curve_2_wd& curve,
                     const Rational& y)
  { return get_x_val_rational_arc_traits(curve,y); }

  //! \brief
  template <typename Number_type>
  Algebraic get_y_val(const Rational_arc_x_monotone_curve_2& curve,
                      const Number_type& x)
  { return get_y_val_rational_arc_traits(curve,x); }

  //! \brief
  template <typename Curve_2>
  Algebraic get_y_val_rational_arc_traits(const Curve_2& curve,
                                          const Algebraic_real_1& _x) {
    Algebraic x = convert_real_to_algebraic(_x);
    return get_y_val_rational_arc_traits(curve,x);
  }

  //! \brief
  template <typename Curve_2>
  Algebraic get_y_val_rational_arc_traits(const Curve_2& curve,
                                          const Algebraic& x) {
    Polynomial core_numer = convert_polynomial(curve.numerator());
    Polynomial core_denom = convert_polynomial(curve.denominator());

    Algebraic y_numer = core_numer.eval(x);
    Algebraic y_denom = core_denom.eval(x);

    Algebraic y(y_numer / y_denom);
    return y;
  }

  //! \brief
  template <typename Curve_2>
  Rational get_x_val_rational_arc_traits(const Curve_2& curve,
                                         const Rational& y) {
    Polynomial core_numer = convert_polynomial(curve.numerator());
    Polynomial core_denom = convert_polynomial(curve.denominator());
    Rational x_0_num = core_numer.getCoeffi(0);
    Rational x_1_num = core_numer.getCoeffi(1);
    Rational x_0_denom = core_denom.getCoeffi(0);
    Rational x_1_denom = core_denom.getCoeffi(1);

    Rational x = (y * x_0_denom - x_0_num) / (x_1_num - x_1_denom * y);

    return x;
  }

  //! \brief
  template <typename Curve_2>
  void get_horizontal_asymptote_y_val(const Curve_2& curve, Algebraic& y) {
    Polynomial core_numer = convert_polynomial(curve.numerator());
    Polynomial core_denom = convert_polynomial(curve.denominator());
    CGAL_assertion(CGAL::degree(curve.numerator()) ==
                   CGAL::degree(curve.denominator()));

    if (CGAL::degree(curve.numerator()) == 1) {
      Integer x_1_num = core_numer.getCoeffi(1);
      Integer x_1_denom = core_denom.getCoeffi(1);
      y = Algebraic(x_1_num) / Algebraic(x_1_denom);
    }
    else {
      //if (core_numer.degree() == 0)
      CGAL_assertion(CGAL::degree(curve.denominator()) == 0);
      Integer x_0_num = core_numer.getCoeffi(0);
      Integer x_0_denom = core_denom.getCoeffi(0);
      y = Algebraic(x_0_num) / Algebraic(x_0_denom);
    }
  }

  /**************************************************************
   * The following function adapts creation of rational segment
   * on plane arr.
   *
   * Input:
   *      source - Segment_2 source point.
   *      target  - Segment_2 end point.
   *
   ***************************************************************/

  //! \brief
  void create_segment_on_plane_arr(Conic_curve_2_wd& cv,
                                   const Rational_point_2& source,
                                   const Rational_point_2& target,
                                   const Rational_segment_3* data) {
    //! \todo Pass a traits object instead of locally creating a new one.
    Conic_traits_arr_on_plane_2_ traits;
    auto ctr_cv = traits.construct_curve_2_object();
    cv = Conic_curve_2_wd(ctr_cv(Rational_segment_2(source, target)), data);
  }

  //! \brief
  void create_segment_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                   const Rational_arc_point_2& source,
                                   const Rational_arc_point_2& target,
                                   const Rational_segment_3* data)
  { cv = Rational_arc_curve_2(source, target); }

  //! \brief
  void create_segment_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                   const Rational_point_2& source,
                                   const Rational_point_2& target,
                                   const Rational_segment_3* data) {
    std::vector<Rational> P(2);
    Rational_arc_traits_arr_on_plane_2 traits_2;
    /* Create vertical segment */
    if (source.x() == target.x()) {
      Rational_arc_point_2 ps(traits_2.construct_point_2_object()(source.x(),source.y()));
      Rational_arc_point_2 pt(traits_2.construct_point_2_object()(target.x(),target.y()));
      // Vertical_segment ver_todo =
      //   traits_2.construct_vertical_segment_object()(ps, pt);
      std::cout << "TOOD create verical segment" << std::endl;
      CGAL_error_msg("");

      // cv = Rational_arc_curve_2(ps,pt); TODO TODO TODO //asafp asaf
    }
    else {
      Algebraic_real_1 xs(source.x());
      Algebraic_real_1 xt(target.x());

      P[1] = (target.y() - source.y())/(target.x() - source.x());
      P[0] = source.y() - P[1] * source.x();

      // std::cout << "Segment = (" << xs << "," << xt << ")   "
      //           << P[0] << " + x * " << P[1] << std::endl;
      cv = traits_2.construct_curve_2_object()(P.begin(), P.end(), xs, xt);
    }
  }

  /**************************************************************
   * The following function adapts creation of rational curve_2
   * on plane arr.
   * The curve is of the type:
   *
   *            V[3]*xy + V[2]*x + v[1]*y + V[0] = 0
   *
   *            y = (-v[0] -v[2]*x)/(v[1] + v[3]*x)
   *
   * Input:
   *      source - Curve source point.
   *      target - Curve end point.
   *      V[]    - Curve coefficients.
   *
   *
   ***************************************************************/

  //! \brief
  void create_curve_on_plane_arr(Conic_curve_2_wd& cv,
                                 const Rational& source_x,
                                 const Rational& target_x,
                                 const Rational coefficients[4],
                                 const Rational_segment_3* data) {
    CGAL_assertion((target_x * coefficients[3] + coefficients[1]) !=
                   Rational(0));
    CGAL_assertion((source_x * coefficients[3] + coefficients[1]) !=
                   Rational(0));

    Rational_point_2 source(source_x,
                            (source_x * (coefficients[2]) + coefficients[0])/
                            (source_x * coefficients[3] + coefficients[1]));

    Rational_point_2 target(target_x,
                            (target_x * (coefficients[2]) + coefficients[0])/
                            (target_x * coefficients[3] + coefficients[1]));

    Conic_point_2 a_source(source.x(), source.y());
    Conic_point_2 a_target(target.x(), target.y());

    create_curve_on_plane_arr(cv, a_source, a_target, coefficients, data);
  }

private:
  //! \brief
  void create_curve_on_plane_arr(Conic_curve_2_wd& cv,
                                 const Conic_point_2& a_source,
                                 const Conic_point_2& a_target,
                                 const Rational coefficients[4],
                                 const Rational_segment_3* data) {
    Algebraic mid_point_y = (a_target.y() + a_source.y())/2;
    Algebraic mid_point_x = (a_target.x() + a_source.x())/2;
    Algebraic mid_point_on_hyp_y =
      ((mid_point_x * (coefficients[2]) + (coefficients[0])) /
       (mid_point_x * coefficients[3] + coefficients[1]));
    Conic_traits_arr_on_plane_2_ traits;
    auto ctr_cv = traits.construct_curve_2_object();
    if (mid_point_y < mid_point_on_hyp_y) {
      cv = Conic_curve_2_wd(ctr_cv(Rational(0), Rational(0),
                                   coefficients[3], -coefficients[2],
                                   coefficients[1], -coefficients[0],
                                   CGAL::CLOCKWISE, a_source, a_target),
                            data);
    }
    else if (mid_point_y > mid_point_on_hyp_y) {
      cv = Conic_curve_2_wd(ctr_cv(Rational(0), Rational(0),
                                   coefficients[3], -coefficients[2],
                                   coefficients[1], -coefficients[0],
                                   CGAL::COUNTERCLOCKWISE, a_source, a_target),
                            data);
    }
    else {
      cv = Conic_curve_2_wd(ctr_cv(Rational(0), Rational(0),
                                   coefficients[3], -coefficients[2],
                                   coefficients[1], -coefficients[0],
                                   CGAL::COLLINEAR, a_source, a_target),
                            data);
    }
  }

public:
  //! \brief
  void create_curve_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                 const Rational& source_x,
                                 const Rational& target_x,
                                 const Rational coefficients[4],
                                 const Rational_segment_3* data) {
    std::vector<Rational> P2(2);
    P2[0] = coefficients[0];
    P2[1] = coefficients[2];

    std::vector<Rational> Q2(2);
    Q2[0] = coefficients[1];
    Q2[1] = coefficients[3];

    Algebraic_real_1 xs(source_x);
    Algebraic_real_1 xt(target_x);

    Rational_arc_traits_arr_on_plane_2 traits_2;
    auto ctr_cv = traits_2.construct_curve_2_object();
    cv = Rational_arc_curve_2_wd(ctr_cv(P2.begin(), P2.end(),
                                        Q2.begin(), Q2.end(), xs, xt), data);
    // cv = Rational_arc_curve_2(P2, Q2,xs,xt);
  }

  //! \brief
  template <typename NT1, typename NT2>
  void create_curve_on_plane_arr(Rational_arc_curve_2_wd& new_cv,
                                 const NT1& source_x,
                                 const NT2& target_x,
                                 const Rational_arc_curve_2& old_cv)
  { create_curve_on_plane_arr_pr(new_cv, source_x, target_x, old_cv); }

  //! \brief
  template <typename NT1,typename NT2>
  void create_curve_on_plane_arr(Rational_arc_curve_2& new_cv,
                                 const NT1& source_x,
                                 const NT2& target_x,
                                 const Rational_arc_x_monotone_curve_2& old_cv)
  { create_curve_on_plane_arr_pr(new_cv, source_x, target_x, old_cv); }

  //! \brief
  /* Curve on arr may be either x monotone or not x monotoe curve. */
  template <typename NT1,typename NT2, typename Curve_on_arr_2>
  void create_curve_on_plane_arr_pr(Rational_arc_curve_2& new_cv,
                                    const NT1& source_x,
                                    const NT2& target_x,
                                    const Curve_on_arr_2& old_cv) {
    Polynomial core_numer = convert_polynomial(old_cv.numerator());
    Polynomial core_denom = convert_polynomial(old_cv.denominator());

    std::vector<Rational> P2(2);
    P2[0] = core_numer.getCoeffi(0);
    P2[1] = core_numer.getCoeffi(1);

    std::vector<Rational> Q2(2);
    Q2[0] = core_denom.getCoeffi(0);
    Q2[1] = core_denom.getCoeffi(1);

    Algebraic_real_1 xs(source_x);
    Algebraic_real_1 xt(target_x);
    Rational_arc_traits_arr_on_plane_2 traits_2;
    auto ctr_cv = traits_2.construct_curve_2_object();
    new_cv = ctr_cv(P2.begin(), P2.end(), Q2.begin(), Q2.end(), xs, xt);
//    new_cv = Rational_arc_curve_2(P2, Q2,xs,xt);
  }

  bool is_vertical(const Rational_arc_curve_2_wd& arc) {
    // arc.source_infinite_in_x() == CGAL::ARR_INTERIOR &&//TODO vertical
    // arc.target_infinite_in_x() == CGAL::ARR_INTERIOR &&
    return (arc.source_x() == arc.target_x());
  }

  //! \brief
  bool is_vertical(const Rational_arc_x_monotone_curve_2& arc) {
    std::cout << "TODO return is vertical" << std::endl;
    CGAL_error_msg("TODO return is vertical");
    // arc.left_infinite_in_x() == CGAL::ARR_INTERIOR &&
    //TODO asafp can not find the function infinite_in_x()
    // arc.right_infinite_in_x() == CGAL::ARR_INTERIOR &&
    return (arc.source_x() == arc.target_x());
  }

  //! \brief
  bool is_vertical(const Conic_x_monotone_curve_2& curve)
  { return (curve.t() == 0 && curve.v() == 0); }

  //! \brief
  void create_curve_on_plane_arr(Rational_arc_curve_2_wd& new_cv,
                                 const Rational_arc_x_monotone_curve_2& old_cv)
  {
    typename Rational_arc_traits_arr_on_plane_2::Is_vertical_2 is_ver_obj;
    Rational_arc_traits_arr_on_plane_2 traits;
    if (is_ver_obj(old_cv)) {
       // Rational_arc_point_2 ps(traits.construct_point_2_object()(old_cv.source.x(),old_cv.source.y()));
       // Rational_arc_point_2 pt(traits.construct_point_2_object()(old_cv.target.x(),old_cv.target.y()));
//       Vertical_segment ver_todo = traits.construct_vertical_segment_object()(ps, pt);
       std::cout << "TODO create vertical line" << std::endl;
       CGAL_error_msg("TODO create vertical line");
//      new_cv = Rational_arc_curve_2(old_cv.source(), old_cv.target()); TODO TODO asafp
      return;
    }

    Polynomial core_numer = convert_polynomial(old_cv.numerator());
    Polynomial core_denom = convert_polynomial(old_cv.denominator());

    std::vector<Rational> P2(2);
    P2[0] = core_numer.getCoeffi(0);
    P2[1] = core_numer.getCoeffi(1);

    std::vector<Rational> Q2(2);
    Q2[0] = core_denom.getCoeffi(0);
    Q2[1] = core_denom.getCoeffi(1);

    Algebraic_real_1 xs(old_cv.source_x());
    Algebraic_real_1 xt(old_cv.target_x());

    auto ctr_cv = traits.construct_curve_2_object();
    new_cv = ctr_cv(P2.begin(), P2.end(), Q2.begin(), Q2.end(), xs, xt);
    // new_cv = Rational_arc_curve_2(P2, Q2, xs, xt);
  }

  //! \brief
  void create_curve_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                 const Rational& source_x,
                                 bool dir_right,
                                 const Rational coefficients[4],
                                 const Rational_segment_3* data) {
    std::vector<Rational> P2(2);
    P2[0] = coefficients[0];
    P2[1] = coefficients[2];

    std::vector<Rational> Q2(2);
    Q2[0] = coefficients[1];
    Q2[1] = coefficients[3];

    Algebraic_real_1 xs(source_x);
    Rational_arc_traits_arr_on_plane_2 traits;
    auto ctr_cv = traits.construct_curve_2_object();
    cv = Rational_arc_curve_2_wd(ctr_cv(P2.begin(), P2.end(),
                                        Q2.begin(), Q2.end(), xs, dir_right),
                                 data);
    // cv = Rational_arc_curve_2(P2, Q2,xs,dir_right);
  }

  //! \brief
  void create_curve_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                 const Rational coefficients[4],
                                 const Rational_segment_3* data) {
    std::vector<Rational> P2(2);
    P2[0] = coefficients[0];
    P2[1] = coefficients[2];

    std::vector<Rational> Q2(2);
    Q2[0] = coefficients[1];
    Q2[1] = coefficients[3];
    Rational_arc_traits_arr_on_plane_2 traits;
    auto ctr_cv = traits.construct_curve_2_object();
    cv = Rational_arc_curve_2_wd(ctr_cv(P2.begin(), P2.end(),
                                        Q2.begin(), Q2.end()), data);
    // cv = Rational_arc_curve_2(P2, Q2);
  }

  //! \brief
  void create_vertical_segment_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                            const Rational_point_2& source) {
    // Rational_arc_traits_arr_on_plane_2 traits_2;
    // Vertical_segment ver_todo = traits_2.construct_vertical_segment_object()(
    // traits_2.construct_point_2_object()(source.x(),source.y()));
    std::cout << "TOOD create vertical segment" << std::endl;
    CGAL_error_msg("");
    // cv = Rational_arc_curve_2(Rational_arc_point_2(source.x(),source.y()));TODO TODO TODO //asafp asaf
  }

  //! \brief
  void create_vertical_segment_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                            const Rational_point_2& source,
                                            bool is_dir_up) {
    // Rational_arc_traits_arr_on_plane_2 traits_2;
    // Vertical_segment ver_todo = traits_2.construct_vertical_segment_object()(
    // traits_2.construct_point_2_object()(source.x(),source.y()), is_dir_up);
     std::cout << "TOOD create vertical segment" << std::endl;
     CGAL_error_msg("");
    // cv = Rational_arc_curve_2(Rational_arc_point_2(source.x(), source.y()),
    //                           is_dir_up); TODO TODO TODO //asafp asaf
  }

  //! \brief
  void create_vertical_segment_on_plane_arr(Conic_curve_2_wd& cv,
                                            const Rational_point_2& source) {
    CGAL_error_msg("Conic arc traits do not support infinite vertical line.");
  }

  //! \brief
  void create_vertical_segment_on_plane_arr(Conic_curve_2_wd& cv,
                                            const Rational_point_2& source,
                                            bool is_dir_up) {
    CGAL_error_msg("Conic arc traits do not support infinite vertical line.");
  }

  //! \brief
  void create_horizontal_curve_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                            const Rational& S2_t,
                                            const Rational& source_x,
                                            bool dir_right,
                                            const Rational_segment_3* data) {
     create_horizontal_curve_on_plane_arr_pr(cv,S2_t,source_x,dir_right, data);
  }

  //! \brief
  void create_horizontal_curve_on_plane_arr_pr
  (Rational_arc_x_monotone_curve_2& cv,
   const Rational& S2_t,
   const Rational& source_x,
   bool dir_right,
   const Rational_segment_3* data) {
     create_horizontal_curve_on_plane_arr_pr(cv,S2_t,source_x,dir_right, data);
  }

private:
  //! \brief
  template <typename Curve_2>
  void
  create_horizontal_curve_on_plane_arr_pr(Curve_2& cv,
                                          const Rational& S2_t,
                                          const Rational& source_x,
                                          bool dir_right,
                                          const Rational_segment_3* data) {
    std::vector<Rational> P2(1);
    P2[0] = S2_t;

    std::vector<Rational> Q2(1);
    Q2[0] = Rational(1);

    Algebraic_real_1 xs(source_x);
    Rational_arc_traits_arr_on_plane_2 traits;
    auto ctr_cv = traits.construct_curve_2_object();
    cv = Curve_2(ctr_cv(P2.begin(), P2.end(), Q2.begin(), Q2.end(),
                        xs, dir_right), data);
    // cv = Rational_arc_curve_2(P2, Q2,xs,dir_right);
  }

public:
  //! \brief
  void create_horizontal_curve_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                            const Rational& S2_t,
                                            const Rational_segment_3* data)
  { create_horizontal_curve_on_plane_arr_pr(cv,S2_t, data); }

  //! \brief
  void create_horizontal_curve_on_plane_arr
  (Rational_arc_x_monotone_curve_2& cv,
   const Rational& S2_t,
   const Rational_segment_3* data)
  { create_horizontal_curve_on_plane_arr_pr(cv,S2_t, data); }

private:
  //! \brief
  template <typename Curve_2>
  void create_horizontal_curve_on_plane_arr_pr
  (Curve_2& cv,
   const Rational& S2_t,
   const Rational_segment_3* data) {
    std::vector<Rational> P2(1);
    P2[0] = S2_t;

    std::vector<Rational> Q2(1);
    Q2[0] = Rational(1);
    Rational_arc_traits_arr_on_plane_2 traits_2;
    auto ctr_cv = traits_2.construct_curve_2_object();
    cv = Curve_2(ctr_cv(P2.begin(), P2.end(), Q2.begin(), Q2.end()), data);
    // cv = Rational_arc_curve_2(P2, Q2);
  }
public:
  //! \brief
  template <typename NT1, typename NT2>
  void create_horizontal_curve_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                            const Rational& S2_t,
                                            const NT1& source_x,
                                            const NT2& target_x,
                                            const Rational_segment_3* data) {
    std::vector<Rational> P2(1);
    P2[0] = S2_t;

    std::vector<Rational> Q2(1);
    Q2[0] = Rational(1);

    Algebraic_real_1 xs(source_x);
    Algebraic_real_1 xt(target_x);
    Rational_arc_traits_arr_on_plane_2 traits_2;
    auto ctr_cv = traits_2.construct_curve_2_object();
    cv = Rational_arc_curve_2_wd(ctr_cv(P2.begin(), P2.end(),
                                        Q2.begin(), Q2.end(), xs, xt), data);
    // cv = Rational_arc_curve_2(P2, Q2, xs, xt);
  }

  //! \brief
  bool is_horizontal(const Rational_arc_curve_2_wd& curve) {
    const int deg_p (CGAL::degree(curve.numerator()));
    const int deg_q (CGAL::degree(curve.denominator()));
    if (deg_p == 0 && deg_q == 0) return true;
    return false;
  }

  //! \brief
  bool is_horizontal(const Rational_arc_x_monotone_curve_2& curve) {
     CGAL_error_msg("Not implemented.");
     return false;
  }

  //! \brief
  bool has_y_value_at_x(const Rational_arc_x_monotone_curve_2& curve,
                        const Rational& x) {
     CGAL_error_msg("Not implemented.");
     return 0;
  }

  //! \brief
  bool has_y_value_at_x(const Conic_x_monotone_curve_2& curve,
                        const Rational& x)
  { return (curve.t() * x + curve.v() != 0); }

  //! \brief
  bool has_x_value_at_y(const Rational_arc_x_monotone_curve_2& curve,
                        const Rational& y) {
     CGAL_error_msg("Not implemented.");
     return 0;
  }

  //! \brief
  bool has_x_value_at_y(const Conic_x_monotone_curve_2& curve,
                        const Rational& y)
  { return ((curve.t() * y + curve.u()) != 0); }

  //! \brief
  bool has_x_value_at_y(const Conic_curve_2_wd& curve, const Rational& y)
  { return ((curve.t() * y + curve.u()) != 0); }

  /* Horizontal line segment*/
  bool is_horizontal(const Conic_curve_2_wd& curve)
  { return (curve.t() == 0 && curve.u() == 0); }

  /* Horizontal line segment*/
  bool is_horizontal(const Conic_x_monotone_curve_2& curve)
  { return (curve.t() == 0 && curve.u() == 0); }

  //! \brief
  Rational get_y_val_of_horizontal_curve(const Conic_curve_2_wd& curve) {
    CGAL_assertion(is_horizontal(curve));
    return (-curve.w()/curve.v());
  }

  //! \brief
  Rational get_y_val_of_horizontal_curve(const Conic_x_monotone_curve_2& curve)
  {
     CGAL_assertion(is_horizontal(curve));
     return (-curve.w()/curve.v());
  }

  //! \brief
  Rational get_y_val_of_horizontal_curve
  (const Rational_arc_x_monotone_curve_2& curve) {
     CGAL_error_msg("Not implemented.");
     return 0;
  }

  /**************************************************************
   * The following function adapts creation of rational curve_2
   * on plane arr.
   *
   * Constructs a circular arc going from p1, the source, through p2, p3
   * and p4 to p5, the target (notice all points have integer coordinates).
   * Precondition: No three points of the five are not collinear.The curve is
   * of the type:
   *
   * Input:
   *
   *
   ***************************************************************/

public:
  //! \brief
  void create_curve_on_plane_arr(Conic_curve_2_wd& cv,
                                 const Rational coefficients[4],
                                 const Rational_segment_3* data) {
    CGAL_error_msg("Unbounded arcs are not supported with the conic traits");
  }

  //! \brief
  void create_curve_on_plane_arr(Conic_curve_2_wd& cv,
                                 const Rational_point_2 point_1,
                                 const Rational_point_2 point_2,
                                 const Rational_point_2 point_3,
                                 const Rational_point_2 point_4,
                                 const Rational_point_2 point_5,
                                 const Rational_segment_3* data) {
    Conic_traits_arr_on_plane_2_ traits;
    auto ctr_cv = traits.construct_curve_2_object();
     cv = Conic_curve_2_wd(ctr_cv(point_1, point_2, point_3,
                                  point_4, point_5), data);
  }

  //! \brief
  void create_curve_on_plane_arr(Rational_arc_curve_2_wd& cv,
                                 const Rational_point_2 point_1,
                                 const Rational_point_2 point_2,
                                 const Rational_point_2 point_3,
                                 const Rational_point_2 point_4,
                                 const Rational_point_2 point_5,
                                 const Rational_segment_3* data) {
    assert(false);

    // typename Rational_kernel::Conic_2   temp_conic;
    // Rational coefficients [4];

    // temp_conic.set (point_1, point_2, point_3, point_4, point_5);

    // /*
    //  * Get the conic coefficients:
    //  * rx^2 + sy^2 + txy + ux + vy + w = 0
    //  * r and s equal to 0.
    //  *
    //  * txy + ux + vy + w = 0
    //  */
    // CGAL_precondition((temp_conic.r() == 0) && (temp_conic.s() == 0));

    // coefficients[3] = temp_conic.t();
    // coefficients[2] = temp_conic.u();
    // coefficients[1] = temp_conic.v();
    // coefficients[0] = temp_conic.w();

    // create_curve_on_plane_arr(cv, point_1.x(), point_5.x(), coefficients, data);
  }

  //! \brief
  void create_horizontal_curve_on_plane_arr(Conic_curve_2_wd& cv,
                                            const Rational& S2_t,
                                            const Rational& source_x,
                                            bool dir_right) {
    CGAL_error_msg("Unbounded arcs are not supported with the conic traits");
  }

  //! \brief
  void create_horizontal_curve_on_plane_arr(Conic_curve_2_wd& cv,
                                            const Rational& S2_t,
                                            const Rational_segment_3* data) {
    CGAL_error_msg("Unbounded arcs are not supported with the conic traits");
  }

  //! \brief
  void create_horizontal_curve_on_plane_arr(Conic_curve_2_wd& cv,
                                            const Rational& S2_t,
                                            const Rational& source_x,
                                            const Rational& target_x,
                                            const Rational_segment_3* data) {
    Rational_point_2 p1(source_x, S2_t);
    Rational_point_2 p2(target_x, S2_t);
    create_segment_on_plane_arr(cv,p1,p2, data);
  }

  /**************************************************************
   * The following function adapts creation of Point_2
   * on plane arr.
   *
   * Input:
   *
   *
   *
   ***************************************************************/
  void
  construct_point(typename Conic_traits_arr_on_plane_2::Point_2& point_2,
                  const Rational& x, const Rational& y) {
    using Point_2 = typename Conic_traits_arr_on_plane_2::Point_2;
    point_2 = Point_2(x,y);
  }

  void
  construct_point(typename Rational_arc_traits_arr_on_plane_2::Point_2& point_2,
                  const Rational& x, const Rational& y) {

    Rational_arc_traits_arr_on_plane_2 traits_2;
    point_2 = traits_2.construct_point_2_object()(x,y);
  }

private:
  //! \brief
  Polynomial convert_polynomial(const Polynomial_1 poly) const {
    // const int    d = CGAL::degree(poly);
    // Rational* coeffs = new Rational[d+1];
    Polynomial core_poly;
    // for (int ii = 0 ; ii <= d; ++ii)
    // {
    //   coeffs[ii] = CGAL::get_coefficient(poly,ii);
    // }

    // core_poly = Polynomial (d,coeffs);

    return core_poly;
  }

public:
#if USE_SQRT_TRAITS
  template <class ET>
  Algebraic convert_real_to_algebraic(const Lazy_exact_nt<ET>& r) const {
    typename CGAL::Coercion_traits<ET, typename Algebraic::ET>::Cast cast;
    return cast(r.exact());
  }
  template <class NT>
  Algebraic convert_real_to_algebraic(const NT& r) const {
    typename CGAL::Coercion_traits<NT, Algebraic>::Cast cast;
    return cast(r);
  }

  template <typename Point_2,typename Number_type>
  void convert_point_y_coordinate (const Point_2& p, const Sqrt_ext& p_x, Number_type& y) const
  {
    y = convert_real_to_algebraic(p.y());
  }

#else

  Algebraic convert_real_to_algebraic(const Algebraic_real_1& r)  const
  {
    typename Traits_algebraic_kernel::Compute_polynomial_1 compute_poly;
    Polynomial_1 poly = compute_poly(r);

    typename Traits_algebraic_kernel::Isolate_1 isolate;
    Bound_pair bound_pair = isolate(r,poly);
    CORE_big_float_interval bound_interval(bound_pair);

    CORE_big_float f = bound_interval.first;
    CORE_big_float s = bound_interval.second;

    if ((f.isExact() == false) || (s.isExact() == false))
    {
      bound_interval.first.makeExact();
      bound_interval.second.makeExact();

      bound_interval.first  -= 0.000000000001;
      bound_interval.second += 0.000000000001;
    }
    Polynomial core_poly = convert_polynomial (poly);
    return Algebraic (core_poly, bound_interval);
  }


  template <typename Point_2,typename Number_type>
  void convert_point_y_coordinate (const Point_2& p, const Algebraic_real_1& p_x, Number_type& y) const
  {
    Algebraic x = convert_real_to_algebraic(p.x());
    Polynomial core_numer = convert_polynomial(p.rational_function().numer());
    Polynomial core_denom = convert_polynomial(p.rational_function().denom());

    Algebraic y_numer = core_numer.eval(x);
    Algebraic y_denom = core_denom.eval(x);

    y = (y_numer / y_denom);
  }
#endif

  template <typename Point_2,typename Number_type>
  void convert_point (const Point_2& p,Number_type& x, Number_type& y) const
  {
    x = convert_real_to_algebraic(p.x());
    convert_point_y_coordinate(p, p.x(), y);
  }


};

template <typename Traits_3_>
class Lines_through_segments_get_algebraic_number_adapt {
private:
  using Traits_3 = Traits_3_;

  using Algebraic = typename Traits_3::Alg_kernel::FT;
  using Integer = CORE::BigInt;
  using Rational_kernel = typename Traits_3::Rational_kernel;
  using Rational = typename Rational_kernel::FT;
#if USE_SQRT_TRAITS
  using AK1 = CGAL::Algebraic_kernel_2_1<Rational>;
#else
  using AK1 = CGAL::Algebraic_kernel_d_1<Integer>;
#endif
  using Rational_arc_traits_arr_on_plane_2 =
    CGAL::Arr_rational_function_traits_2<AK1>;

  using Algebraic_real_1 =
    typename Rational_arc_traits_arr_on_plane_2::Algebraic_real_1;

private:
  template <typename NT>
  Algebraic get_algebraic_number(const NT &x) { return Algebraic(x); }

  Algebraic get_algebraic_number(const Algebraic_real_1 &x) {
    Lines_through_segments_traits_on_plane_adapt<Traits_3>
      traits_on_plane_adapt;
    return traits_on_plane_adapt.convert_real_to_algebraic(x);
  }

public:
  template <typename NT>
  Algebraic operator()(const NT &x) { return get_algebraic_number(x); }
};

} //namespace CGAL

#endif /*LINES_THROUGH_SEGMENTS_TRAITS_2_ADAPT_H*/
