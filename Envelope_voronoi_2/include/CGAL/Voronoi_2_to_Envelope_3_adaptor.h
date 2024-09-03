// Copyright (c) 2008  Tel-Aviv University (Israel).
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
// Author(s)     : Ophir Setter       <ophirset@post.tau.ac.il>

#ifndef CGAL_VORONOI_2_TO_ENVELOPE_3_ADAPTOR_H
#define CGAL_VORONOI_2_TO_ENVELOPE_3_ADAPTOR_H

namespace CGAL {

template <typename EnvelopeVoronoiTraits_>
class Voronoi_2_to_Envelope_3_adaptor : public EnvelopeVoronoiTraits_ {
  public:
  using Envelope_voronoi_traits = EnvelopeVoronoiTraits_;
  using Base = Envelope_voronoi_traits;

  using Point_2 = typename Envelope_voronoi_traits::Point_2;
  using X_monotone_curve_2 =
    typename Envelope_voronoi_traits::X_monotone_curve_2;
  using Multiplicity = typename    Envelope_voronoi_traits::Multiplicity;

  using Xy_monotone_surface_3 = typename Envelope_voronoi_traits::Site_2;
  using Surface_3 = Xy_monotone_surface_3;

public:
  Voronoi_2_to_Envelope_3_adaptor<Envelope_voronoi_traits>
    (const Envelope_voronoi_traits & traits) : Base (traits)
  {}

  Voronoi_2_to_Envelope_3_adaptor<Envelope_voronoi_traits> () {}

  class Make_xy_monotone_3 {
  public:
    template <typename OutputIterator>
    OutputIterator operator()(const Surface_3& s, bool /* is_lower */,
                              OutputIterator o) const {
      *o++ = s;
      return o;
    }
  };

  Make_xy_monotone_3 make_xy_monotone_3_object() const
  { return Make_xy_monotone_3(); }

  class Compare_z_at_xy_3 {
  protected:
    const Envelope_voronoi_traits * _p_traits;

  public:
    Compare_z_at_xy_3(const Envelope_voronoi_traits * p_traits) :
      _p_traits(p_traits)
    {}

    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    { return _p_traits->compare_distance_at_point_2_object() (h1, h2, p); }

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    {
      Arr_parameter_space sbx =
        _p_traits->parameter_space_in_x_2_object() (cv, ARR_MIN_END);
      Arr_parameter_space sby =
        _p_traits->parameter_space_in_y_2_object() (cv, ARR_MIN_END);

      if (sbx == ARR_INTERIOR && sby == ARR_INTERIOR) {
        Point_2 p = _p_traits->construct_min_vertex_2_object() (cv);
        Comparison_result res = _p_traits->compare_distance_at_point_2_object()
          (h1, h2, p);
        if (res != EQUAL) return res;
      }

      Arr_parameter_space tbx =
        _p_traits->parameter_space_in_x_2_object() (cv, ARR_MAX_END);
      Arr_parameter_space tby =
        _p_traits->parameter_space_in_y_2_object() (cv, ARR_MAX_END);
      if (tbx == ARR_INTERIOR && tby == ARR_INTERIOR) {
        Point_2 p = _p_traits->construct_max_vertex_2_object() (cv);
        Comparison_result res = _p_traits->compare_distance_at_point_2_object()
          (h1, h2, p);
        if (res != EQUAL)
          return res;
      }

      Point_2 p = _p_traits->construct_point_on_x_monotone_2_object() (cv);
      return _p_traits->compare_distance_at_point_2_object() (h1, h2, p);
    }

    Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    { return _p_traits->compare_dominance_2_object() (h1, h2); }
  };

  Compare_z_at_xy_3 compare_z_at_xy_3_object() const
  { return Compare_z_at_xy_3(this); }

  class Compare_z_at_xy_above_3 {
  protected:
    const Envelope_voronoi_traits * _p_traits;

  public:
    Compare_z_at_xy_above_3 (const Envelope_voronoi_traits * p_traits) :
      _p_traits (p_traits)
    {}

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    { return _p_traits->compare_distance_above_2_object() (h1, h2, cv); }
  };

  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object() const
  { return Compare_z_at_xy_above_3(this); }

  class Compare_z_at_xy_below_3 {
  protected:
    Compare_z_at_xy_above_3 _cmp_above;

  public:
    Compare_z_at_xy_below_3(const Compare_z_at_xy_above_3& cmp_above) :
      _cmp_above (cmp_above)
    {}

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const
    { return CGAL::opposite(_cmp_above(cv, h1, h2)); }

  };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object() const
  { return Compare_z_at_xy_below_3(compare_z_at_xy_above_3_object()); }

  class Construct_projected_boundary_2 {
  public:

    template <class OutputIterator>
      OutputIterator operator()(const Xy_monotone_surface_3& s,
                                OutputIterator o) const
    { return o; }
  };

  Construct_projected_boundary_2 construct_projected_boundary_2_object() const
  { return Construct_projected_boundary_2(); }

  class Construct_projected_intersections_2 {
  protected:
    const Envelope_voronoi_traits * _p_traits;

  public:
    Construct_projected_intersections_2(const Envelope_voronoi_traits * p_traits)
      : _p_traits (p_traits)
    {}

    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const {
      using X_curve_list = std::list<X_monotone_curve_2>;
      using Intersection_curve = std::pair<X_monotone_curve_2, Multiplicity>;

      X_curve_list x_curves;
      _p_traits->construct_bisector_2_object()(s1, s2,
                                               std::back_inserter(x_curves));

      for (auto it = x_curves.begin(); it != x_curves.end(); ++it)
        *o++ = Intersection_curve(*it, 1);
      return o;
    }
  };

  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const
  { return Construct_projected_intersections_2(this); }

};

} //namespace CGAL


#endif
