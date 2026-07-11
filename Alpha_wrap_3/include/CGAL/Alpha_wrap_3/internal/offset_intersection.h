// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_OFFSET_INTERSECTION_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_OFFSET_INTERSECTION_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/number_utils.h>

#include <boost/algorithm/clamp.hpp>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

template <typename AABBTree,
          typename AABBTraversalTraits>
struct AABB_tree_oracle_helper;

template <typename AABBTree, typename AABBTraversalTraits>
struct AABB_distance_oracle
{
  using FT = typename AABBTree::FT;
  using Point_3 = typename AABBTree::Point;

  using AABB_helper = AABB_tree_oracle_helper<AABBTree, AABBTraversalTraits>;

  AABB_distance_oracle(const AABBTree& tree) : tree(tree) { }

  FT operator()(const Point_3& p) const
  {
    return approximate_sqrt(AABB_helper::squared_distance(p, tree));
  }

public:
  const AABBTree& tree;
};

// @todo even with EPECK, the precision cannot be 0 (otherwise it will not converge),
// thus exactness is pointless. Might as well use a cheap kernel (e.g. SC<double>), as long
// as there exists a mechanism to catch when the cheap kernel fails to converge (iterations?
// see also Tr_3::locate() or Mesh_3::Robust_intersection_traits_3.h)
template <class Kernel, class DistanceOracle>
class Offset_intersection
{
  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;

public:
  Offset_intersection(const DistanceOracle& oracle,
                      const FT& off,
                      const FT& prec,
                      const FT& lip)
    : dist_oracle(oracle), offset(off), precision(prec), lipschitz(lip)
  { }

  bool first_intersection(const Point_3& s,
                          const Point_3& t,
                          Point_3& output_pt)
  {
    source = s;
    target = t;

#define CGAL_AW3_OUTPUT_SPHERE_MARCHING_STEPS
#ifdef CGAL_AW3_OUTPUT_SPHERE_MARCHING_STEPS
    if (ray == 0) {
      std::remove(output_file.c_str());
      output_header();
    }
    output_step(s, "source", ray, "init");
    output_step(t, "target", ray, "init");
#endif

    bool output;

    std::array<std::string, 4> vars = {
      "DICHOTOMY", "SPHERE_MARCHING_OLD", "SPHERE_MARCHING", "RELAXED"
    };
    std::array<int,4> vals;
    transform(vars.begin(), vars.end(), vals.begin(), getenv);
    int max = *std::max_element(vals.begin(), vals.end());

    if (max == 0) {
      // default
      output = trace("SPHERE_MARCHING", s, t, output_pt);
    } else {
      // Run the algos with positive environnment variables.
      // Run the algo with highest value last, so that its output gets used.
      for (int i = 1; i <= max; i++) {
        for (const std::string& var : vars) {
          if (i == getenv(var)) {
            std::cout << "Tracing:" << var << std::endl;
            output = trace(var, s, t, output_pt);
          }
        }
      }
    }

#ifdef CGAL_AW3_OUTPUT_SPHERE_MARCHING_STEPS
    output_step(output_pt, "intersection", ray, "init");
    ray++;
#endif

    return output;
  }


private:
  Point_3 source;
  Point_3 target;
  FT seg_length;
  Vector_3 seg_unit_v;
  DistanceOracle dist_oracle;
  FT offset;
  FT precision;
  FT lipschitz;
#ifdef CGAL_AW3_OUTPUT_SPHERE_MARCHING_STEPS
  inline static int ray = 0;
  std::string output_file = "steps.csv";
#endif

  static int getenv(std::string var) {
    char* val = std::getenv(var.c_str());
    if (val) {
      return std::atoi(val);
    } else {
      return 0;
    }
  }

  bool trace(std::string algo,
             const Point_3& s,
             const Point_3& t,
             Point_3& output_pt) {
    if (algo == "DICHOTOMY") {
      return dichotomy(s, t, output_pt);
    } else if (algo == "SPHERE_MARCHING_OLD") {
      return sphere_tracing_old(s, t, output_pt);
    } else if (algo == "SPHERE_MARCHING") {
      return sphere_tracing(s, t, output_pt);
    } else if (algo == "RELAXED") {
      return relaxed_sphere_tracing(s, t, output_pt);
    }
  }

#ifdef CGAL_AW3_OUTPUT_SPHERE_MARCHING_STEPS
  void output_header() {
    std::ofstream stream(output_file, std::ios_base::app);
    stream.precision(17);
    stream << "ray" << ","
             << "algo" << ","
             << "type" << ","
             << "x" << ","
             << "y" << ","
             << "z" << "\n";
    stream.close();
  }

  void output_step(Point_3 point, std::string type, int ray, std::string algo) {
    std::ofstream stream(output_file, std::ios_base::app);
    stream.precision(17);
    stream << ray << ","
             << algo << ","
             << type << ","
             << point.x() << ","
             << point.y() << ","
             << point.z() << "\n";
    stream.close();
  }
#endif

  bool dichotomy (const Point_3& s, const Point_3& t, Point_3& output_pt) {
    seg_length = approximate_sqrt(squared_distance(s, t));
    seg_unit_v = (t - s) / seg_length;
    const Point_2 p0 { 0, dist_oracle(source) };
    const Point_2 p1 { seg_length, dist_oracle(target) };
    return recursive_dichotomic_search(p0, p1, output_pt);
  }

  template <class Point>
  bool recursive_dichotomic_search(const Point_2& s, const Point_2& t,
                                   Point& output_pt)
  {

#ifdef CGAL_AW3_OUTPUT_SPHERE_MARCHING_STEPS
    const Point_3 segment_data { s.x(), s.y(), t.y() };
    output_step(segment_data, "segment", ray, "dichotomy");
#endif

    if(CGAL::abs(s.x() - t.x()) < precision)
    {
      if(CGAL::abs(s.y() - offset) < precision)
      {
        const FT x_clamp = boost::algorithm::clamp<FT>(s.x(), FT{0}, seg_length);
        output_pt = source + (seg_unit_v * x_clamp);
        return true;
      }

      return false;
    }

    const bool sign_s = (s.y() > offset);
    const bool sign_t = (t.y() > offset);
    const FT gs_a = (sign_s) ? -lipschitz : lipschitz;
    const FT gs_b = s.y() - (gs_a * s.x());
    const FT gt_a = (sign_t) ? lipschitz : -lipschitz;

    const FT gt_b = t.y() - (gt_a * t.x());
    FT ms = (offset - gs_b) / gs_a;
    FT mt = (offset - gt_b) / gt_a;

    // early exit if there is no intersection
    if(sign_s == sign_t)
    {
      FT ui = (gt_b - gs_b) / (gs_a - gt_a);
      const FT gs_ui = (gs_a * ui) + gs_b;
      if((sign_s && (gs_ui > offset)) || (!sign_s && (gs_ui < offset)))
      {
        if(CGAL::abs(s.y() - offset) < precision)
        {
          const FT x_clamp = boost::algorithm::clamp<FT>(s.x(), FT{0}, seg_length);
          output_pt = source + (seg_unit_v * x_clamp);
          return true;
        }
        else if(CGAL::abs(t.y() - offset) < precision)
        {
          const FT x_clamp = boost::algorithm::clamp<FT>(t.x(), FT{0}, seg_length);
          output_pt = source + (seg_unit_v * x_clamp);
          return true;
        }

        return false;
      }
      else
      {
        ms = boost::algorithm::clamp<FT>(ms, FT{0}, seg_length);
        ui = boost::algorithm::clamp<FT>(ui, FT{0}, seg_length);
        mt = boost::algorithm::clamp<FT>(mt, FT{0}, seg_length);
        const Point_2 ms_pt { ms, dist_oracle(source + (seg_unit_v * ms)) };
        const Point_2 ui_pt { ui, dist_oracle(source + (seg_unit_v * ui)) };
        const Point_2 mt_pt { mt, dist_oracle(source + (seg_unit_v * mt)) };

        if(CGAL::abs(ms_pt.y() - offset) < precision)
        {
          const FT x_clamp = boost::algorithm::clamp<FT>(ms_pt.x(), FT{0}, seg_length);
          output_pt = source + (seg_unit_v * x_clamp);
          return true;
        }
        else if(CGAL::abs(ui_pt.y() - offset) < precision)
        {
          const FT x_clamp = boost::algorithm::clamp<FT>(ui_pt.x(), FT{0}, seg_length);
          output_pt = source + (seg_unit_v * x_clamp);
          return true;
        }
        else if(CGAL::abs(mt_pt.y() - offset) < precision)
        {
          const FT x_clamp = boost::algorithm::clamp<FT>(mt_pt.x(), FT{0}, seg_length);
          output_pt = source + (seg_unit_v * x_clamp);
          return true;
        }

        return (recursive_dichotomic_search(ms_pt, ui_pt, output_pt) ||
                recursive_dichotomic_search(ui_pt, mt_pt, output_pt));
      }
    }
    else // there is an intersection
    {
      if(CGAL::abs(mt - ms) <= precision) // linear approximation
      {
        const FT fsft_a = (t.y() - s.y()) / (t.x() - s.x());
        const FT fsft_b = s.y() - fsft_a * s.x();
        FT m_fsft;
        if(fsft_a == FT{0})
        {
          if(CGAL::abs(s.y() - offset) < precision)
            m_fsft = s.x();
          else
            return false;
        }
        else
        {
          m_fsft = (offset - fsft_b) / fsft_a;
        }
        m_fsft = boost::algorithm::clamp<FT>(m_fsft, FT{0}, seg_length);
        output_pt = source + (seg_unit_v * m_fsft);
        return true;
      }
      else
      {
        FT m = (ms + mt) / FT{2};
        ms = boost::algorithm::clamp<FT>(ms, FT{0}, seg_length);
        m = boost::algorithm::clamp<FT>(m, FT{0}, seg_length);
        mt = boost::algorithm::clamp<FT>(mt, FT{0}, seg_length);

        const Point_2 ms_pt { ms, dist_oracle(source + (seg_unit_v * ms)) };
        const Point_2 m_pt { m, dist_oracle(source + (seg_unit_v * m)) };
        const Point_2 mt_pt { mt, dist_oracle(source + (seg_unit_v * mt)) };

        return (recursive_dichotomic_search(ms_pt, m_pt, output_pt) ||
                recursive_dichotomic_search(m_pt, mt_pt, output_pt));
      }
    }
  }

  bool sphere_tracing_old(const Point_3& s,
                          const Point_3& t,
                          Point_3& output_pt)
  {
#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
    static int calls = -1;
    ++calls;
    std::cout << "Sphere march between " << s << " and " << t << " (#" << calls << ")" << std::endl;
#endif

    CGAL_precondition(s != t);

    Point_3 current_pt = s;
    Point_3 closest_point = dist_oracle.tree.closest_point(current_pt, t);
    FT current_dist = approximate_sqrt(squared_distance(current_pt, closest_point)) - offset;

    const FT sq_seg_length = squared_distance(s, t);
    const FT seg_length = approximate_sqrt(sq_seg_length);
    if(is_zero(seg_length))
    {
      closest_point = dist_oracle.tree.closest_point(t);
      current_dist = approximate_sqrt(squared_distance(t, closest_point)) - offset;
      output_pt = t;
      return (CGAL::abs(current_dist) < precision);
    }

    const Vector_3 seg_unit_v = (t - s) / seg_length;

    for(;;)
    {
#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
      std::cout << "current point " << current_pt << std::endl;
      std::cout << "current dist " << current_dist << std::endl;
#endif

      if(CGAL::abs(current_dist) < precision)
      {
#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
        std::cout << "Intersection point " << current_pt << std::endl;
#endif
        output_pt = current_pt;
        return true;
      }

      current_pt = current_pt + (current_dist * seg_unit_v);

#ifdef CGAL_AW3_OUTPUT_SPHERE_MARCHING_STEPS
      output_step(current_pt, "step", ray, "sphere-marching-old");
#endif

      if(squared_distance(s, current_pt) > sq_seg_length)
      {
#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
    std::cout << "No intersection " << std::endl;
#endif
         return false;
      }
      // previous closest point is a good hint for the next closest point
      closest_point = dist_oracle.tree.closest_point(current_pt, closest_point);
      current_dist = approximate_sqrt(squared_distance(current_pt, closest_point)) - offset;
    }

    return false;
  }

  bool sphere_tracing(const Point_3& s,
                      const Point_3& t,
                      Point_3& output_pt)
  {
#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
    static int calls = -1;
    ++calls;
    std::cout << "Sphere march between " << s << " and " << t << " (#" << calls << ")" << std::endl;
#endif

    CGAL_precondition(s != t);

    Point_3 current_pt = s;
    FT current_dist;

    const FT sq_seg_length = squared_distance(s, t);
    const FT seg_length = approximate_sqrt(sq_seg_length);
    if(is_zero(seg_length))
    {
      current_dist = dist_oracle(current_pt) - offset;
      output_pt = t;
      return (CGAL::abs(current_dist) < precision);
    }

    const Vector_3 direction = (t - s) / seg_length;

    while(squared_distance(s, current_pt) <= sq_seg_length)
    {
      current_dist = dist_oracle(current_pt) - offset;

#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
      std::cout << "current point " << current_pt << std::endl;
      std::cout << "current dist " << current_dist << std::endl;
#endif

      if(CGAL::abs(current_dist) < precision)
      {
#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
      std::cout << "Intersection point " << current_pt << std::endl;
#endif
        output_pt = current_pt;
        return true;
      }

      current_pt = current_pt + (current_dist * direction);

#ifdef CGAL_AW3_OUTPUT_SPHERE_MARCHING_STEPS
      output_step(current_pt, "step", ray, "sphere-marching");
#endif
    }

#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
    std::cout << "No intersection " << std::endl;
#endif
    return false;
  }

  bool relaxed_sphere_tracing(const Point_3& s,
                              const Point_3& t,
                              Point_3& output_pt)
  {
#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
    static int calls = -1;
    ++calls;
    std::cout << "Sphere march between " << s << " and " << t << " (#" << calls << ")" << std::endl;
#endif

    CGAL_precondition(s != t);

    Point_3 current_pt = s;
    FT current_dist;
    FT omega = 1.2;

    const FT sq_seg_length = squared_distance(s, t);
    const FT seg_length = approximate_sqrt(sq_seg_length);
    if(is_zero(seg_length))
    {
      current_dist = dist_oracle(current_pt) - offset;
      output_pt = t;
      return (CGAL::abs(current_dist) < precision);
    }

    const Vector_3 direction = (t - s) / seg_length;

    FT current_step;
    FT previous_dist = 0;
    bool sorFail = false;

    current_dist = dist_oracle(current_pt) - offset;
    current_step = current_dist;

    int steps = 0;

    while(squared_distance(s, current_pt) <= sq_seg_length * omega * omega || sorFail)
    {
#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
      std::cout << "current point " << current_pt << std::endl;
      std::cout << "current dist " << current_dist << std::endl;
#endif

      sorFail = omega > 1 &&
        (CGAL::abs(current_dist) + CGAL::abs(previous_dist)) < current_step;
      if (sorFail) {
#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
        std::cout << "Overstepped abs(" << current_dist << ") + abs(" << previous_dist << ") < " << current_step << std::endl;
#endif

        current_step -= omega * current_step;
        omega = 1;

#ifdef CGAL_AW3_OUTPUT_SPHERE_MARCHING_STEPS
        output_step(current_pt, "overstep", ray, "relaxed");
#endif
      } else {
        current_step = current_dist * omega;

#ifdef CGAL_AW3_OUTPUT_SPHERE_MARCHING_STEPS
        if (steps > 0) {
          if (omega > 1) {
            output_step(current_pt, "relaxed", ray, "relaxed");
          } else {
            output_step(current_pt, "step", ray, "relaxed");
          }
        }
#endif
      }

      if(CGAL::abs(current_dist) < precision && !sorFail)
      {
#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
        std::cout << "Intersection point " << current_pt << std::endl;
#endif
        output_pt = current_pt;
        return true;
      }

      current_pt = current_pt + (current_step * direction);

      previous_dist = current_dist;
      current_dist = dist_oracle(current_pt) - offset;
      steps++;
    }

#ifdef CGAL_AW3_DEBUG_SPHERE_MARCHING
    std::cout << "No intersection. Current point " << current_pt << " with " << squared_distance(s, current_pt) << " <= " << sq_seg_length << " * " << omega << std::endl;

#endif
    return false;
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_OFFSET_INTERSECTION_H
