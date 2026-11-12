
#ifndef CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_GEOMETRY_UTILS_H
#define CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_GEOMETRY_UTILS_H

#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>

#include <cmath>
#include <vector>

namespace CGAL {
namespace internal {
namespace Splat_surface_reconstruction_3 {

template <typename Vector_3, typename FT>
Vector_3 normalize_vec(const Vector_3& n) {
  const double len2 = CGAL::to_double(n.squared_length());
  if (len2 <= 0.0) {
    return CGAL::NULL_VECTOR;
  }
  const Vector_3 normalized = n / std::sqrt(len2);
  CGAL_assertion(std::abs(normalized.squared_length() - FT(1.0)) < 1e-8);
  return normalized;
}



/**
 * @brief Tests whether two 2D segments intersect in their interiors.
 *
 * Endpoint touching is not considered a strict intersection. When
 * `allow_overlap` is false, overlapping segments are considered intersecting.
 *
 * @param a First endpoint of the first segment.
 * @param b Second endpoint of the first segment.
 * @param c First endpoint of the second segment.
 * @param d Second endpoint of the second segment.
 * @param allow_overlap Whether overlapping segments should be ignored.
 *
 * @return `true` if the segments violate the local intersection constraint.
*/
template <typename Point_2, typename Kernel>
bool intersect_2d(const Point_2& a,
                  const Point_2& b,
                  const Point_2& c,
                  const Point_2& d,
                  bool allow_overlap = false) {
  const CGAL::Segment_2<Kernel> s1(a, b);
  const CGAL::Segment_2<Kernel> s2(c, d);

  auto inter = CGAL::intersection(s1, s2);
  if (!inter) {
    return false;
  }

  Point_2 ip;
  if (CGAL::assign(ip, *inter)) {
    // Allow touching at endpoints.
    return !(ip == a || ip == b || ip == c || ip == d);
  }

  // If the intersection is not a point, it is an overlapping segment.
  if (!allow_overlap) {
    return true;
  }
  return false;
};

template <typename Point_3, typename Vector_3, typename FT>
bool point_in_cylinder(const Point_3& p, const Point_3& c, const Vector_3& n, FT radius, FT height) {
  Vector_3 pc = p - c;
  FT dist_along_normal = pc * n;
  if (dist_along_normal < -height || dist_along_normal > height) {
    return false;
  }
  Vector_3 proj = pc - dist_along_normal * n;
  return std::sqrt(proj.squared_length()) <= radius;
}

/**
 * @brief Finds the circle defined by the intersection of two spheres centered at the parent points with the given radius.
 *
 * @param parent_a First point.
 * @param parent_b Second point.
 * @param radius Sphere radius.
 * @param circle_center Output circle center.
 * @param circle_normal Output circle normal.
 * @param circle_radius Output circle radius.
 *
 * @return True if a valid circle is found, false otherwise.
 */
template <typename Point_3, typename Vector_3, typename FT>
bool find_circle_between_points(
  const Point_3& parent_a,
  const Point_3& parent_b,
  const FT radius,
  Point_3& circle_center,
  Vector_3& circle_normal,
  FT& circle_radius) {

  const Vector_3 ab = parent_b - parent_a;
  const FT d2 = ab.squared_length();
  if (d2 <= FT(0)) {
    return false;
  }

  // Circle radius from the two-parent construction.
  const FT circle_r2 = radius * radius - d2 / FT(4);
  if (circle_r2 <= FT(0)) {
    return false;
  }

  circle_radius = CGAL::sqrt(circle_r2);

  // Circle center = midpoint of the two parents.
  circle_center = Point_3(
    (parent_a.x() + parent_b.x()) / FT(2),
    (parent_a.y() + parent_b.y()) / FT(2),
    (parent_a.z() + parent_b.z()) / FT(2)
  );

  // Normal of the circle plane: along the parent-parent direction.
  const double ab_len2 = CGAL::to_double(d2);
  if (ab_len2 <= 0.0) {
    return false;
  }
  circle_normal = ab * (1.0 / std::sqrt(ab_len2));

  return true;
}

/**
 * @brief Computes an orthonormal tangent frame for a normal.
 *
 * The returned vectors span the plane orthogonal to the input normal.
 *
 * @param n Surface normal.
 *
 * @return Two approximately unit-length tangent vectors.
*/
template <typename Point_3, typename Vector_3, typename FT>
std::vector<Vector_3> compute_local_tangent_frame(const Vector_3& n) {
  // Find a vector that is not parallel to n to construct the tangent frame.
  Vector_3 temp;
  if (std::abs(CGAL::to_double(n.x())) < 0.9)
  temp = Vector_3(1,0,0);
  else
  temp = Vector_3(0,1,0);

  Vector_3 u = CGAL::cross_product(n, temp);
  Vector_3 v = CGAL::cross_product(n, u);

  // Each vector should be orthogonal to the others.
  CGAL_assertion(u*v < 1e-8);
  CGAL_assertion(n*u < 1e-8);
  CGAL_assertion(n*v < 1e-8);

  // None should be the zero vector.
  CGAL_assertion(u.squared_length() > 0);
  CGAL_assertion(v.squared_length() > 0);
  CGAL_assertion(n.squared_length() > 0);

  // Normalize u and v to have unit length.
  const double len_u = std::sqrt(CGAL::to_double(u.squared_length()));
  const double len_v = std::sqrt(CGAL::to_double(v.squared_length()));
  if (len_u > 0) u = u / len_u;
  if (len_v > 0) v = v / len_v;

  return {u, v};
}

template <typename Point_3, typename Vector_3, typename FT>
std::vector<Point_3> circle_intersects_splat(
  const Point_3& circle_center,
  const Vector_3& circle_normal,
  const FT circle_radius,
  const Point_3& splat_center,
  const Vector_3& splat_normal,
  const FT splat_radius)
{
  std::vector<Point_3> result;

  // Build an orthonormal basis (u, v) for the circle plane.
  std::vector<Vector_3> frame = compute_local_tangent_frame<Point_3, Vector_3, FT>(circle_normal);
  const Vector_3& u = frame[0];
  const Vector_3& v = frame[1];

  // Normalize the splat normal.
  Vector_3 ns = splat_normal;
  const double ns_len2 = CGAL::to_double(ns.squared_length());
  if (ns_len2 <= 0.0)
    return result;

  ns = ns * (1.0 / std::sqrt(ns_len2));

  // Circle points:
  //   x(t) = circle_center + circle_radius * (cos(t) * u + sin(t) * v)
  //
  // Splat plane equation:
  //   dot(ns, x - splat_center) = 0
  //
  // This becomes:
  //   A cos(t) + B sin(t) = C
  const FT A = circle_radius * (u * ns);
  const FT B = circle_radius * (v * ns);
  const FT C = (splat_center - circle_center) * ns;

  const FT R2 = A * A + B * B;
  if (R2 <= FT(0))
    return result;

  const FT R = CGAL::sqrt(R2);

  if (std::abs(CGAL::to_double(C)) > CGAL::to_double(R))
    return result;

  const FT phi = std::atan2(CGAL::to_double(B), CGAL::to_double(A));
  const FT delta = std::acos(CGAL::to_double(C / R));

  const FT t1 = phi + delta;
  const FT t2 = phi - delta;

  const auto make_point = [&](FT t) -> Point_3
  {
    const FT ct = std::cos(CGAL::to_double(t));
    const FT st = std::sin(CGAL::to_double(t));

    return Point_3(
      circle_center.x() + circle_radius * (ct * u.x() + st * v.x()),
      circle_center.y() + circle_radius * (ct * u.y() + st * v.y()),
      circle_center.z() + circle_radius * (ct * u.z() + st * v.z()));
  };

  const FT splat_r2 = splat_radius * splat_radius;

  const Point_3 p1 = make_point(t1);
  const Point_3 p2 = make_point(t2);

  // Keep only points that lie inside the splat disk.
  if (CGAL::squared_distance(p1, splat_center) <= splat_r2)
    result.push_back(p1);

  // Avoid duplicate if the two solutions coincide numerically.
  if (CGAL::squared_distance(p2, splat_center) <= splat_r2 &&
      CGAL::squared_distance(p2, p1) > FT(1e-12))
    result.push_back(p2);

  return result;
}

template <typename Point_3, typename Vector_3, typename FT>
bool triangle_intersects_aabb(
  const Point_3& t1,
  const Point_3& t2,
  const Point_3& t3,
  const Point_3& aabb_min,
  const Point_3& aabb_max)
{
  using Kernel = typename CGAL::Kernel_traits<Point_3>::Kernel;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Iso_cuboid_3 = typename Kernel::Iso_cuboid_3;

  const Triangle_3 triangle(t1, t2, t3);
  const Iso_cuboid_3 box(aabb_min, aabb_max);

  return CGAL::do_intersect(triangle, box);
}

template <typename Point_3, typename Vector_3, typename FT>
bool line_intersects_aabb(
  const Point_3& l1,
  const Point_3& l2,
  const Point_3& aabb_min,
  const Point_3& aabb_max)
{
  using Kernel = typename CGAL::Kernel_traits<Point_3>::Kernel;
  using Line_3 = typename Kernel::Line_3;
  using Iso_cuboid_3 = typename Kernel::Iso_cuboid_3;

  return CGAL::do_intersect(
    Line_3(l1, l2),
    Iso_cuboid_3(aabb_min, aabb_max));
}

template <typename Point_3, typename Vector_3, typename FT>
std::vector<Point_3> circle_intersects_triangle(
  const Point_3& circle_center,
  const Vector_3& circle_normal,
  const FT circle_radius,
  const Point_3& t1,
  const Point_3& t2,
  const Point_3& t3)
{
  using Kernel = typename CGAL::Kernel_traits<Point_3>::Kernel;

  using Segment_3 = typename Kernel::Segment_3;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Plane_3 = typename Kernel::Plane_3;

  std::vector<Point_3> result;

  const Plane_3 circle_plane(circle_center, circle_normal);
  const Triangle_3 triangle(t1, t2, t3);

  const auto intersection = CGAL::intersection(circle_plane, triangle);

  if (!intersection)
    return result;

  const FT radius_squared = circle_radius * circle_radius;

  auto add_point = [&](const Point_3& p)
  {
    if (std::find(result.begin(), result.end(), p) == result.end())
      result.push_back(p);
  };

  // Check whether a point lies on the circle circumference.
  auto check_point = [&](const Point_3& p)
  {
    if ((p - circle_center).squared_length() == radius_squared)
      add_point(p);
  };

  // Intersect a segment with the circle circumference.
  auto check_segment = [&](const Segment_3& segment)
  {
    const Point_3& a = segment.source();
    const Vector_3 direction = segment.target() - a;
    const Vector_3 offset = a - circle_center;

    const FT A = direction.squared_length();
    const FT B = FT(2) * (offset * direction);
    const FT C = offset.squared_length() - radius_squared;

    const FT discriminant = B * B - FT(4) * A * C;

    if (discriminant < FT(0))
      return;

    const FT sqrt_discriminant = CGAL::sqrt(discriminant);
    const FT denominator = FT(2) * A;

    const FT t1 = (-B - sqrt_discriminant) / denominator;
    const FT t2 = (-B + sqrt_discriminant) / denominator;

    if (t1 >= FT(0) && t1 <= FT(1))
      add_point(a + direction * t1);

    if (t2 >= FT(0) && t2 <= FT(1))
      add_point(a + direction * t2);
  };

  std::visit([&](const auto& geometry)
  {
    using T = std::decay_t<decltype(geometry)>;

    if constexpr (std::is_same_v<T, Point_3>)
    {
      check_point(geometry);
    }
    else if constexpr (std::is_same_v<T, Segment_3>)
    {
      check_segment(geometry);
    }
    else if constexpr (std::is_same_v<T, Triangle_3>)
    {
      check_segment(Segment_3(geometry.vertex(0), geometry.vertex(1)));
      check_segment(Segment_3(geometry.vertex(1), geometry.vertex(2)));
      check_segment(Segment_3(geometry.vertex(2), geometry.vertex(0)));
    }
  }, *intersection);

  return result;
}


template <typename Point_3, typename Vector_3, typename FT>
std::vector<Point_3> line_intersects_sphere(
  const Point_3& line_point_1,
  const Point_3& line_point_2,
  const Point_3& sphere_center,
  FT sphere_radius)
{
  std::vector<Point_3> intersections;

  const Vector_3 direction = line_point_2 - line_point_1;
  const Vector_3 offset = line_point_1 - sphere_center;

  const FT a = direction.squared_length();

  // Degenerate segment.
  if (a == FT(0))
    return intersections;

  const FT b = FT(2) * (offset * direction);
  const FT c = offset.squared_length() - sphere_radius * sphere_radius;

  const FT discriminant = b * b - FT(4) * a * c;

  if (discriminant < FT(0))
    return intersections;

  const FT sqrt_discriminant = CGAL::sqrt(discriminant);
  const FT denominator = FT(2) * a;

  const FT t1 = (-b - sqrt_discriminant) / denominator;
  const FT t2 = (-b + sqrt_discriminant) / denominator;

  // Keep intersections lying on the segment.
  if (t1 >= FT(0) && t1 <= FT(1))
    intersections.push_back(line_point_1 + t1 * direction);

  // Avoid adding the same point twice for tangential intersections.
  if (t2 >= FT(0) && t2 <= FT(1) && t2 != t1)
    intersections.push_back(line_point_1 + t2 * direction);

  return intersections;
}


template <typename Point_3, typename Vector_3, typename FT, typename Kernel>
bool triple_plane_intersection(
  const Point_3& p1, const Vector_3& n1, const FT& r1,
  const Point_3& p2, const Vector_3& n2, const FT& r2,
  const Point_3& p3, const Vector_3& n3, const FT& r3,
  FT minimum_angle,
  Point_3& intersection_point)
{
  // Reject invalid normals.
  if (n1 == CGAL::NULL_VECTOR ||
      n2 == CGAL::NULL_VECTOR ||
      n3 == CGAL::NULL_VECTOR)
  {
    return false;
  }

  // Reject intersections where the angle between any pair of normals is too small.
  if (std::abs(CGAL::approximate_angle(n1, n2)) < minimum_angle ||
      std::abs(CGAL::approximate_angle(n1, n3)) < minimum_angle ||
      std::abs(CGAL::approximate_angle(n2, n3)) < minimum_angle)
  {
    return false;
  }

  // Construct the three supporting planes.
  const typename Kernel::Plane_3 plane1(p1, n1);
  const typename Kernel::Plane_3 plane2(p2, n2);
  const typename Kernel::Plane_3 plane3(p3, n3);

  // Calculate their common intersection.
  const auto result = CGAL::intersection(plane1, plane2, plane3);

  if (!result) {
    return false;
  }

  // We only accept a unique intersection point.
  const Point_3* point = std::get_if<Point_3>(&*result);

  if (point == nullptr) {
    return false;
  }

  // Check whether the intersection lies inside all three splat disks.
  if (CGAL::squared_distance(*point, p1) > r1 * r1 ||
      CGAL::squared_distance(*point, p2) > r2 * r2 ||
      CGAL::squared_distance(*point, p3) > r3 * r3)
  {
    return false;
  }

  intersection_point = *point;
  return true;
}

} // namespace Splat_surface_reconstruction_3
} // namespace internal
} // namespace CGAL

#endif // CGAL_INTERNAL_SPLAT_SURFACE_RECONSTRUCTION_GEOMETRY_UTILS_H
