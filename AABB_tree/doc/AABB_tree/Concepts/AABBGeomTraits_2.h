
/*!
\ingroup PkgAABB_treeConcepts
\cgalConcept

The concept `AABBGeomTraits_2` defines the requirements for the first template parameter of the class `CGAL::AABB_traits<AABBGeomTraits, AABBPrimitive>`. It provides predicates and constructors to detect and compute intersections between query objects and the primitives stored in the AABB tree. In addition, it contains predicates and constructors to compute distances between a point query and the primitives stored in the AABB tree.

\cgalHasModel Any Kernel is a model of this traits concept.

\sa `CGAL::AABB_traits<AABBGeomTraits,AABBPrimitive>`

*/

class AABBGeomTraits_2 {
public:

/// \name Types
/// @{

/*!
Circle type, that should be consistent with the distance function chosen for the distance queries, namely the `Squared_distance_2` functor.
*/
typedef unspecified_type Circle_2;

/*!
Point type.
*/
typedef unspecified_type Point_2;

/*!
A functor object to detect intersections between two geometric objects.
Provides the operators:
`bool operator()(const Type_1& type_1, const Type_2& type_2);`
where `Type_1` and `Type_2` are relevant types
among `Ray_2`, `Segment_2`, `Line_2`, `Triangle_2`, `Plane_2` and `Bbox_2`. Relevant herein means that a line primitive (ray, segment, line) is tested against a planar or solid primitive (plane, triangle, box), and a solid primitive is tested against another solid primitive (box against box). The operator returns `true` iff `type_1` and `type_2` have a non empty intersection.
*/
typedef unspecified_type Do_intersect_2;

/*!
A functor object to construct the intersection between two geometric objects.
This functor must support the result_of protocol, that is the return
type of the `operator()(A, B)` is `CGAL::cpp11::result<Intersect_2(A,B)>`.

Provides the operators:
`CGAL::cpp11::result<Intersect_2(A,B)> operator()(const A& a, const B& b);`
where `A` and `B` are any relevant types among `Ray_2`, `Segment_2`, `Line_2`,
`Triangle_2`, `Plane_2` and `Bbox_2`.
Relevant herein means that a line primitive (ray, segment, line) is tested
against a planar or solid primitive (plane, triangle, box).
A model of `Kernel::Intersect_2` fulfills those requirements.
*/
typedef unspecified_type Intersect_2;

/*!
A functor object to construct the circle centered at one point and passing through another one. Provides the operator:
`Circle_2 operator()(const Point_2& p, const Point_2 & q);` which returns the circle centered at `p` and passing through `q`.
*/
typedef unspecified_type Construct_circle_2;

/*!
A functor object to compute the point on a geometric primitive which is closest from a query. Provides the operator:
`Point_2 operator()(const Point_2& p, const Type_2& type_2);` where `Type_2` is any type among `Segment_2` and `Triangle_2`. The operator returns the point on `type_2` which is closest to `p`.
*/
typedef unspecified_type Compute_closest_point_2;

/*!
A functor object to detect if a point lies inside a circle or not.
Provides the operator:
`bool operator()(const Circle_2& s, const Point_2& p);` which returns `true` iff the closed volume bounded by `s` contains `p`.
*/
typedef unspecified_type Has_on_bounded_side_2;

/*!
A functor object to compute the squared radius of a circle. Provides the operator:
`FT operator()(const Circle_2& s);` which returns the squared radius of `s`.
*/
typedef unspecified_type Compute_squared_radius_2;

/*!
A functor object to compute the squared distance between two points. Provides the operator:
`FT operator()(const Point_2& p, const Point_2& q);}` which returns the squared distance between \a p and \a q.
*/
typedef unspecified_type Compute_squared_distance_2;


/// @}

/// \name Operations
/// @{

/*!
Returns the intersection detection functor.
*/
Do_intersect_2 do_intersect_2_object();

/*!
Returns the intersection constructor.
*/
Intersect_2 intersect_2_object();

/*!
Returns the distance comparison functor.
*/
Construct_circle_2 construct_circle_2_object();

/*!
Returns the closest point constructor.
*/
Compute_closest_point_2 compute_closest_point_2_object();

/*!
Returns the closest point constructor.
*/
Has_on_bounded_side_2 has_on_bounded_side_2_object();

/*!
Returns the squared radius functor.
*/
Compute_squared_radius_2 compute_squared_radius_2_object();

/*!
Returns the squared distance functor.
*/
Compute_squared_distance_2 compute_squared_distance_2_object();

/// @}

}; /* end AABBGeomTraits */

