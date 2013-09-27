/*!
\ingroup PkgExtremePointsDConcepts
\cgalConcept

Requirements of the traits class to be used with the functions `CGAL::extreme_points_d`, `CGAL::extreme_points_d_dula_helgason`, `CGAL::extreme_points_d_simple` and the class `CGAL::Extreme_points_d`.

\cgalHasModel `CGAL::Extreme_points_traits_d<Point>`

*/
class  ExtremePointsTraits_d {
public:

/// The point type on which the extreme point algorithm operates
typedef Hidden_type Point;
/// The `RingNumberType` used for the homogeneous coordinates of the input points
typedef Hidden_type RT;
/// Binary predicate object type comparing `Point` s
/// lexicographically. Must provide `bool operator()(Point p, Point q)` where `true`
/// is returned iff \f$ p \f$ is lexicographically smaller than \f$ q \f$.
typedef Hidden_type Less_lexicographically;

/// Function object type that provides `RandomAccessIterator operator()(Point &p)`, which 
/// returns a random access iterator over the homogeneous coordinates of \f$ p\f$ pointing to the zeroth 
/// homogeneous coordinate \f$ h_0\f$ of \f$ p\f$. The value type of this random access iterator (i.e., the
/// type of the homogeneous coordinates) must be `RT`.
typedef Hidden_type Homogeneous_begin;
}
