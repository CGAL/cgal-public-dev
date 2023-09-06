/*!
\ingroup PkgCollisions3Concepts
\cgalConcept

The concept `CollisionTestBoundary` defines the boundary needed to test for a collision 
between two geometric objects. The test boundary is constructed by defining a parameterization,
referred to as a collision function, over all arbitrary vectors connecting the two 
geometric objects throughout their motion, and then extracting the boundary of the 
collision function's codomain. The collision function evaluates to the origin if and only 
if a zero-length vector connects the two objects at some point during their motion. The 
number of such roots has the same parity as the number of ray-intersections with the 
boundary for any ray emanating from the origin, and an odd number of roots corresponds 
to a collision.


\cgalHasModel `CGAL::Collisions::internal::Segment_3_Segment_3_collision_test_boundary<Kernel>`
\cgalHasModel `CGAL::Collisions::internal::Point_3_Triangle_3_collision_test_boundary<Kernel>`

*/

template <class Kernel>
class CollisionTestBoundary {
public:


/*!
  A functor object with a vector-valued operator that evaluates to 0 if and 
  only if a potential collision occurs. Provides the operator:
  `Point_3 operator()(const FT& t, const FT& u, const FT& v);`
*/
typedef unspecified_type Collision_function;

/*!
  Returns the number of intersections between a given ray and the facets of the
  test boundary. An odd number of intersections indicates a collision has occurred.
*/
std::size_t num_ray_intersections(Kernel::Ray_3 r);

}; /* end CollisionTestBoundary */

