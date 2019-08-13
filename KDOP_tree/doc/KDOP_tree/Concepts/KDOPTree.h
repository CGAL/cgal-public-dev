/*!
 * \ingroup PkgKDOPTreeConcepts
 * \cgalConcept
 *
 * The concept 'KDOPTree' provides functions of a k-dop tree, including
 * building up the tree, traversing the tree, intersection computation,
 * etc.
 *
 * \cgalHasModel 'CGAL::KDOP_tree<KDOPTraits>'
 *
 */

class KDOPTree
{
public:
  /// \name Types
  /// @{

  /// Number type of the geometry kernel
  typedef unspecified_type FT;

  /// Type of 3D point
  typedef unspecified_type Point;

  /// Type of input primitive
  typedef unspecified_type Primitive;

  /// Identifier for a primitive in the tree
  typedef unspecified_type Primitive_id;

  /// Unsigned integer size type
  typedef unspecified_type size_type;

  /// Type of k-dop
  typedef unspecified_type Kdop;

  /// Type of direction
  typedef unspecified_type Direction_type;

  /// 3D point and primitive id type
  typedef unspecified_type Point_and_primitive_id;

  /// An alias to `KDOPTraits::Intersection_and_primitive_id<Query>`
  typedef unspecified_type Intersection_and_primitive_id;

  /// @}

public:
  /// \name Tree build
  /// @{

  /// Build the k-dop tree with a binary tree and compute k-dops in the process.
  unspecified_type build();

  /// Add a primitive to the k-dop tree
  unspecified_type insert(const Primitive& p);

  /// Add a sequence of primitives to the k-dop tree
  unspecified_type insert(unspecified_type first, unspeficified_type beyond);

  /// Set user-defined directions
  unspecified_type set_kdop_directions(const unspecified_type& directions);

  /// @}

public:
  /// \name Intersection tests
  /// @{

  /// Check if a query intersects the k-dop tree of the primitives.
  unspecified_type do_intersect(const unspecified_type& query);

  /// Return the number of intersected primitives by the query.
  unspecified_type num_of_intersected_primitives(const unspecified_type& query);

  /// Return all intersected primitives by the query.
  unspecified_type all_intersected_primitives(const unspecified_type& query, unspecified_type out);

  /// Return the id of the intersected primitive which is encountered first in
  /// the tree traversal.
  unspecified_type any_intersected_primitive(const unspecified_type& query);

  /// Return the id of the first intersected primitive closest to the source
  /// point of the ray query
  unspecified_type first_intersected_primitive(const unspecified_type& query);

  /// @}

  /// \name Intersections
  /// @{

  /// Return the list of all intersections by the query.
  unspecified_type all_intersections(const unspecified_type& query, unspecified_type out);

  /// Return the intersection encountered first in the tree traversal.
  unspecified_type any_intersection(const unspecified_type& query);

  /// Return the first intersection and primitive id closest to the source point
  /// of the ray query.
  unspecified_type first_intersection(const unspecified_type& query, const unspecified_type& skip);

  /// @}

  /// \name Distance queries
  /// @{

  /// Return the point in all input primitives closest to the query.
  unspecified_type closest_point(const unspecified_type& query);

  /// Return the squared distance between the query point and the closest point.
  unspecified_type squared_distance(const unspecified_type& query);

  /// Return the point and the primitive id closest to the query.
  unspecified_type closest_point_and_primitive(const unspecified_type& query);

  /// @}

};
