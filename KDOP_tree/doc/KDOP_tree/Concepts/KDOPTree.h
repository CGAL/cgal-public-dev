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
 * \sa `CGAL::KDOP_tree<KDOPTraits>`
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

  /// 3D point and primitive id type
  typedef unspecified_type Point_and_primitive_id;

  /// An alias to `KDOPTraits::Intersection_and_primitive_id<Query>`
  struct Intersection_and_primitive_id {
    typedef unspecified_type Type;
  };

  /// @}

public:
  /// \name Tree build
  /// @{

  /// Build the k-dop tree with a binary tree (or an octree)
  /// without computing k-dops in the process.
  /// \todo Different from AABB tree structure, the k-dops are computed in the process of tree build.
  unspecified_type build();

  /// Add a primitive to the k-dop tree
  unspecified_type insert(const Primitive& p);

  /// Add a sequence of primitives to the k-dop tree
  unspecified_type insert(unspecified_type first, unspeficified_type beyond);

  /// Compute k-dops of the tree
  /// \todo The k-dops are computed after build() by traversing the tree.
  unspecified_type compute_kdop(unspecified_type first, unspecified_type beyond);

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

  /// Return the first intersected primitive id closest to the source point of
  /// the ray query
  unspecified_type first_intersection_primitive(const unspecified_type& query, const unspecified_type& skip);

  /// @}

};
