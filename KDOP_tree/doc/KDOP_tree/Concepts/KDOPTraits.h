/*!
\ingroup PkgKDOPTreeConcepts
\cgalConcept

The concept `KDOPTraits` provides the geometric primitive types and methods for the class `CGAL::KDOP_tree<KDOPTraits>`.

\cgalHasModel `CGAL::KDOP_traits<KDOPGeomTraits,KDOPPrimitive>`

\cgalRefines `SearchGeomTraits_3`

\sa `CGAL::KDOP_traits<KDOPGeomTraits,KDOPPrimitive>`
\sa `CGAL::KDOP_tree<KDOPTraits>`
\sa `KDOPPrimitive`
\sa `KDOPKdop`

 */
class KDOPTraits {
public:

  /// \name Types
  /// @{

  /// Type of geometry traits (kernel)
  typedef unspecified_type Geom_traits;

  /// Number type of the geometry kernel
  typedef unspecified_type FT;

  /// Type of k-DOP traits
  typedef unspecified_type KT;

  /// Type of primitives
  typedef unspecified_type Primitive;

  /// 3D point and primitive id type
  typedef unspecified_type Point_and_primitive_id;

  /*!
Type of a 3D point.
   */
  typedef unspecified_type Point_3;

  /// Type of a sphere
  typedef unspecified_type Sphere_3;

  /// Type of a bounding box
  typedef unspecified_type Bounding_box;

  /*!
K-dop type.
   */
  typedef unspecified_type Kdop;

  /// Type of support heights of a k-DOP
  typedef unspecified_type Array_height;

  /// @}

  /// \name Splitting
  /// During the construction of the KDOP tree, the primitives are
  /// splitted in the binary way.
  /// @{

  /*!
A functor object to split a range of primitives. Provides the operator:
 `void operator()(InputIterator first, InputIterator beyond);`
 %Iterator type `InputIterator` must be a model of RandomAccessIterator
 and have `Primitive` as value type. The operator is used for determining
 the primitives assigned to the children nodes of a given node,
 assuming that the goal is to split at the centroid of the collection of
 primitives. The primitives assigned to this node are passed as argument
 to the operator.
   */
  typedef unspecified_type Split_primitives;

  /*!
  A functor object to compute the bounding box of a set of primitives. Provides the operator:
  `Bounding_box operator()(Input_iterator begin, Input_iterator beyond);` %Iterator type `InputIterator` must have `Primitive` as value type.
  */
  typedef unspecified_type Compute_bbox;

  /// @}


  /// \name K-dop computation
  /// Different from the implementation of AABB tree, the K-dop computation is
  /// after the whole splitting is completed so that the K-dop can be computed
  /// recursively in a bottom-top manner, without duplicate computation as in
  /// AABB tree where the bounding box is computed after each splitting.
  /// @{
  /*!
A functor object to compute the kdop (support heights) of a set of primitives.
Provides the operator: `Kdop operator()(Input_iterator begin, Input_iterator beyond);`
 %Iterator type `InputIterator` must have `Primitive` as value type. The operation
 contains the union of k-dops of the children of a node in the tree.
   */
  typedef unspecified_type Compute_kdop;

  /// @}

  /// \name Intersections
  /// @{

  /*!
A functor object to compute intersection predicates between the query and the nodes of the tree. Provides the operators:
- `bool operator()(const Query & q, const Kdop & kdop, const Array_height & support_heights);` which returns `true` iff the query intersects the kdop
- `bool operator()(const Query & q, const Primitive & primitive);` which returns `true` iff the query intersects the primitive
   */
  typedef unspecified_type Do_intersect;

  /*!
A functor object to compute the intersection of a query and a primitive. Provides the operator:
`boost::optional<Intersection_and_primitive_id<Query>::%Type > operator()(const Query & q, const Primitive& primitive);` which returns the intersection as a pair composed of an object and a primitive id, iff the query intersects the primitive.
   */
  typedef unspecified_type Intersection;

  /// @}

  /// \name Distance queries
  /// @{

  /// A functor object to compute closest point from the query on a primitive. Provides the operator:
  /// `Point_3 operator()(const Query& query, const Primitive& primitive, const Point_3 & closest);` which returns the closest point to `query`, among `closest` and all points of the primitive.
  typedef unspecified_type Closest_point;

  /// A functor object to compare distance between the point and the nodes of the tree. Provides the operators:
  ///- `bool operator()(const Query & query, const Kdop & kdop, const Array_height & support_heights, const Point & closest);` which returns `true` iff the k-DOP is closer to `query` than `closest` is
  typedef unspecified_type Compare_distance;

  /// @}

  /// \name Operations
  /// @{

  /*!
Return the primitive splitting functor.
   */
  Split_primitives split_primitives_object();

  /// Return the bounding box functor.
  Compute_bbox compute_bbox_object();

  /*!
Return the support heights constructor.
   */
  Compute_kdop compute_kdop_object();

  /*!
Return the intersection detection functor.
   */
  Do_intersect do_intersect_object();

  /*!
Return the intersection constructor.
   */
  Intersection intersection_object();

  /// Return the distance comparison functor
  Compare_distance compare_distance_object();

  /// Return the closest point computation functor
  Closest_point closest_point_object();

  /// @}


}; /* end KDOPTraits */
