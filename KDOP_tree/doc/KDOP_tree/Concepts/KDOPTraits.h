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

/*!
Type of a 3D point.
*/
typedef unspecified_type Point_3;

/*!
Type of primitive.
Must be a model of the concepts `KDOPPrimitive` or `KDOPPrimitiveWithSharedData`.
*/
typedef unspecified_type Primitive;

/*!
K-dop type.
*/
typedef unspecified_type Kdop;

/*!
3D Point and Primitive Id type
*/
typedef std::pair<Point_3, Primitive::Id> Point_and_primitive_id;

/// @}

/// \name Splitting
/// During the construction of the KDOP tree, the primitives are
/// splitted in some way.
/// \todo Split into an octree or a binary tree.
/// \todo Splitting scheme based on the centroid of the collection of primitives.
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

/// @}


/// \name K-dop computation
/// Different from the implementation of AABB tree, the K-dop computation is
/// after the whole splitting is completed so that the K-dop can be computed
/// recursively in a bottom-top manner, without duplicate computation as in
/// AABB tree where the bounding box is computed after each splitting.
/// \todo A recursive approach to computing K-dop tree.
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
- `bool operator()(const Query & q, const Kdop & kdop);` which returns `true` iff the query intersects the kdop
- `bool operator()(const Query & q, const Primitive & primitive);` which returns `true` iff the query intersects the primitive
*/
typedef unspecified_type Do_intersect;

/*!
A functor object to compute the intersection of a query and a primitive. Provides the operator:
`boost::optional<Intersection_and_primitive_id<Query>::%Type > operator()(const Query & q, const Primitive& primitive);` which returns the intersection as a pair composed of an object and a primitive id, iff the query intersects the primitive.
*/
typedef unspecified_type Intersection;

/*!
A functor object to compare two points. Provides the operator:
`bool operator()(const Point_3& p, const Point_3& q);}` which returns `true` if `p` is equal to `q`.
*/
typedef unspecified_type Equal_3;
/// @}

/// \name Operations
/// @{

/*!
Returns the primitive splitting functor.
*/
Split_primitives split_primitives_object();

/*!
Returns the support heights constructor.
*/
Compute_kdop compute_kdop_object();

/*!
Returns the intersection detection functor.
*/
Do_intersect do_intersect_object();

/*!
Returns the intersection constructor.
*/
Intersection intersection_object();

/*!
Returns the equal functor.
*/
Equal_3 equal_3_object();

/// @}


}; /* end KDOPTraits */
