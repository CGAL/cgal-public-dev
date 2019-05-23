/*!
\ingroup PkgKDOPTreeConcepts
\cgalConcept

The concept `KDOPPrimitive` describes the requirements for the primitives stored in the KDOP tree data structure. The concept encapsulates a type for the input datum (a geometric object) and an identifier (id) type through which those primitives are referred to. The concept `KDOPPrimitive` also refines the concepts DefaultConstructible and Assignable.

\sa `CGAL::KDOP_tree<KDOPTraits>`
\sa `KDOPPrimitiveWithSharedData`

\cgalHeading{Example}

The `Primitive` type can be, e.g., a wrapper around a `Handle`. Assume for instance that the input objects are the triangle faces of a mesh stored as a `CGAL::Polyhedron_3`. The `Datum` would be a `Triangle_3` and the `Id` would be a polyhedron `Face_handle`. Method `datum()` can return either a `Triangle_3` constructed on the fly from the face handle or a `Triangle_3` stored internally. This provides a way for the user to trade memory for efficiency.

\cgalHasModel `CGAL::KDOP_primitive<Id,ObjectPropertyMap,PointPropertyMap,Tag_false,CacheDatum>`
\cgalHasModel `CGAL::KDOP_segment_primitive<Iterator,CacheDatum>`
\cgalHasModel `CGAL::KDOP_triangle_primitive<Iterator,CacheDatum>`
\cgalHasModel `CGAL::KDOP_halfedge_graph_segment_primitive<HalfedgeGraph,VertexPointPMap,Tag_false,CacheDatum>`
\cgalHasModel `CGAL::KDOP_face_graph_triangle_primitive<FaceGraph,VertexPointPMap,Tag_false,CacheDatum>`
*/

class KDOPPrimitive {
public:

/// \name Types
/// @{

/*!
3D point type.
*/
typedef unspecified_type Point;

/*!
Type of input datum.
*/
typedef unspecified_type Datum;

/*!
Point reference type returned by the function `point()`. It is convertible to the type `Point`.
 */
typedef unspecified_type Point_reference;

/*!
Datum reference type returned by the function `datum()`. It is convertible to the type `Datum`.
*/
typedef unspecified_type Datum_reference;

/*!
Type of identifiers through which the input objects are referred to. It must be a model of the concepts DefaultConstructible and Assignable.
*/
typedef unspecified_type Id;

/// @}

/// \name Operations
/// @{

/*!
Returns the datum (geometric object) represented by the primitive.
*/
Datum_reference datum();

/*!
Returns the corresponding identifier. This identifier is only used as a reference for the objects in the output of the `KDOP_tree` methods.
*/
Id id();

/*!
Returns a 3D point located on the geometric object represented by the primitive. This function is used to sort the primitives during the KDOP tree construction as well as to construct the search KD-tree internal to the KDOP tree used to accelerate distance queries.
*/
Point_reference reference_point();

/// @}

}; /* end KDOPPrimitive */
