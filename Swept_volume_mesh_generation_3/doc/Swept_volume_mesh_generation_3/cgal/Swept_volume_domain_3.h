/*!
\ingroup PkgSVMG_3Classes


The class `Swept_volume_domain_3` is a model of the concept 'MeshDomain_3', it 
describes the knowledge required to generate a mesh approximating a swept volume,
given a generator and trajectory.

Using the provided meshing criteria 'Swept_volume_facet_criteria_3'
the generated mesh is guaranteed to be conservative (i.e., does not intersect the actual)
swept volume and obeys a user defined a-priori geometric bound in terms of the one-sided Hausdorff distance.

*/

template<typename GeometryTraits>
class Swept_volume_domain_3 {
public:

/// \name Types required by concept
/// @{

/*!
Geometric traits class. This type is defined to ensure compatibility with
`CGAL::Kernel_traits<T>`.
*/
typedef GeometryTraits R;

/*!
Point type.
*/
typedef GeometryTraits::Point_3 Point_3;

/*!
Segment type.
*/
typedef GeometryTraits::Segment_3 Segment_3;

/*!
Ray type.
*/
typedef GeometryTraits::Ray_3 Ray_3;

/*!
Line type.
*/
typedef GeometryTraits::Line_3 Line_3;

/*!
A type to distinguish
`MeshDomain_3` models from `MeshDomainWithFeatures_3` models.
*/
typedef CGAL::Tag_false Has_features;

/*!
Type of indices for subdomains of the
input domain. Must be a model of CopyConstructible,
Assignable, DefaultConstructible and EqualityComparable.
The default constructed value must match the label of the exterior of
the domain (which contains at least the unbounded component).
*/
typedef unspecified_type Subdomain_index;

/*!
Type of indices for surface patches
(boundaries and interfaces) of the
input domain. Must be a model of CopyConstructible,
Assignable, DefaultConstructible and EqualityComparable.
The default constructed value must be the index value assigned
to a non surface facet.
*/
typedef unspecified_type Surface_patch_index;

/*!
Type of indices to be stored at mesh vertices
to characterize the lowest dimensional face of the input complex
on which the vertex lies. Must be a model of CopyConstructible,
Assignable, DefaultConstructible and EqualityComparable.

*/
typedef unspecified_type Index;

/*!
Return type of `Construct_intersection` queries.
`int` represents the
dimension of the lower dimensional face of the input complex on which the intersection
point lies and `Index` is the index of this face.
*/
typedef CGAL::cpp11::tuple<Point_3, Index, int> Intersection;

/*!
A function object to construct
a set of initial points on the surface of the domain. Provides the
following operators:

`template<typename OutputIterator>`
<br>
`OutputIterator operator()(OutputIterator pts)`

`template<typename OutputIterator>`
<br>
`OutputIterator operator()(int n, OutputIterator pts)`

Those two operators output a set of (`n`) surface points to the
output iterator `pts`, as objects of type `std::pair<Point_3,
Index>`. If `n` is not given, the functor must provide enough
points to initialize the mesh generation process.
*/
typedef unspecified_type Construct_initial_points;

/*!
A function object to query whether a point is in
the input domain or not. In the positive case, it outputs the
subdomain which includes the query point. Provides the operator:

`boost::optional<Subdomain_index> operator()(Point_3 p)`
*/
typedef unspecified_type Is_in_domain;

/*!
A function object which answers
intersection queries between the surface patches of the domain and
objects of type `Segment_3`, `Ray_3` or
`Line_3`. Provides the operators:

`boost::optional<Surface_patch_index> operator()(Segment_3 s)`

`boost::optional<Surface_patch_index> operator()(Ray_3 r)`

`boost::optional<Surface_patch_index> operator()(Line_3 l)`

The return type of the operators tell whether or not the query intersects a
surface patch. In the positive case, it provides (through operator*()) the
`Surface_patch_index` of one of the intersected surface patches.
*/
typedef unspecified_type Do_intersect_surface;

/*!
A function object to construct the
intersection between an object of type `Segment_3`, `Ray_3` or
`Line_3` and an interface. Provides the operators:

`Intersection operator()(Segment_3 s)`

`Intersection operator()(Ray_3 r)`

`Intersection operator()(Line_3 l)`
\pre do_intersect_surface(s/r/l) == true
*/
typedef unspecified_type Construct_intersection;

/// @}

/// \name Additional Types
/// @{

/*!
Affine Transformation Type
*/
typedef  GeometryTraits::Aff_transformation_3 Aff_transformation_3;

///*!
//3D Polyhedral Surfaces Type
//*/
//typedef   CGAL::Polyhedron_3 Polyhedron_3;


/*!
A model of `MeshCriteria_3`. 
*/
template< typename MeshCriteria_3>
class Swept_volume_criteria_3; 

/*!
The returned Swept_volume_criteria_3 object first applies the criteria 
that are given. In case the given entety (face or cell) 
is not already classified as bad, some additional criteria are 
applied. 
This ensures that  the generated mesh is conservative, i.e., 
the swept volume is enclosed by the output mesh, and that the 
one-sided Hausdorff distance of the generated mesh to the swept volume is 
upper bounded by the user defined tolerance (already given in the 
constructor of 'Swept_volume_domain_3').
*/
template< typename MeshCriteria_3>
Swept_volume_criteria_3<MeshCriteria_3> 
swept_volume_criteria_3_object(const MeshCriteria_3& criteria);


/// @}


/// \name Constructors
/// The following are constructors
/// @{
/*!
The swept object (generator) must be given as an indexed face set, that is, a range of vertices and a range of triples.
Each triple defines the indices of one triangle, the indices reference to the range of vertices.

The trajectory is expected to be a sequence of rigid body transformations and has to be a 'Range' of 'Aff_transformation_3'. 
The bound \f$ \epsilon\f$ determines the geometric fidelity of the final swept volume, the one-sided Hausdorff
error between the actual SV and its approximation is guaranteed to be smaller than \f$ \epsilon \f$.

With downstepping enabled, firstly a coarser approximation is computed, and then refined. This happens without loss of geometric 
guarantees. It is a trade off between running time and memory consumption.
*/

public Swept_volume_domain_3(range<Point_3> vertices, range< CGAL::cpp11::triple<int, int, int> > indices, range<Aff_transformation_3> trajectory, double epsilon, bool downstep = false);





/// @}



/// \name Operations
/// The following functions give access to the function objects:
/// @{

/*!

*/
Construct_initial_points construct_initial_points_object();

/*!

*/
Is_in_domain is_in_domain_object();

/*!

*/
Do_intersect_surface do_intersect_surface_object();

/*!

*/
Construct_intersection construct_intersection_object();

/*!
Returns in vertices and indices the computed mesh appoximating the swept volume.
*/
void getIndexedFaceSet(range<Point_3> &vertices, range< CGAL::cpp11::triple<int, int, int> > &indices);




/// @}

/// \name Index Conversion
/// These methods are designed to convert indices:
/// @{

/*!
Returns
the index to be stored at a vertex lying on the surface patch identified by `surface_patch_index`.
*/
Index index_from_surface_patch_index(Surface_patch_index surface_patch_index);

/*!
Returns
the index to be stored at a vertex lying in the subdomain identified by `subdomain_index`.
*/
Index index_from_subdomain_index(Subdomain_index subdomain_index);

/*!
Returns the `Surface_patch_index` of the surface patch
where lies a vertex with dimension 2 and index `index`.
*/
Surface_patch_index surface_patch_index(Index index);

/*!
Returns the index
of the subdomain containing a vertex with dimension 3 and index `index`.
*/
Subdomain_index subdomain_index(Index index);

/// @}

}; /* end MeshDomain_3 */
