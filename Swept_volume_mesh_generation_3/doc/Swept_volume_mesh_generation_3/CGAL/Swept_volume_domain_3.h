namespace CGAL{
/*!
\ingroup PkgSVMG_3

\brief The class `Swept_volume_domain_3` is a model of the concept 
`MeshDomain_3`, it describes the knowledge required to generate a 
mesh approximating a swept volume, given a generator and trajectory.

Using the provided meshing criteria (see below) the generated mesh is 
guaranteed to be conservative (i.e., does not intersect the actual)
swept volume and obeys a user defined a-priori geometric bound in terms 
of the one-sided Hausdorff distance.
*/

template<typename GeometryTraits>
class Swept_volume_domain_3{
public:

/// \name Types required by concept
/// @{

/*!
`Swept_volume_domain_3` is not a model of `MeshDomainWithFeatures_3`.
*/
typedef CGAL::Tag_false Has_features;

/*!
Type used to index subdomains. 
*/
typedef int Subdomain_index;

/*!
Type to index surface patches. 
*/
typedef std::pair<int> Surface_patch_index;

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
typedef   CGAL::Polyhedron_3<GeometryTraits> Polyhedron_3;


/*!
A model of `MeshCriteria_3`. 
*/
template< typename MeshCriteria_3>
class Swept_volume_criteria_3; 

/*!
The returned Swept_volume_criteria_3 object first applies the criteria 
that are given. In case the given entity (face or cell) 
is not already classified as bad, the additional criteria discussed 
in cgal:SV-vDHS-2012 are applied. 
This ensures that the generated mesh is conservative, i.e., 
the swept volume is enclosed by the output mesh, and that the 
one-sided Hausdorff distance of the generated mesh to the swept volume is 
upper bounded by the user defined tolerance (already given in the 
constructor of `Swept_volume_domain_3`).

The termination of the refinement process only depends on the criteria 
as the additional criteria do not apply as soon as the current mesh 
is within the tolerance area. 
*/
template< typename MeshCriteria_3>
Swept_volume_criteria_3<MeshCriteria_3> 
swept_volume_criteria_3_object(const MeshCriteria_3& criteria);


/// @}


/// \name Constructors
/// The following are constructors
/// @{
/*!
\brief Construction from indexed face set and trajectory. 

The swept object (generator) must be given as an indexed face set, that is, a range of vertices and a range of 3-tuples.
Each tuple defines the indices of one triangle, the indices reference to the range of vertices.

The trajectory is expected to be a sequence of rigid body transformations and has to be a `Range` of `Aff_transformation_3`. 
The bound \f$ \epsilon\f$ determines the geometric fidelity of the final swept volume, the one-sided Hausdorff
error between the actual SV and its approximation is guaranteed to be smaller than \f$ \epsilon \f$.

With downstepping enabled, firstly a coarser approximation is computed, and then refined. This happens without loss of geometric 
guarantees. It is a trade off between running time and memory consumption.
*/

public Swept_volume_domain_3(range<Point_3> vertices, range< CGAL::cpp11::tuple<int, int, int> > indices, range<Aff_transformation_3> trajectory, double epsilon, bool downstep = false);


/*!
\brief Construction from polyhedron and trajectory. 

The swept object (generator) is given by a polyhedron with triangular facets.

The trajectory is expected to be a sequence of rigid body transformations and has to be a `Range` of `Aff_transformation_3`. 
The bound \f$ \epsilon\f$ determines the geometric fidelity of the final swept volume, the one-sided Hausdorff
error between the actual SV and its approximation is guaranteed to be smaller than \f$ \epsilon \f$.

With downstepping enabled, firstly a coarser approximation is computed, and then refined. This happens without loss of geometric 
guarantees. It is a trade off between running time and memory consumption.

\pre polyhedron.is_pure_triangle()
*/

public Swept_volume_domain_3(const Polyhedron& polyhedron, range<Aff_transformation_3> trajectory, double epsilon, bool downstep = false);






/// @}




}; /* end Swept_volume_domain_3 */

}
