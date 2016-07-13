/*!
 * \ingroup PkgTopologicalInvariantConcepts
 * \cgalConcept
 * 
 * the concept `Morphism` define a morphism from a `Topological_surface` source to a `Topological_surface` target.
*/
class Morphism{
public:
	/// \name Source Type
	/// @{
	
	typedef unspecified_type Source_generalized_graph;
	typedef unspecified_type Source_halfedge_handle;
	/// @}
	
	/// \name Target Type
	/// @{
	typedef unspecified_type Target_generalized_graph;
	typedef unspecified_type Target_halfedge_handle;
	/// @}
	
	/// \name Access Member Functions
	/// @{
	
	/*!
	 * return the morphism of a dart.
	 */
	Target_halfedge_handle operator() (Source_halfedge_handle source) const;
	
	/*!
	 * return the generalized map source.
	 */
	Source_generalized_graph& source();
	
	/*!
	 * return the generalized map target.
	 */
	Target_generalized_graph& target();
	
	/// @}
};
