/*!
 * \ingroup PkgTopologicalInvariantConcepts
 * \cgalConcept
 * 
 * the concept `Morphism` define a morphism from a `GeneralizedMap` source to a `GeneralizedMap` target.
*/
class Morphism{
public:
	/// \name Source Type
	/// @{
	
	typedef unspecified_type Source_generalized_graph;
	typedef unspecified_type Source_halfedge_descriptor;
	/// @}
	
	/// \name Target Type
	/// @{
	typedef unspecified_type Target_generalized_graph;
	typedef unspecified_type Target_halfedge_descriptor;
	/// @}
	
	/// \name Access Member Functions
	/// @{
	
	/*!
	 * return the morphism of a dart.
	 */
	Target_halfedge_descriptor operator() (Source_halfedge_descriptor source) const;
	
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
