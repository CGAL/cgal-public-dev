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
	
	/*!
	%Source_generalized_map type, a model of the `GeneralizedMap` concept.
	*/
	typedef unspecified_type Source_generalized_map;
	
	/*!
	%Source_dart_handle handle type, equal to `Source_generalized_map::Dart_handle`.
	*/
	typedef Source_generalized_map::Dart_handle Source_dart_handle;
	/// @}
	
	/// \name Target Type
	/// @{
	/*!
	%Target_generalized_map type, a model of the `GeneralizedMap` concept.
	*/
	typedef unspecified_type Target_generalized_map;
	
	/*!
	%Target_dart_handle handle type, equal to `Target_generalized_map::Dart_handle`.
	*/
	typedef Target_generalized_map::Dart_handle Target_dart_handle;
	/// @}
	
	/// \name Access Member Functions
	/// @{
	
	/*!
	 * return the morphism of a dart.
	 */
	Target_dart_handle operator() (Source_dart_handle source) const;
	
	/*!
	 * return the generalized map source.
	 */
	Source_generalized_map& source();
	
	/*!
	 * return the generalized map target.
	 */
	Target_generalized_map& target();
	
	/// @}
};
