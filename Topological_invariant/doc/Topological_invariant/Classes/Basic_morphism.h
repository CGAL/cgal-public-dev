namespace CGAL
{
	/*!
	 * \ingroup PkgTopologicalInvariantClasses
	 * \cgalModels `::Morphism`
	 * 
	 * `Basic_morphism` define a morphism between to `Topological_surface`.
	 * 
	 * \tparam Source_surface must be a model of the `Topological_surface` concept.
	 * \tparam Target_surface must be a model of the `Topological_surface` concept.
	 * \tparam Alloc has to match the standard allocator requirements. The `rebind` mechanism  `Alloc` will be used to create appropriate allocators internally with value type `Dart`.the default value is CGAL_ALLOCATOR(int)` from the `<CGAL/memory.h>` header file.
	 * 
	 */
	template< typename Source_generalized_graph , typename Target_generalized_graph, typename Alloc >
	class Basic_morphism{
	public:
		/// \name Construction
		/// @{
	
		/*!
		* construct a new morphism with gsource as source and gtarget as target.
		*/
		Basic_morphism(Source_generalized_graph gsource, Target_generalized_graph gtarget);
		
		/// @}
		/// \name Access Member Functions
		/// @{
		
		/*!
		* return true iff a Morphism exist for all dart of source.
		*/
		bool is_complete() const;
		
		/*!
		* Returns true iff the Morphism is valid.
		* 
		* a morphisme \f$ M \f$ is valid if for all dart \f$d_s \in source \f$ and \f$d_t \in target\f$ with \f$ d_t = M(d_s)\f$ : \f$ \alpha_i(M(d_s)) = M(\alpha_id_t)) \f$
		*/
		bool 	is_valid () const;
		
		/// @}
		/// \name Modifier
		///@{
		
		/*!
		* define the morphism \f$ Morphism(source) = target \f$.
		*/
		void add(Source_halfedge_handle source, Target_halfedge_handle target);
		/// @}
	};
}
