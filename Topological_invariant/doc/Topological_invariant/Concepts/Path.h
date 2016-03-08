/*!
\ingroup PkgTopologicalInvariantConcepts
\cgalConcept
*/
class Path
{
public:
	/// \name Type
	/// @{
	
	/*!
     * Arc_occurence is a part of a path. Each path is associated with a dart
     * Arc_occurence is a model of ArcOccurence.
     */
	typedef unspecified_type Arc_occurence;
	
	/// Arc_occurence handle type, equal to Arc_occurence::Arc_occurence_handle. 
	typedef unspecified_type Arc_occurence_handle;
	
	/// Dart_handle of a GeneralizedMapGraph.
	typedef unspecified_type Dart_handle;
	
	/// Range of all Arc_occurence of the path.
	typedef unspecified_type Arc_occurence_range;
	
	/// Range of all Arc_occurence of the path in a revert order.
	typedef unspecified_type Arc_occurence_revert_range
	
	/// @}
	/// \name Access
	/// @{
	
	///return the number of Arc_occurence in the path.
	size_type number_of_arc_occurence() const;
	
	/// Return true iff the path is a loop. The source vertex is also the destination vertex. 
	bool is_loop() const;
	
	/// @}
	/// \name Range 
	/// @{
	
	/// return a range of all Arc_occurence in the path.
	Arc_occurence_range arc_occurences();
	
	/// return a range of  all Arc_occurence in the path in a revert order.
	Arc_occurence_revert_range arc_occurences_revert();
	
	/// @}
	/// \name Modifier
	/// @{
	
	/// add a copy of all Arc_occurence of the path p at the begining of this path. 
	Arc_occurence push_front(const Path& p);
	
	/// add a new Arc_occurence a the begining of the path
	Arc_occurence push_front(Dart_handle d);
	
	/// add a copy of all Arc_occurence of the path p at the end of this path. 
	Arc_occurence push_back(const Path& p);
	
	/// add a new Arc_occurence a the end of the path
	Arc_occurence push_back(Dart_handle d);
	
	/// remove the first Arc_occurence.
	void pop_front();
	
	/// remove the last Arc_occurence.
	void pop_back();
	
	/// @}
	/// \name opperation
	/// @{
	
	/*! 
	\cgalFigureBegin{figPath2,path2.png}
	A example of result.
	\cgalFigureEnd
	 */
	Arc_occurence detour_face(Arc_occurence_handle a);	
	
	/// @}
};
