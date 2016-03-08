
/*!
\ingroup PkgTopologicalInvariantConcepts
\cgalConcept

\cgalHasModel \ref  CGAL::Arc_occurence "CGAL::Arc_occurence<Generalized_surface>"

*/
class ArcOccurence
{
public:
	/// a handle to a Path.
	typedef unspecified_type Path_handle;
	
	/// halfedge_descriptor type of a GeneralizedSurfaceGraph
	typedef unspecified_type halfedge_descriptor;
	
	/// Arc_occurence handle type. 
	typedef unspecified_type Arc_occurence_handle;
	
	/// return the path.
	Path_handle path();
	
	/// return the next ArcOccurence.
	Arc_occurence_handle next();
	
	/// return the previous ArcOccurence.
	Arc_occurence_handle previous();
	
	/// return the associated halfedge.
	halfedge_descriptor halfedge();
    
     /// return the left occurrence in the order around the arc. 
    Arc_occurence left();
    
    /// return the right occurrence in the order around the arc. 
    Arc_occurence right();
    
    /// left move of a in the order around the arc. 
    Arc_occurence move_left();
    
    /// right move of a in the order around the arc. 
    Arc_occurence move_right();
	
private:
	
};