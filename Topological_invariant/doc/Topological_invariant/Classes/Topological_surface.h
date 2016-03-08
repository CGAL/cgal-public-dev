namespace CGAL {
/*!
\ingroup PkgTopologicalInvariantConcepts
\cgalConcept

\tparam Item
\tparam HalfedgeGraph model of HalfedgeGraph
\tparam Allocator has to match the standard allocator requirements.

*/
template<Item, HalfedgeGraph, Allocator>
class Topological_surface
{
public:
    typedef Map::Halfedge_descriptor Halfedge_descriptor;
    
    ///model of Path
    typedef unspecified_type Path_handle;
    
    ///model of Graph_on_surface
    typedef unspecified_type Graph_handle;
    
    ///model of Arc_occurence
    typedef unspecified_type Arc_occurence_handle;
    
    ///order around a dart
    typedef unspecified_type Arc_occurence_Order_range;
    
    ///Contructor
    Topological_surface(HalfedgeGraph& surface);
    
    ///range of all Arc_occurence of a halfedge from left to right
    Arc_occurence_Order_range Order(Halfedge_descriptor) ;
    
    /// create a new path
    Path_handle create_path(); 
    
    ///erase a path
    void erase_path(Path_handle);
    
    ///create a embeded graph
    Graph_handle create_graph();
    
    ///erase a embeded graph
    void erase_graph(Graph_handle);
}; 
}