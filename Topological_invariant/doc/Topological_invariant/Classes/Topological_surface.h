namespace CGAL {
/*!
\ingroup PkgTopologicalInvariantConcepts

\cgalModels `Generalized_map`

\tparam Item
\tparam GMap model of Generalized_map
\tparam Allocator has to match the standard allocator requirements.
*/



template<Item, GMap, Allocator>
class Topological_surface
{
public:
    /// \name Constructor
    /// @{
    TopoSurface(const GMap& gmap);
    /// @}
    
    /// \name Types
    /// @{
    typedef unspecified_type Halfedge;
    typedef unspecified_type Halfedge_handle;
    typedef unspecified_type Halfedge_const_handle;
    /// @}
    
    /// \name Path and embeded graph types
    /// @{
    
    ///model of Path
    typedef unspecified_type Path_handle;
    
    ///model of Graph_on_surface
    typedef unspecified_type Graph_handle;
    
    ///model of Arc_occurence
    typedef unspecified_type Arc_occurence_handle;
    
    ///order around a dart
    typedef unspecified_type Arc_occurence_Order_range;
    
    /// @}
    
    /// \name Range types
    /// @{
    typedef unspecified_type Halfedge_range;
    template<unsigned int i>
    using Halfedge_of_cell = unspecified_type;
    template<unsigned int i>
    using One_halfedge_per_cell = unspecified_type;
    /// @}
    
    /// \name Accesss
    /// @{
    Halfedge_handle vertex_next(Halfedge_handle he);
    Halfedge_handle opposite();
    Halfedge_handle signature();
    Dart_handle Dart(Halfedge_handle, bool);
    Halfedge_handle halfedge(Dart);
    bool sign(Dart_handle);
    Halfedge_handle face_next(Halfedge_handle he);
    /// @}
    
    /// \name Range access
    /// @{
    Halfedge_range halfedges();
    template<unsigned int i>
    Halfedge_of_cell<i> halfedge_of_cell();
    template<unsigned int i>
    One_halfedge_per_cell<i> one_halfedge_per_cell();
    /// @}
    
    /// \name Operation
    /// @{
    Halfedge_handle edge_contraction(Halfedge_handle he);
    Halfedge_handle edge_deletrion(Halfedge_handle he)
    Halfedge_handle edge_sudivision(Halfedge_handle he);
    Halfedge_handle face_subdivision(Halfedge_handle a, Halfedge_handle b);
    void cut(Path_handle path);
    void flip(Halfedge_handle he);
    /// @}
    
    /// \name Path and embeded graph
    /// @{
    
    ///range of all Arc_occurence of a halfedge from left to right
    Arc_occurence_Order_range Order(Halfedge_handle) ;
    
    /// create a new path
    Path_handle create_path(); 
    
    ///erase a path
    void erase_path(Path_handle);
    
    ///create a embeded graph
    Graph_handle create_graph();
    
    ///erase a embeded graph
    void erase_graph(Graph_handle);
    
    ///@}
}; 
}