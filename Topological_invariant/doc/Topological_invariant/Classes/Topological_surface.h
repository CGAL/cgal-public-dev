namespace CGAL {
/*!
\ingroup PkgTopologicalInvariantConcepts

\cgalModels `Generalized_map`

\tparam Item model,of GeneralizedMapItems
\tparam Allocator has to match the standard allocator requirements.
\tparam GMap model of Generalized_map
*/



template<Item, Allocator, GMap>
class Topological_surface
{
public:
    /// \name Constructor
    /// @{
    
    /// Default constructor: Creates an empty surface.
    Topological_surface();
    
    ///Construct a new surface from another one.
    Topological_surface(const GMap& gmap);
    
    /// @}
    /// \name Types
    /// @{
    
    /// Halfedge type, a model of the TopologicalHalfedge concept. 
    typedef unspecified_type Halfedge;
    
    /// Halfedge handle type, equal to Halfedge::Halfedge_handle. 
    typedef unspecified_type Halfedge_handle;
    
    /// Halfedge const handle type, equal to Halfedge::Halfedge_const_handle. 
    typedef unspecified_type Halfedge_const_handle;
    
    /// @}
    /// \name Path and embedded graph types
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
    
    /*!Range of all the halfedges of the surface.
     * 
     * This type is a model of Range concept, its iterator type is bidirectional and its value type is Halfedge. */
    typedef unspecified_type Halfedge_range;
    
    /// Range of all the halfedges of an <I>i</I>-cell. This type is a model of `Range` concept, its iterator type is forward and its value type is Halfedge. 
    template<unsigned int i>
    using Halfedge_of_cell_range = unspecified_type;
    
    /// Range of one halfedge of each i-cell of the surface. 
    template<unsigned int i>
    using One_halfedge_per_cell = unspecified_type;
    
    /// @}
    /// \name Accesss
    /// @{
    
    /// return the next Halfedge around the vertex.
    Halfedge_handle vertex_next(Halfedge_handle he);
    
    /// return the opposite halfedge
    Halfedge_handle opposite(Halfedge_handle he);
    
    /// return the signature of a halfedge. True if the edge is flipped. This value is calculated with the sign of the darts of the edge.
    bool signature(Halfedge_handle he);
    
    /// return the dart with the signature s of halfedge he
    Dart_handle Dart(Halfedge_handle, bool s);
    
    /// return the halfedge of a dart
    Halfedge_handle halfedge(Dart);
    
    /// return the sign of a dart.  
    bool signature(Dart_handle d) const;
    
    /// return the sign of a dart.  
    bool signature(Dart_const_handle d) const;
    
    /// return the next Halfedge around the face.
    Halfedge_handle face_next(Halfedge_handle he);
    
    /// @}
    /// \name Range access
    /// @{
    
    /// Returns a range of all the halfedges in the surface. 
    Halfedge_range halfedges();
    
    /// Returns a range of all the halfedges of the i-cell containing *he. 
    template<unsigned int i>
    Halfedge_of_cell<i> halfedge_of_cell(Halfedge_handle he);
    
    /// Returns a range of one halfedge of each i-cell in the generalized map. 
    template<unsigned int i>
    One_halfedge_per_cell<i> one_halfedge_per_cell();
    
    /// @}
    /// \name Operation
    /// @{
    
    /// flip the vertex of he
    void flip(Halfedge_handle he);
    
    ///update all signature
    void update_signature();
    
    /// @}
    /// \name Path and embedded graph
    /// @{
    
    ///range of all Arc_occurence of a halfedge from left to right
    Arc_occurence_Order_range Order(Halfedge_handle) ;
    
    /// create a new path
    Path_handle create_path(); 
    
    ///erase a path
    void erase_path(Path_handle);
    
    ///create a embedded graph
    Graph_handle create_graph();
    
    ///erase a embedded graph
    void erase_graph(Graph_handle);
    
    ///@}
}; 
}