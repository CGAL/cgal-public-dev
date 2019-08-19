/*!
\ingroup PkgShape_regularizationConcepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Regularization::Shape_regularization` 
to access the neighbors of an item.

\cgalHasModel 
- `CGAL::Regularization::Delaunay_neighbor_query_2`
*/
class NeighborQuery {

public:

  /*!  
    fills `neighbors` with the indices of all items, which are connected to the 
    item with the index `query_index`.

    `CGAL::Regularization::Shape_regularization` calls this function each time when 
    a new query item is selected.
  */
  void operator()(
    const std::size_t query_index, 
    std::vector<std::size_t>& neighbors) {
        
  }
};
