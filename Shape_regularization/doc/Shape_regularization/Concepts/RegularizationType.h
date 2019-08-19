/*!
\ingroup PkgShape_regularizationConcepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Regularization::Shape_regularization` 
to apply regularization.

\cgalHasModel 
- `CGAL::Regularization::Angle_regularization_2`, 
- `CGAL::Regularization::Ordinate_regularization_2`
*/
class RegularizationType {

public:

  /*!  
  returns the bound of the item with index `i`.

  `CGAL::Regularization::Shape_regularization` calls this function for each item 
  that participates in the regularization process.
  */
  typename GeomTraits::FT bound(const std::size_t i) const { 

  }


  /*!  
    returns the target value between two neighboring items with indices `i` and `j`.

    `CGAL::Regularization::Shape_regularization` calls this function for each pair 
    of neighboring items that participate in the regularization process.
  */
  typename GeomTraits::FT target_value(const std::size_t i, const std::size_t j) {

  }

  /*!
    applies the results from the QP solver to the initial items.

    `CGAL::Regularization::Shape_regularization` calls this function one time, after
    the QP problem has been solved during the regularization process.
  */
  void update(
    const std::vector<GeomTraits::FT> & result) {
    
  }
};
