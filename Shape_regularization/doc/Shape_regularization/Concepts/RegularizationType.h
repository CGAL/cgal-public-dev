/*!
\ingroup PkgShape_regularizationConcepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Regularization::Shape_regularization` 
to do regularization.

\cgalHasModel 
- `CGAL::Regularization::Angle_regularization_2`, 
- `CGAL::Regularization::Ordinate_regularization_2`

*/
class RegularizationType {

public:

  /*!  
  here will be my text
  */
  typename GeomTraits::FT bound(const std::size_t i) const { 

  }


  /*!  
    here is my text
  */
  typename GeomTraits::FT target_value(const std::size_t i, const std::size_t j) {

  }

  /*!
    Here will be my text.
  */
  void update(
    const std::vector<std::size_t>& indices) {
    
  }
};
