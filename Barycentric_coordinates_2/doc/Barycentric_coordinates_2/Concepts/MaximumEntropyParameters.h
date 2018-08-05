/*!
\ingroup PkgBarycentric_coordinates_2Concepts
\cgalConcept

Requirements of the template parameter `Parameters` for Maximum entropy coordinates class from the namespace `CGAL::Barycentric_coordinates`.

  \cgalHasModel `CGAL::Barycentric_coordinates::Maximum_entropy_parameters`

*/

class MaximumEntropyParameters {

public:

/// \name Creation
/// @{

/// Creates a class that posseses max number of iterations and tolerance of Maximum entropy coordinates.
Maximum_entropy_parameters(const Traits &barycentric_traits);

/// @}

/// \name Creation
/// @{

/// Creates a class that posseses max number of iterations and tolerance of Maximum entropy coordinates.
Maximum_entropy_parameters(const size_t max_iter_num, const FT tolerance, const Traits &barycentric_traits);

/// @}


/// \name Functions
/// @{

/*!
	A function that returns the max number of iterations for Maximum entropy coordinates solver.
	`Type_of_algorithm` determines different iteratation steps and tolerance.
*/
size_t max_number_of_iterations(const Type_of_algorithm type_of_algorithm);

FT tolerance(const Type_of_algorithm type_of_algorithm)；

/// @}

/// \name Functions
/// @{

/*!
	A function that returns the tolerance for Maximum entropy coordinates solver.
	`Type_of_algorithm` determines different iteratation steps and tolerance.
*/
FT tolerance(const Type_of_algorithm type_of_algorithm)；

/// @}

} /* end MaximumEntropyParameters */
