/*!
\ingroup PkgBarycentric_coordinates_2Concepts
\cgalConcept

Requirements of the template parameter `Solver` for Maximum entropy coordinates class from the namespace `CGAL::Barycentric_coordinates`.

  \cgalHasModel `CGAL::Barycentric_coordinates::Maximum_entropy_newton_solver`

*/

class MaximumEntropySolver {

public:

/// \name Creation
/// @{

/// Creates a class that implements maximum entropy numerical newton solver.
/// The polygon is given by a range of vertices of the type `Traits::Point_2` stored in a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a>.
Maximum_entropy_newton_solver(const std::vector<typename Traits::Point_2> &vertices, const Traits &barycentric_traits);

/// @}


/// \name Functions
/// @{

/*!
	A function that computes maximum entropy coordinates iteratively based on a numerical newton solver.
	`Type_of_algorithm` determines different iteratation steps and tolerance. The computation results are stored in lambda vector.
*/
void solve(FT_vector &lambda, const Matrix &vtilde, const FT_vector &m, const Type_of_algorithm type_of_algorithm);

/// @}

} /* end MaximumEntropySolver */
