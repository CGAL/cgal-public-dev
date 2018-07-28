/*!
\ingroup PkgBarycentric_coordinates_2Concepts
\cgalConcept

Requirements of the template parameter `Prior` for Maximum entropy coordinates class from the namespace `CGAL::Barycentric_coordinates`.

  \cgalHasModel `CGAL::Barycentric_coordinates::Maximum_entropy_prior_function_type_one`
  \caglHasModel `CGAL::Barycentric_coordinates::Maximum_entropy_prior_function_type_two`

*/

class MaximumEntropyPrior {

public:

/// \name Creation
/// @{

/// Creates a class that implements maximum entropy prior functions for any query point inside arbitrary simple polygons.
/// The polygon is given by a range of vertices of the type `Traits::Point_2` stored in a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a>.
Maximum_entropy_prior_function_type_one(const std::vector<typename Traits::Point_2> &vertices, const Traits &barycentric_traits);

/// @}

/// \name Creation
/// @{

/// Creates a class that implements maximum entropy prior functions for any query point inside arbitrary simple polygons.
/// The polygon is given by a range of vertices of the type `Traits::Point_2` stored in a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a>.
Maximum_entropy_prior_function_type_two(const std::vector<typename Traits::Point_2> &vertices, const Traits &barycentric_traits);

/// @}

/// \name Functions
/// @{

/*!
	A function that computes maximum entropy prior functions for any query point inside arbitrary simple polygons.
	Prior functions are computed with respect to a query point of the type `Traits::Point_2` and stored in the output iterator `Range`.
*/
template<class Range>
    void compute(typename Traits::Point_2 query_point, Range &m);

/// @}

} /* end MaximumEntropyPrior */
