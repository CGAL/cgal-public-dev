
namespace CGAL {

/*!
\ingroup PkgTopologicalInvariantClasses

The class `Generalized_map_min_items` is a model of the `GeneralizedMapItems`
concept. It defines the type of darts which is a
`Dart<d,CMap>`. The `Generalized_map_min_items` has a
template argument for the dimension of the combinatorial map.
In this class, no attribute is enabled.

\cgalModels `GeneralizedMapItems`

\cgalHeading{Example}

The following example shows the implementation of the
`Generalized_map_min_items` class.

\code{.cpp}
template <unsigned int d>
struct Generalized_map_min_items
{
   template < class Refs >
   struct Dart_wrapper
   {
     typedef CGAL::GMap_dart< d, Refs > Dart;
     typedef CGAL::cpp11::tuple<> Attributes;
   };
};
\endcode

\sa `Generalized_map<d,Items,Alloc>`
\sa `Generalized_dart<d,CMap>`

*/
template< typename d >
class Generalized_map_min_items {
public:

/// @}

}; /* end Combinatorial_map_min_items */
} /* end namespace CGAL */
