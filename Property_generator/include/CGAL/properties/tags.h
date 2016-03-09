#ifndef CGAL_PROPERTY_TAGS
#define CGAL_PROPERTY_TAGS

namespace CGAL {
//namespace Properties {

/******************************************************************************/
/*!
\ingroup Property_generator_tags
\brief Tag to disable finiteness tests for function objects in the
`Properties` namespace.

Sometimes \cgal has geometric objects which may neighbour, or themselves
represent, infinite objects. This often means that property functions
operating on those objects require special cases to deal with them, requiring
conditionals. When a user can be sure that this is not possible, they may
provide the `No_finite_test_tag` to avoid these comparisons.
*/

struct No_finite_test_tag {};

/******************************************************************************/
/*!
\ingroup Property_generator_tags
\brief Tag to enable finitness tests in property functors.
This is generally used as the default setting.

Sometimes \cgal has geometric objects which may neighbour, or themselves
represent, infinite objects. This often means that property functions
operating on those objects require special cases to deal with them, requiring
conditionals. This is the default behaviour, and is explicitly activated 
by the `Finite_test_tag`.

*/

struct Finite_test_tag {};

/******************************************************************************/

//} // namespace Properties
} // namespace CGAL

#endif