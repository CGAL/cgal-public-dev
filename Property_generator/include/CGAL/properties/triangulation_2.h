/******************************************************************************/
#ifndef TRIANGULATION_2_PROPERTIES
#define TRIANGULATION_2_PROPERTIES
/******************************************************************************/
// This is the main documentation grouping for Triangulation_2 properties.

namespace CGAL
{
namespace Properties
{
namespace Triangulation_2
{

/*!

\ingroup PkgProperty_generator_properties

\defgroup Triangulation_2_properties Triangulation_2
@{

\brief A namespace to group together functors computing properties of
triangulations.

`#include <CGAL/properties/triangulation_2.h>`

The function objects in this namespace are provided to compute common properties
of 2D triangulations. In general, they are intended for statistics purposes, and
so the values computed are computed using double approximations when required.

By default, the implementations have been written to deal sensibly with input
that may contain infinite vertices. Thus, by default they will perform infinite
tests which may be unnecesary if the user can guarantee that the input is
finite. Tags are provided to disable these tests if required. 

@}

*/

}  // namespace Triangulation_2
}  // namespace Properties
}  // namespace CGAL

/******************************************************************************/
// Include properties for each part of triangulation seperately.

#include "./triangulation_2_edges.h"
#include "./triangulation_2_vertices.h"
#include "./triangulation_2_faces.h"

/******************************************************************************/
#endif
/******************************************************************************/
