/******************************************************************************/
#ifndef TRIANGULATION_2_EDGE_PROPERTIES
#define TRIANGULATION_2_EDGE_PROPERTIES
/******************************************************************************/

#include <cmath>
#include <limits>
#include <CGAL/properties/tags.h>
#include <CGAL/properties/triangulation_2_faces.h>  // For face area.

/******************************************************************************/

namespace CGAL
{
namespace Properties
{
namespace Triangulation_2
{

/******************************************************************************/
// Doxygen group for file.

/*!

\ingroup Triangulation_2_properties

\brief Function objects operating on `Triangulation_2` edges.
\defgroup Edge_properties Edges
@{

*/

/******************************************************************************/

// Function objects

/******************************************************************************/

/*!
  Function object to compute the approximate length of an edge to double
  accuracy.

  By default the function object checks that the supporting vertices of the edge
  are both finite, returning
  `Geom_traits::FT(std::numeric_limits<double>::infinity())` when the edge is
  infinite. Thus the result in this case may result in undefined behaviour if
  `Geom_traits::FT` does not support infinities. This checking can be disabled
  by supplying the tag `CGAL::No_finite_test_tag`, in which case the constructor
  requires no argument. If needed `CGAL::to_double()` will be called to use STL
  functions taking a `double` as input.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Finite_test_tag` or
  `CGAL::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Length
{
 public:

  /// Result type of function object.
  #ifdef DOXYGEN_RUNNING
    typedef Geom_traits::FT result_type;
  #else
    typedef typename Triangulation_2::Geom_traits::FT result_type;
  #endif

  /// Argument type of function object.
  #ifdef DOXYGEN_RUNNING
    typedef Edge argument_type;
  #else
    typedef typename Triangulation_2::Edge argument_type;
  #endif


  /// Constructor.
  Length(Triangulation_2 const&);

  /*!
    Operator to compute the edge length.

    \pre The `Edge` provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  result_type operator()(argument_type) const;
};

/******************************************************************************/

/*!
  Function object to compute the neighbor area of an edge.
  We define the neighbor area to be the sum of the areas of the two triangular
  faces incident to the given edge.

  By default the function object checks that the neighboring faces of the edge
  are finite, returning
  `Geom_traits::FT(std::numeric_limits<double>::infinity())` when one is
  infinite. Thus the result in this case may result in undefined behaviour if
  `Geom_traits::FT` does not support infinities. This checking can be disabled
  by supplying the tag `CGAL::No_finite_test_tag`. If needed `CGAL::to_double()`
  will be called to use STL functions taking a `double` as input.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Finite_test_tag` or
  `CGAL::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Neighbor_area
{
 public:

  /// Result type of function object.
  #ifdef DOXYGEN_RUNNING
    typedef Geom_traits::FT result_type;
  #else
    typedef typename Triangulation_2::Geom_traits::FT result_type;
  #endif

  /// Argument type of function object.
  #ifdef DOXYGEN_RUNNING
    typedef Edge argument_type;
  #else
    typedef typename Triangulation_2::Edge argument_type;
  #endif

  /// Constructor.
  Neighbor_area(Triangulation_2 const&);

  /*!
    Operator to compute neighbor area.

    \pre The `Edge` provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  result_type operator()(argument_type) const;
};

/******************************************************************************/

/*!
  Function object to compute the approximate dual length of an edge, to double
  accuracy. We define the dual length to be the length of the line segment
  between the circumcenters of the two triangular faces adjacent to the given
  edge, or zero when the triangulation is degenerate.

  \note
  Iterating over finite edges in a Triangulation_2 will always result in
  infinite dual edges.

  By default the function object checks for infinite vertices, returning
  `Geom_traits::FT(std::numeric_limits<double>::infinity())` when the dual
  length is infinite. Thus the result in this case may result in undefined
  behaviour if `Geom_traits::FT` does not support infinities. This checking can
  be disabled by supplying the tag `CGAL::No_finite_test_tag`. If needed
  `CGAL::to_double()` will be called to use STL functions taking a `double` as
  input.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Finite_test_tag` or
  `CGAL::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Dual_length
{
 public:

  /// Result type of function object.
  #ifdef DOXYGEN_RUNNING
    typedef Geom_traits::FT result_type;
  #else
    typedef typename Triangulation_2::Geom_traits::FT result_type;
  #endif

  /// Argument type of function object.
  #ifdef DOXYGEN_RUNNING
    typedef Edge argument_type;
  #else
    typedef typename Triangulation_2::Edge argument_type;
  #endif

  /// Constructor.
  Dual_length(Triangulation_2 const&);

  /*!
    Operator to compute dual length.

    \pre The `Edge` provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  result_type operator()(argument_type) const;
};

/******************************************************************************/

// Free functions

/******************************************************************************/

/*!
  Construct a functor to compute edge lengths. The tag is optional and defaults
  to `Finite_test_tag`.

  @param  tr_2  The Triangulation_2 to be associated with the
          functor.
  @tparam Triangulation_2 The Triangulation_2 type.
  @tparam Tag either `CGAL::Finite_test_tag` or
  `CGAL::No_finite_test_tag`.
*/
template <typename Triangulation_2, typename Tag>
Length<Triangulation_2, Tag> make_length(const Triangulation_2& tr_2, Tag);
template <typename Triangulation_2>
Length<Triangulation_2, Finite_test_tag> make_length(
    const Triangulation_2& tr_2);

/******************************************************************************/

/*!
  Construct a functor to compute the neighbor area of an edge. The tag is
  optional and defaults to `Finite_test_tag`.

  @param  tr_2            The Triangulation_2 to be associated with the
                       functor.
  @tparam Triangulation_2 The Triangulation_2 type.
  @tparam Tag either `CGAL::Finite_test_tag` or
  `CGAL::No_finite_test_tag`.
*/
template <typename Triangulation_2, typename Tag>
Neighbor_area<Triangulation_2, Tag> make_neighbor_area(
    const Triangulation_2& tr_2,
    Tag);

template <typename Triangulation_2>
Neighbor_area<Triangulation_2, Finite_test_tag> make_neighbor_area(
    const Triangulation_2& tr_2);

/******************************************************************************/

/*!
  Construct a functor to compute the dual length of an edge.
  @param  tr_2            The Triangulation_2 to be associated with the
                          functor.
  @tparam Triangulation_2 The Triangulation_2 type.
  @tparam Tag either `CGAL::Finite_test_tag` or
  `CGAL::No_finite_test_tag`.
*/
template <typename Triangulation_2, typename Tag>
Dual_length<Triangulation_2, Tag> make_dual_length(const Triangulation_2& tr_2,
                                                   Tag);

template <typename Triangulation_2>
Dual_length<Triangulation_2, CGAL::Finite_test_tag>
    make_dual_length(const Triangulation_2& tr_2);

/******************************************************************************/
// End of documentation and declaration.

/*!

@}

*/

/******************************************************************************/
// Implementations
/******************************************************************************/

// forward-declare to break circular dependences.

template <typename T, typename Tag>
class Area;

/******************************************************************************/

template <typename Triangulation_2>
class Length<Triangulation_2, Finite_test_tag>
{
  typedef typename Triangulation_2::Vertex_handle Vertex_handle_;
  typedef typename Triangulation_2::Edge Edge_;
  typedef typename Triangulation_2::Face_handle Face_handle_;

  const Triangulation_2& tr;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;

  Length(const Triangulation_2& tr) : tr(tr) {};

  result_type operator()(Edge_ e) const
  {
    return operator()(e.first, e.second);
  }

  result_type operator()(Face_handle_ f, unsigned short i) const
  {
    Vertex_handle_ v1 = f->vertex(f->cw(i));
    Vertex_handle_ v2 = f->vertex(f->ccw(i));

    if ( tr.is_infinite(v1) || tr.is_infinite(v2) )
      return std::numeric_limits<double>::infinity();

    return sqrt(to_double((v1->point() - v2->point()).squared_length()));
  }
};

/******************************************************************************/

// Specialisation to disable finiteness tests.
template <typename Triangulation_2>
class Length<Triangulation_2, No_finite_test_tag>
{
  typedef typename Triangulation_2::Face_handle Face_handle_;
  typedef typename Triangulation_2::Vertex_handle Vertex_handle_;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;

  result_type operator()(typename Triangulation_2::Edge e) const
  {
    return operator()(e.first, e.second);
  }

  result_type operator()(Face_handle_ f, unsigned short i) const
  {

    Vertex_handle_ v1 = f->vertex(f->cw(i));
    Vertex_handle_ v2 = f->vertex(f->ccw(i));

    return sqrt(to_double((v1->point() - v2->point()).squared_length()));
  }
};

/******************************************************************************/

template <typename Triangulation_2>
class Neighbor_area<Triangulation_2, Finite_test_tag>
{
  typedef typename Triangulation_2::Edge Edge_;
  typedef typename Triangulation_2::Face_handle Face_handle_;

  const Triangulation_2& tr;

  // Functor to compute face area.
  Area<Triangulation_2, Finite_test_tag> area_functor;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;

  Neighbor_area(const Triangulation_2& tr) : tr(tr), area_functor(tr)
  {
  }

  result_type operator()(Edge_ e) const
  {
    return operator()(e.first, e.second);
  }

  result_type operator()(Face_handle_ f, int i) const
  {
    // NOTE: code is duplicated below.

    // By convention.
    if (tr.dimension() < 2)
      return 0;

    return area_functor(f) + area_functor(f->neighbor(i));
  }
};

/******************************************************************************/
// Specialisation to disable finiteness tests.

template <typename Triangulation_2>
class Neighbor_area<Triangulation_2, No_finite_test_tag>
{
  typedef typename Triangulation_2::Face_handle Face_handle_;

  const Triangulation_2& tr;

  // Functor to compute face areas.
  Area<Triangulation_2, No_finite_test_tag> area_functor;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;

  Neighbor_area(const Triangulation_2& tr) : tr(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Edge e) const
  {
    return operator()(e.first, e.second);
  }
  result_type operator()(Face_handle_ f, unsigned short i) const
  {
    // NOTE: code is duplicated above.
    
    // By convention.
    if (tr.dimension() < 2)
      return 0;

    return area_functor(f) + area_functor(f->neighbor(i));
  }
};

/******************************************************************************/
// Limit this function to use within this unit alone.

namespace
{
// Internal funtion to compute circumradius.
// We take two faces, as it makes the code slightly more efficient.
template <typename Triangulation_2>
static inline typename Triangulation_2::Geom_traits::FT dual_length_internal(
    const typename Triangulation_2::Face_handle& f1,
    const typename Triangulation_2::Face_handle& f2)
{
  typedef typename Triangulation_2::Geom_traits::Point_2 Point;

  const Point& p1 = f1->vertex(0)->point();
  const Point& p2 = f1->vertex(1)->point();
  const Point& p3 = f1->vertex(2)->point();

  const Point& q1 = f2->vertex(0)->point();
  const Point& q2 = f2->vertex(1)->point();
  const Point& q3 = f2->vertex(2)->point();

  const Point& c1 = circumcenter(p1, p2, p3);
  const Point& c2 = circumcenter(q1, q2, q3);
  return std::sqrt(to_double((c1 - c2).squared_length()));
}
}

/******************************************************************************/

template <typename Triangulation_2>
class Dual_length<Triangulation_2, Finite_test_tag>
{
  typedef typename Triangulation_2::Edge Edge_;
  typedef typename Triangulation_2::Face_handle Face_handle_;

  // Allow us to keep a reference to the triangulation.
  const Triangulation_2& tr;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;

  Dual_length(const Triangulation_2& tr) : tr(tr) {};

  result_type operator()(Edge_ e) const
  {
     if (tr.dimension() < 2)
      return 0;


    Face_handle_ f1 = e.first;
    Face_handle_ f2 = f1->neighbor(e.second);

    if (tr.is_infinite(f1) || tr.is_infinite(f2))
      return std::numeric_limits<double>::infinity();

    return dual_length_internal<Triangulation_2>(f1, f2);
  }
};

/******************************************************************************/
// Specialisation to disable finiteness tests.

template <typename Triangulation_2>
class Dual_length<Triangulation_2, No_finite_test_tag>
{
  typedef typename Triangulation_2::Face_handle Face_handle_;

  // Allow us to keep a reference to the triangulation.
  const Triangulation_2& tr;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;

  Dual_length(const Triangulation_2& tr) : tr(tr) {};

  result_type operator()(typename Triangulation_2::Edge e) const
  {
    if (tr.dimension() < 2)
      return 0;

    Face_handle_ f1 = e.first;
    Face_handle_ f2 = f1->neighbor(e.second);
    return dual_length_internal<Triangulation_2>(f1, f2);
  }
};

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Length<Triangulation_2, finite_test_tag> make_length(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Length<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Length<Triangulation_2, CGAL::No_finite_test_tag> make_length(
    const Triangulation_2& tr_2,
    CGAL::No_finite_test_tag)
{
  return Length<Triangulation_2, CGAL::No_finite_test_tag>();
}

/******************************************************************************/

template <typename Triangulation_2>
Length<Triangulation_2, CGAL::Finite_test_tag> make_length(
    const Triangulation_2& tr_2)
{
  return Length<Triangulation_2, CGAL::Finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Neighbor_area<Triangulation_2, finite_test_tag> make_neighbor_area(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Neighbor_area<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Neighbor_area<Triangulation_2, CGAL::Finite_test_tag>
    make_neighbor_area(const Triangulation_2& tr_2)
{
  return Neighbor_area<Triangulation_2, CGAL::Finite_test_tag>(
      tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Dual_length<Triangulation_2, finite_test_tag> make_dual_length(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Dual_length<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Dual_length<Triangulation_2, CGAL::Finite_test_tag>
    make_dual_length(const Triangulation_2& tr_2)
{
  return Dual_length<Triangulation_2, CGAL::Finite_test_tag>(tr_2);
}

/******************************************************************************/

}  // namespace Triangulation_2
}  // namespace Properties
}  // namespace CGAL

/******************************************************************************/
#endif
/******************************************************************************/