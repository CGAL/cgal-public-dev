/******************************************************************************/
#ifndef TRIANGULATION_2_VERTEX_PROPERTIES
#define TRIANGULATION_2_VERTEX_PROPERTIES
/******************************************************************************/

#include <cmath>
#include <limits>
#include <CGAL/number_utils.h>
#include <CGAL/properties/tags.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/properties/triangulation_2_faces.h>  // For face area and angles.
#include <CGAL/number_utils.h>

/******************************************************************************/

namespace CGAL
{
namespace Properties
{
namespace Triangulation_2
{

/******************************************************************************/
// This is the main Doxygen group for this file.

/*!

\ingroup Triangulation_2_properties

\brief Function objects operating on `Triangulation_2` vertices.
\defgroup Vertex_properties Vertices
@{

*/

/******************************************************************************/

// Function objects

/******************************************************************************/

/*!
  Function object to compute the degree of a `Triangulation_2` `Vertex_handle`.
  This function object is simply a wrapper provided only for consistency
  purposes.

  @tparam Triangulation_2 The triangulation type.
*/

template <typename Triangulation_2>
class Degree
{
 public:

  /// Result type of the function object.
  typedef unsigned result_type;

  /// Argument type of function object.
  #ifdef DOXYGEN_RUNNING
    typedef Vertex_handle argument_type;
  #else
    typedef typename Triangulation_2::Vertex_handle argument_type;
  #endif

  /// Constructor
  Degree();

  /// Operator to compute degree.
  result_type operator()(argument_type) const;
};

/******************************************************************************/

/*!
  \brief Function object to compute the dual area associated with a
  `Vertex_handle`.

  We define the dual area to be the area of the convex polygon formed by
  connecting up the circumcentres of the triangular faces adjacent to this
  vertex. This area coincides with the area of the Voronoi cell.
  When the dimension is degenerate, we define the dual area to be zero.

  By default this function object checks for infinite faces, returning
  `Geom_traits::FT(std::numeric_limits<double>::infinity())`
  when the dual area  is infinite. Thus the result in this case may result
  in undefined behaviour if `Geom_traits::FT` does not support infinities.
  The infinity checking can be disabled by supplying the tag
  `CGAL::No_finite_test_tag`.

  @tparam Triangulation_2 The Triangulation type.
*/

template <typename Triangulation_2, typename tag = Finite_test_tag>
class Dual_area
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
    typedef Vertex_handle argument_type;
  #else
    typedef typename Triangulation_2::Vertex_handle argument_type;
  #endif

  /// Constructor.
  Dual_area(Triangulation_2 const&);

  /*!
    Operator to compute dual area.

    \pre The `Vertex_handle` provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  result_type operator()(argument_type) const;
};

/******************************************************************************/

/*!
  Function object to compute the star area of a vertex.

  We define the star area to be the sum of the areas of the triangular faces
  that contain this vertex, or zero for degenerate triangulations. 

  By default this function object checks for infinite faces, returning
  `Geom_traits::FT(std::numeric_limits<double>::infinity())`
  when the star area is infinite. Thus the result in this case may result
  in undefined behaviour if `Geom_traits::FT` does not support infinities.
  The infinity checking can be disabled by supplying the tag
  `CGAL::No_finite_test_tag`.

  \pre The `Vertex_handle` provided to the operator must be associated with the
  `Triangulation_2` provided on construction.

  @tparam Triangulation_2 The Triangulation type.
*/

template <typename Triangulation_2, typename tag = Finite_test_tag>
class Star_area
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
    typedef Vertex_handle argument_type;
  #else
    typedef typename Triangulation_2::Vertex_handle argument_type;
  #endif

  /// Constructor.
  Star_area(Triangulation_2 const&);

  /*!
    Operator to compute star area.

    \pre The Vertex_handle provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  result_type operator()(argument_type) const;
};

/******************************************************************************/

/*!
  Function object to compute an approximation to the link length of a
  `Triangulation_2` vertex to double accuracy.

  We define the link length to be the sum of the lengths of the edges not
  incident to the input vertex, but on a triangular face containing it.
  We return zero for degenerate triangulations.

  By default this function object checks for infinite edges, returning
  `Geom_traits::FT(std::numeric_limits<double>::infinity())` when the link
  length is infinite. Thus the result in this case may result in undefined
  behaviour if `Geom_traits::FT` does not support infinities. The infinity
  checking can be disabled by supplying the tag `CGAL::No_finite_test_tag`. If
  needed, `CGAL::to_double()` will be called to use STL functions taking a
  `double`.

  @tparam Triangulation_2 The Triangulation type.
*/

template <typename Triangulation_2, typename tag = Finite_test_tag>
class Link_length
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
    typedef Vertex_handle argument_type;
  #else
    typedef typename Triangulation_2::Vertex_handle argument_type;
  #endif

  /// Constructor.
  Link_length(Triangulation_2 const&);

  /*!
    Operator to compute link length.

    \pre The `Vertex_handle` provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  result_type operator()(argument_type) const;
};

/******************************************************************************/

/*!
  Function object to compute the approximate maximum star angle in radians of a
  vertex to double accuracy.

  We define the maximum star angle to be the largest angle between two adjacent
  edges exiting this vertex, or zero for degenerate triangulations.

  By default this function object checks for infinite faces. The infinity
  checking can be disabled by supplying the tag `CGAL::No_finite_test_tag`. If
  needed, `CGAL::to_double()` will be called to use STL functions taking a
  `double`.

  @tparam Triangulation_2 The Triangulation type.
*/

template <typename Triangulation_2, typename tag = Finite_test_tag>
class Max_star_angle
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
    typedef Vertex_handle argument_type;
  #else
    typedef typename Triangulation_2::Vertex_handle argument_type;
  #endif

  /// Constructor.
  Max_star_angle(Triangulation_2 const&);

  /// Operator to compute maximum star angle in degree.
  result_type operator()(argument_type) const;
};

/******************************************************************************/

/*!
  Function object to compute the approximate minimum star angle in radians of
  this vertex to double accuracy

  We define the minimum star angle to be the smallest angle between two adjacent
  edges exiting this vertex, or zero for degenerate triangulations.

  By default this function object checks for infinite faces. The infinity
  checking can be disabled by supplying the tag `CGAL::No_finite_test_tag`. If
  needed, `CGAL::to_double()` will be called to use STL functions taking a
  `double`.

  @tparam Triangulation_2 The Triangulation type.
*/

template <typename Triangulation_2, typename tag = Finite_test_tag>
class Min_star_angle
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
    typedef Vertex_handle argument_type;
  #else
    typedef typename Triangulation_2::Vertex_handle argument_type;
  #endif

  /// Constructor.
  Min_star_angle(Triangulation_2 const&);

  /*!
    Operator to compute minimum star angle in radians.

    \pre The Vertex_handle provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  result_type operator()(argument_type) const;
};

/******************************************************************************/

// Free functions

/******************************************************************************/

/*!
  Construct a functor to compute the min_star_angle for Triangulation_2
  Vertex handles.  The tag is optional and defaults to `Finite_test_tag`.

  @param  tr_2            The Triangulation_2 to be associated with the
                          functor.
  @tparam Triangulation_2 The Triangulation_2 type.
*/

template <typename Triangulation_2, typename tag>
Min_star_angle<Triangulation_2, tag> make_min_star_angle(
    const Triangulation_2& tr_2,
    tag);

// Default implementation when no tag is provided.
template <typename Triangulation_2>
Min_star_angle<Triangulation_2, Finite_test_tag> make_min_star_angle(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a functor to compute the max_star_angle for Triangulation_2
  Vertex handles.  The tag is optional and defaults to `Finite_test_tag`.

  @param  tr_2            The Triangulation_2 to be associated with the
                         functor.
  @tparam Triangulation_2 The Triangulation_2 type.
*/

template <typename Triangulation_2, typename tag>
Max_star_angle<Triangulation_2, tag> make_max_star_angle(
    const Triangulation_2& tr_2,
    tag);

// Default implementation when no tag is provided.
template <typename Triangulation_2>
Max_star_angle<Triangulation_2, Finite_test_tag> make_max_star_angle(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a functor to compute the link_length for Triangulation_2 Vertex
  handles.  The tag is optional and defaults to `Finite_test_tag`.

  @param  tr_2            The Triangulation_2 to be associated with the
                         functor.
  @tparam Triangulation_2 The Triangulation_2 type.
*/

template <typename Triangulation_2, typename tag>
Link_length<Triangulation_2, tag> make_link_length(const Triangulation_2& tr_2,
                                                   tag);

// Default implementation when no tag is provided.
template <typename Triangulation_2>
Link_length<Triangulation_2, Finite_test_tag> make_link_length(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a functor to compute the star_area for Triangulation_2 Vertex
  handles.  The tag is optional and defaults to `Finite_test_tag`.

  @param  tr_2            The Triangulation_2 to be associated with the
                         functor.
  @tparam Triangulation_2 The Triangulation_2 type.
*/

template <typename Triangulation_2, typename tag>
Star_area<Triangulation_2, tag> make_star_area(const Triangulation_2& tr_2,
                                               tag);

template <typename Triangulation_2>
Star_area<Triangulation_2, Finite_test_tag> make_star_area(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a functor to compute the dual_area for Triangulation_2 Vertex
  handles.  The tag is optional and defaults to `Finite_test_tag`.

  @param  tr_2 The Triangulation_2 to be associated with the functor.
  @tparam Triangulation_2 The Triangulation_2 type.
*/

template <typename Triangulation_2, typename tag>
Dual_area<Triangulation_2, tag> make_dual_area(const Triangulation_2& tr_2,
                                               tag);

template <typename Triangulation_2>
Dual_area<Triangulation_2, Finite_test_tag> make_dual_area(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a functor to compute the degree for Triangulation_2 Vertex
  handles.

  @tparam Triangulation_2 The Triangulation_2 type.
*/

template <typename Triangulation_2>
Degree<Triangulation_2> make_degree(const Triangulation_2&);

/******************************************************************************/
// End of documentation and declarations.

/*!

@}

*/

/******************************************************************************/
// Implementations                                                            //
/******************************************************************************/

// forward-declare to break circular dependences.

template <typename T, typename U>
class area;

template <typename T, typename U>
class angle;

/******************************************************************************/

// No partial specialisation for the degree function object, provide
// implementation for primary template.

template <typename Tr_2>
Degree<Tr_2>::Degree()
{
}

/******************************************************************************/

template <typename Tr_2>
unsigned Degree<Tr_2>::operator()(typename Tr_2::Vertex_handle v) const
{
  return v->degree();
}

/******************************************************************************/

template <typename Triangulation_2>
class Dual_area<Triangulation_2, Finite_test_tag>
{
  typedef typename Triangulation_2::Vertex_handle Vertex_handle;
  typedef typename Triangulation_2::Geom_traits Gt;
  typedef typename Triangulation_2::Face_circulator Face_circulator;

  const Triangulation_2& tr;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;
  typedef typename Triangulation_2::Vertex_handle argument_type;


  Dual_area(const Triangulation_2& tr) : tr(tr)
  {
  }

  // NOTE: Code duplicated below.
  result_type operator()(typename Triangulation_2::Vertex_handle v) const
  {
    // The measure of lower dimensional subspaces is zero
    // under the two dimensional measure.
    if (tr.dimension() < 2)
      return 0;

    // The dual area of the infinite vertex is infinite.
    if (tr.is_infinite(v))
      return std::numeric_limits<double>::infinity();

    Polygon_2<Gt> p;

    Face_circulator f = tr.incident_faces(v), done(f);
    
    // A vertex with no neighbours has infinite dual area. 
    if (f != 0)
    {      
      do
      {      
        if (tr.is_infinite(f))
          return std::numeric_limits<double>::infinity();
        p.push_back(tr.circumcenter(f));
      }
      while (++f != done);
    }

    return p.area();
  }
};

/******************************************************************************/

template <typename Triangulation_2>
class Dual_area<Triangulation_2, No_finite_test_tag>
{
  typedef typename Triangulation_2::Face_circulator Face_circulator;
  const Triangulation_2& tr;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;
  typedef typename Triangulation_2::Vertex_handle argument_type;


  // We need a reference to the triangulation to compute this property.
  Dual_area(const Triangulation_2& tr) : tr(tr)
  {
  }

  // NOTE: Code duplicated above.
  result_type operator()(typename Triangulation_2::Vertex_handle v) const
  {
    // The measure of lower dimensional subspaces is zero
    // under the two dimensional measure.
    if (tr.dimension() < 2)
      return 0;

    Polygon_2<typename Triangulation_2::Geom_traits> p;
    Face_circulator f = tr.incident_faces(v), done(f);
    if (f != 0)
    {
      do
        p.push_back(tr.circumcenter(f));
      while (++f != done);
    }

    return abs(p.area());
  }
};

/******************************************************************************/
// With finite tests.

template <typename Triangulation_2>
class Star_area<Triangulation_2, Finite_test_tag>
{
  typedef typename Triangulation_2::Vertex_handle Vertex_handle;
  typedef typename Triangulation_2::Face_circulator Face_circulator;

  const Triangulation_2& tr;
  Area<Triangulation_2, Finite_test_tag> area_functor;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;
  typedef typename Triangulation_2::Vertex_handle argument_type;

  Star_area(const Triangulation_2& tr) : tr(tr), area_functor(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Vertex_handle v) const
  {

    if (tr.dimension() < 2)
      return 0;

    // NOTE: Code duplicated below.
    result_type area = 0.;
    Face_circulator f = tr.incident_faces(v), done(f);
    if (f != 0)
    {
      do
      {
        area += area_functor(f);
      } while (++f != done);
    }

    return area;
  }
};

/******************************************************************************/
// No finite tests.

template <typename Triangulation_2>
class Star_area<Triangulation_2, No_finite_test_tag>
{
  typedef typename Triangulation_2::Face_circulator Face_circulator;

  const Triangulation_2& tr;
  Area<Triangulation_2, No_finite_test_tag> area_functor;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;
  typedef typename Triangulation_2::Vertex_handle argument_type;  

  Star_area(const Triangulation_2& tr) : tr(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Vertex_handle v) const
  {
    if (tr.dimension() < 2)
      return 0;

    // NOTE: Code duplicated above.
    result_type area = 0;
    Face_circulator f = tr.incident_faces(v), done(f);
    if (f != 0)
    {
      do
      {
        area += area_functor(f);
      } while (++f != done);
    }
    return area;
  }
};

/******************************************************************************/
// With finite tests.

template <typename Triangulation_2>
class Link_length<Triangulation_2, Finite_test_tag>
{
  typedef typename Triangulation_2::Vertex_handle Vertex_handle;
  typedef typename Triangulation_2::Face_circulator Face_circulator;

  const Triangulation_2& tr;
  Length<Triangulation_2> length_functor;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;
  typedef typename Triangulation_2::Vertex_handle argument_type;


  Link_length(const Triangulation_2& tr) : tr(tr), length_functor(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Vertex_handle v) const
  {
    // NOTE: Code duplicated below.

    if (tr.dimension() < 2)
      return 0;
    
    result_type length = 0.;
    Face_circulator f = tr.incident_faces(v), done(f);
    if (f != 0)
    {
      do{
        length += length_functor(f, f->index(v));
      }
      while (++f != done);
    }
    return length;
  }
};

/******************************************************************************/

// Specialisation to disable finiteness tests.
template <typename Triangulation_2>
class Link_length<Triangulation_2, No_finite_test_tag>
{
  typedef typename Triangulation_2::Face_circulator Face_circulator_;
  
  const Triangulation_2& tr;
  Length<Triangulation_2, No_finite_test_tag> Length_functor;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;
  typedef typename Triangulation_2::Vertex_handle argument_type;

  Link_length(const Triangulation_2& tr) : tr(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Vertex_handle v) const
  {

    // NOTE: Code duplicated above.

    if (tr.dimension() < 2)
      return 0;

    result_type length = 0;
    Face_circulator_ f = tr.incident_faces(v), done(f);
    if (f != 0)
    {
      do
        length += Length_functor(f, f->index(v));
      while (++f != done);
    }
    return length;
  }
};

/******************************************************************************/

template <typename Triangulation_2>
class Max_star_angle<Triangulation_2, Finite_test_tag>
{
  typedef typename Triangulation_2::Vertex_handle Vertex_handle_;
  typedef typename Triangulation_2::Face_circulator Face_circulator_;

  const Triangulation_2& tr;

  Angle<Triangulation_2, Finite_test_tag> angle_functor;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;
  typedef typename Triangulation_2::Vertex_handle argument_type;

  Max_star_angle(const Triangulation_2& tr) : tr(tr), angle_functor(tr)
  {
  }

  result_type operator()(Vertex_handle_ v) const
  {
    if (tr.dimension() < 2)
      return 0;

    // NOTE: Code duplicated below.
    result_type max = 0.;
    Face_circulator_ f = tr.incident_faces(v), done(f);
    if (f != 0)
    {
      do
      {
        result_type value = angle_functor(f, f->index(v));
        if (value > max)
          max = value;
      } while (++f != done);
    }

    return max;
  }
};

/******************************************************************************/
// Specialisation to disable finiteness tests.

template <typename Triangulation_2>
class Max_star_angle<Triangulation_2, No_finite_test_tag>
{
  typedef typename Triangulation_2::Face_circulator Face_circulator_;
  const Triangulation_2& tr;
  Angle<Triangulation_2, No_finite_test_tag> angle_functor;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;
  typedef typename Triangulation_2::Vertex_handle argument_type;


  Max_star_angle(const Triangulation_2& tr) : tr(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Vertex_handle v) const
  {
    // NOTE: Code duplicated above.
    if (tr.dimension() < 2)
      return 0;

    result_type max = 0;
    Face_circulator_ f = tr.incident_faces(v), done(f);
    if (f != 0)
    {
      do
      {
        result_type value = angle_functor(f, f->index(v));
        if (value > max)
          max = value;
      } while (++f != done);
    }

    return max;
  }
};

/******************************************************************************/

template <typename Triangulation_2>
class Min_star_angle<Triangulation_2, Finite_test_tag>
{
  typedef typename Triangulation_2::Vertex_handle Vertex_handle_;
  typedef typename Triangulation_2::Face_circulator Face_circulator_;
  const Triangulation_2& tr;
  Angle<Triangulation_2, Finite_test_tag> angle_functor;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;
  typedef typename Triangulation_2::Vertex_handle argument_type;


  Min_star_angle(const Triangulation_2& tr) : tr(tr), angle_functor(tr)
  {
  }

  result_type operator()(Vertex_handle_ v) const
  {
    // NOTE: Code duplicated below.
    if (tr.dimension() < 2)
      return 0;

    result_type min = CGAL_PI;
    Face_circulator_ f = tr.incident_faces(v), done(f);
    if (f != 0)
    {
      do
      {
        result_type value = angle_functor(f, f->index(v));
        if (value < min)
          min = value;
      } while (++f != done);
    }
    return min;
  }
};

/******************************************************************************/
// Specialisation to disable finiteness tests.

template <typename Triangulation_2>
class Min_star_angle<Triangulation_2, No_finite_test_tag>
{
  typedef typename Triangulation_2::Face_circulator Face_circulator;

  const Triangulation_2& tr;
  Angle<Triangulation_2, No_finite_test_tag> angle_functor;

 public:
  typedef typename Triangulation_2::Geom_traits::FT result_type;
  typedef typename Triangulation_2::Vertex_handle argument_type;


  Min_star_angle(const Triangulation_2& tr) : tr(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Vertex_handle v) const
  {
    // NOTE: Code duplicated above.
    if (tr.dimension() < 2)
      return 0;

    result_type min = CGAL_PI;
    Face_circulator f = tr.incident_faces(v), done(f);
    if (f != 0)
    {
      do
      {
        result_type value = angle_functor(f, f->index(v));

        if (value < min)
          min = value;
      } while (++f != done);
    }

    return min;
  }
};

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Min_star_angle<Triangulation_2, finite_test_tag> make_min_star_angle(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Min_star_angle<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Min_star_angle<Triangulation_2, Finite_test_tag> make_min_star_angle(
    const Triangulation_2& tr_2)
{
  return Min_star_angle<Triangulation_2, Finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Max_star_angle<Triangulation_2, finite_test_tag> make_max_star_angle(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Max_star_angle<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Max_star_angle<Triangulation_2, Finite_test_tag> make_max_star_angle(
    const Triangulation_2& tr_2)
{
  return Max_star_angle<Triangulation_2, Finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Link_length<Triangulation_2, finite_test_tag> make_link_length(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Link_length<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Link_length<Triangulation_2, Finite_test_tag> make_link_length(
    const Triangulation_2& tr_2)
{
  return Link_length<Triangulation_2, Finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Star_area<Triangulation_2, finite_test_tag> make_star_area(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Star_area<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Star_area<Triangulation_2, Finite_test_tag> make_star_area(
    const Triangulation_2& tr_2)
{
  return Star_area<Triangulation_2, Finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Dual_area<Triangulation_2, finite_test_tag> make_dual_area(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Dual_area<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Dual_area<Triangulation_2, Finite_test_tag> make_dual_area(
    const Triangulation_2& tr_2)
{
  return Dual_area<Triangulation_2, Finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Degree<Triangulation_2> make_degree(const Triangulation_2& tr_2)
{
  return Degree<Triangulation_2>();
}

/******************************************************************************/

}  // namespace Triangulation_2
}  // namespace Properties
}  // namespace CGAL

/******************************************************************************/
#endif
/******************************************************************************/
