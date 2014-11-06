/******************************************************************************/
#ifndef TRIANGULATION_2_FACE_PROPERTIES
#define TRIANGULATION_2_FACE_PROPERTIES
/******************************************************************************/

#include <cmath>
#include <limits>
#include <boost/math/constants/constants.hpp>
#include <CGAL/properties/triangulation_2_edges.h>  // For edge lengths.
#include <CGAL/properties/tags.h>

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

\brief Function objects operating on `Triangulation_2` faces.
\defgroup Face_properties Faces
@{

*/

/******************************************************************************/

// Function objects

/******************************************************************************/

/*!
  Function object to compute the area of a face in a `Triangulation_2`.

  By default the function object checks for infinite vertices, but this checking
  can be disabled by supplying the tag `CGAL::Properties::No_finite_test_tag`.
  In this case, the constructor no longer takes any arguments.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Area
{
 public:
  /// Constructor.
  Area(Triangulation_2 const&);

  /*!
    Operator to compute the area.

    \pre The Vertex_handle provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  double operator()(typename Triangulation_2::Face_handle) const;
};

/******************************************************************************/

/*!
  Function object to compute the aspect ratio of a Triangulation_2 face.

  We define the aspect ratio to be longest edge length of the face divided by
  the length of the shortest edge length of the face.

  By default the function object checks for infinite vertices, but this checking
  can be removed by supplying the tag `CGAL::Properties::No_finite_test_tag`.
  In this case, the constructor no longer takes any arguments.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/
template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Aspect_ratio
{
 public:
  /// Constructor.
  Aspect_ratio(Triangulation_2 const&);

  /*!
    Operator to compute the aspect ratio.

    \pre The Vertex_handle provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  double operator()(typename Triangulation_2::Face_handle) const;
};

/******************************************************************************/

/*!
  Function object to compute the circumradius of a Triangulation_2 face.

  We define the circumradius of a face to be the circumradius of the
  associated circumcircle of the triangle defined by the face.

  By default the function object checks for finite vertices, but this checking
  can be removed by supplying the tag `CGAL::Properties::No_finite_test_tag`.
  In this case, the constructor no longer takes any arguments.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/
template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Circumradius
{
 public:
  /// Constructor.
  Circumradius(Triangulation_2 const&);

  /*!
    Operator to compute the circumradius.

    \pre The Vertex_handle provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  double operator()(typename Triangulation_2::Face_handle) const;
};

/******************************************************************************/

/*!
  Function object to compute the angle within a  Triangulation_2 face.

  By default the function object checks for finite vertices, but this checking
  can be removed by supplying the tag `CGAL::Properties::No_finite_test_tag`.
  In this case, the constructor no longer takes any arguments.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Angle
{
 public:
  /// Constructor.
  Angle(Triangulation_2 const&);

  /*!
    Operator to compute the angle.

    \pre The Vertex_handle provided to the operator must be associated with
    the `Triangulation_2` provided on construction.

    \pre The index provided must be either 0, 1 or 2.
  */
  double operator()(typename Triangulation_2::Face_handle, unsigned) const;
};

/******************************************************************************/

/*!
  Function object to compute the minimum angle of a Triangulation_2 face.

  We return the minimum internal angle of the input face, or zero
  if one of the vertices is at infinity.

  By default the function object checks for finite vertices, but this checking
  can be removed by supplying the Tag `CGAL::Properties::No_finite_test_tag`.
  In this case, the constructor no longer takes any arguments.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Min_angle
{
 public:
  /// Constructor.
  Min_angle(Triangulation_2 const&);

  /*!
    Operator to compute the minimum angle.

    \pre The Vertex_handle provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  double operator()(typename Triangulation_2::Face_handle) const;
};

/******************************************************************************/

/*!
  Function object to compute the max angle of a Triangulation_2 face.

  Return the maximum internal angle of a face, or \f$\pi\f$ if one
  of the vertices is infinite.

  By default the function object checks for finite vertices, but this checking
  can be removed by supplying the tag `CGAL::Properties::No_finite_test_tag`.
  In this case, the constructor no longer takes any arguments.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Max_angle
{
 public:
  /// Constructor.
  Max_angle(Triangulation_2 const&);

  /*!
    Operator to compute the maximum angle.

    \pre The Vertex_handle provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  double operator()(typename Triangulation_2::Face_handle) const;
};

/******************************************************************************/

/*!
  Function object to compute the min edge length of a Triangulation_2 face.

  We define the minimum edge length as the length of the shortest edge
  contained within the input face.

  By default the function object checks for finite vertices, but this checking
  can be removed by supplying the tag `CGAL::Properties::No_finite_test_tag`.
  In this case, the constructor no longer takes any arguments.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/
template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Min_edge_length
{
 public:
  /// Constructor.
  Min_edge_length(Triangulation_2 const&);

  /*!
    Operator to compute the minimum edge length.

    \pre The Vertex_handle provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  double operator()(typename Triangulation_2::Face_handle) const;
};

/******************************************************************************/

/*!
  Function object to compute the max edge length of a Triangulation_2 face.

  We define the maximum edge length to be the length of the longest
  edge contained within this face, or infinity if one of the vertices
  on this face is infinite.

  By default the function object checks for finite vertices, but this checking
  can be removed by supplying the tag `CGAL::Properties::No_finite_test_tag`.
  In this case, the constructor no longer takes any arguments.

  @tparam Triangulation_2 The Triangulation type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag = Finite_test_tag>
class Max_edge_length
{
 public:
  /// Constructor.
  Max_edge_length(Triangulation_2 const&);

  /*!
    Operator to compute the maximum edge length.

    \pre The Vertex_handle provided to the operator must be associated with
    the `Triangulation_2` provided on construction.
  */
  double operator()(typename Triangulation_2::Face_handle) const;
};

/******************************************************************************/

// Free functions

/******************************************************************************/

/*!
  Construct a function object to compute the XYZ for Triangulation_2 Face
  handles.

  @tparam Triangulation_2 The Triangulation_2 type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag>
Area<Triangulation_2, Tag> make_area(const Triangulation_2&, Tag);

template <typename Triangulation_2>
Area<Triangulation_2, No_finite_test_tag> make_area(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a function object to compute the circumradius for Triangulation_2
  Face handles.

  @tparam Triangulation_2 The Triangulation_2 type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag>
Circumradius<Triangulation_2, Tag> make_circumradius(
    const Triangulation_2&,
    Tag);

template <typename Triangulation_2>
Circumradius<Triangulation_2, Finite_test_tag> make_circumradius(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a function object to compute the aspect_ratio for Triangulation_2
  Face handles.

  @tparam Triangulation_2 The Triangulation_2 type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag>
Aspect_ratio<Triangulation_2, Tag> make_aspect_ratio(
    const Triangulation_2&,
    Tag);

template <typename Triangulation_2, typename Tag>
Aspect_ratio<Triangulation_2, Finite_test_tag> make_aspect_ratio(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a function object to compute the min_edge_length for Triangulation_2
  Face handles.

  @tparam Triangulation_2 The Triangulation_2 type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag>
Min_edge_length<Triangulation_2, Tag> make_min_edge_length(
    const Triangulation_2&,
    Tag);

template <typename Triangulation_2>
Min_edge_length<Triangulation_2, Finite_test_tag> make_min_edge_length(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a function object to compute the max_edge_length for Triangulation_2
  Face handles.

  @tparam Triangulation_2 The Triangulation_2 type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag>
Max_edge_length<Triangulation_2, Tag> make_max_edge_length(
    const Triangulation_2&,
    Tag);

template <typename Triangulation_2>
Max_edge_length<Triangulation_2, Finite_test_tag> make_max_edge_length(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a function object to compute the min_angle for Triangulation_2 Face
  handles.

  @tparam Triangulation_2 The Triangulation_2 type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag>
Min_angle<Triangulation_2, Tag> make_min_angle(const Triangulation_2&,
                                                       Tag);

template <typename Triangulation_2>
Min_angle<Triangulation_2, Finite_test_tag> make_min_angle(
    const Triangulation_2&);

/******************************************************************************/

/*!
  Construct a function object to compute the max_angle for Triangulation_2 Face
  handles.

  @tparam Triangulation_2 The Triangulation_2 type.
  @tparam Tag either `CGAL::Properties::Finite_test_tag` or `CGAL::Properties::No_finite_test_tag`.
*/

template <typename Triangulation_2, typename Tag>
Max_angle<Triangulation_2, Tag> make_max_angle(const Triangulation_2&,
                                                       Tag);

template <typename Triangulation_2>
Max_angle<Triangulation_2, Finite_test_tag> make_max_angle(
    const Triangulation_2&);

/******************************************************************************/
// End of documentation and declarations.

/*!

@}

*/

/******************************************************************************/
// Implementations                                                            //
/******************************************************************************/

// forward-declare to break circular dependences.

template <typename T, typename Tag>
class Length;

/******************************************************************************/

template <typename Triangulation_2>
class Area<Triangulation_2, Finite_test_tag>
{
  const Triangulation_2& tr;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  Area(const Triangulation_2& tr) : tr(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Face_handle f) const
  {
    if (tr.is_infinite(f))
      return std::numeric_limits<double>::infinity();

    // Note, code duplicated below.
    return CGAL::area(
        f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
  }
};

/******************************************************************************/

// Specialisation to disable finiteness tests.
template <typename Triangulation_2>
class Area<Triangulation_2, No_finite_test_tag>
{
 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  double operator()(typename Triangulation_2::Face_handle f) const
  {
    // Note, code duplicated above.
    return CGAL::area(
        f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
  }
};

/******************************************************************************/

// Limit this function to use within this unit alone.
namespace
{
// Internal funtion to compute circumradius.
template <typename Triangulation_2>
static inline double circumradius_internal(
    const typename Triangulation_2::Face_handle& f)
{
  typedef typename Triangulation_2::Geom_traits::Point_2 Point;
  const Point& p0 = f->vertex(0)->point();
  const Point& p1 = f->vertex(1)->point();
  const Point& p2 = f->vertex(2)->point();
  // The circumradius is the distance from the centre to any vertex:
  double r2 = (circumcenter(p0, p1, p2) - p0).squared_length();
  return std::sqrt(r2);
}
}  // internal namespace

/******************************************************************************/

template <typename Triangulation_2>
class Circumradius<Triangulation_2, Finite_test_tag>
{
  const Triangulation_2& tr;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  Circumradius(const Triangulation_2& tr) : tr(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Face_handle f) const
  {
    // Infinite faces have infinite circumradius.
    if (tr.is_infinite(f))
      return std::numeric_limits<double>::infinity();
    return circumradius_internal<Triangulation_2>(f);
  }
};

/******************************************************************************/

// Specialisation to disable finiteness tests.
template <typename Triangulation_2>
class Circumradius<Triangulation_2, No_finite_test_tag>
{
 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  double operator()(typename Triangulation_2::Face_handle f) const
  {
    return circumradius_internal<Triangulation_2>(f);
  }
};

/******************************************************************************/

// Limit this function to use within this unit alone.
namespace
{
// Internal function to compute aspect ratio.
template <typename Triangulation_2>
static inline double aspect_ratio_internal(
    const typename Triangulation_2::Face_handle& f)
{
  typedef typename Triangulation_2::Geom_traits::Point_2 Point_2;
  double min = std::numeric_limits<double>::infinity();
  double max = 0.;
  for (unsigned short i = 0; i < 3; ++i)
  {
    const Point_2& p1 = f->vertex(i)->point();
    const Point_2& p2 = f->vertex((i + 1) % 3)->point();

    double value = (p1 - p2).squared_length();
    if (value < min)
      min = value;
    if (value > max)
      max = value;
  }
  return std::sqrt(max) / std::sqrt(min);
}
}

/******************************************************************************/

template <typename Triangulation_2>
class Aspect_ratio<Triangulation_2, Finite_test_tag>
{
  // Allow us to keep a reference to the triangulation.
  const Triangulation_2& tr;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  Aspect_ratio(const Triangulation_2& tr) : tr(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Face_handle f) const
  {
    // Infinite faces have infinite aspect ratio.
    if (tr.is_infinite(f))
      return std::numeric_limits<double>::infinity();

    return aspect_ratio_internal<Triangulation_2>(f);
  }
};

/******************************************************************************/

template <typename Triangulation_2>
class Aspect_ratio<Triangulation_2, No_finite_test_tag>
{
 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  double operator()(typename Triangulation_2::Face_handle f) const
  {
    return aspect_ratio_internal<Triangulation_2>(f);
  }
};

/******************************************************************************/

// Limit this function to use within this unit alone.
namespace
{
// Internal funtion to compute angle.
template <typename Triangulation_2>
static inline double angle_internal(
    const typename Triangulation_2::Face_handle& f,
    unsigned i)
{
  typedef typename Triangulation_2::Geom_traits::Point_2 Point;
  typedef typename Triangulation_2::Geom_traits::Vector_2 Vector;

  Point p1 = f->vertex(i)->point();
  Point p2 = f->vertex((i + 1) % 3)->point();
  Point p3 = f->vertex((i + 2) % 3)->point();

  Vector v1 = (p2 - p1);
  Vector v2 = (p3 - p1);

  v1 = v1 / std::sqrt(v1.squared_length());
  v2 = v2 / std::sqrt(v2.squared_length());

  return std::acos((v1 * v2));
}
}

/******************************************************************************/

template <typename Triangulation_2>
class Angle<Triangulation_2, Finite_test_tag>
{
  const Triangulation_2& tr;

  typedef typename Triangulation_2::Face_handle Face_handle_;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  Angle(const Triangulation_2& tr) : tr(tr)
  {
  }

  result_type operator()(Face_handle_ f, unsigned i) const
  {
    // There are only three angles in a triangle.
    CGAL_precondition(0 <= i && i < 3);

    // Face has an infinite vertex.
    if (tr.is_infinite(f))
    {
      if (tr.is_infinite(f->vertex(i)))
        return 0;
      else
        return boost::math::constants::pi<double>();
    }
    return angle_internal<Triangulation_2>(f, i);
  }
};

/******************************************************************************/
// Specialisation to disable finiteness tests.

template <typename Triangulation_2>
class Angle<Triangulation_2, No_finite_test_tag>
{
 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  double operator()(typename Triangulation_2::Face_handle f, unsigned i) const
  {
    return angle_internal<Triangulation_2>(f, i);
  }
};

/******************************************************************************/

// Limit this function to use within this unit alone.
namespace
{
// Internal funtion to compute min_angle.
template <typename Triangulation_2, typename finite_test_tag>
static inline double min_angle_internal(
    const typename Triangulation_2::Face_handle& f,
    const angle<Triangulation_2, finite_test_tag>& angle_functor)
{
  double min = std::numeric_limits<double>::infinity();
  for (unsigned short i = 0; i < 3; ++i)
  {
    double value = angle_functor(f, i);
    if (value < min)
      min = value;
  }

  return min;
}
}

/******************************************************************************/

template <typename Triangulation_2>
class Min_angle<Triangulation_2, Finite_test_tag>
{
  const Triangulation_2& tr;

  Angle<Triangulation_2> angle_functor;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  min_angle(const Triangulation_2& tr) : tr(tr), angle_functor(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Face_handle f) const
  {
    double min = std::numeric_limits<double>::infinity();
    for (unsigned short i = 0; i < 3; ++i)
    {
      double value = angle_functor(f, i);
      if (value < min)
        min = value;
    }
    return min;
  }
};

/******************************************************************************/
// Specialisation to disable finiteness tests.

template <typename Triangulation_2>
class Min_angle<Triangulation_2, No_finite_test_tag>
{
  // We cache the angle function object locally to save wasted constructions.
  Angle<Triangulation_2, No_finite_test_tag> angle_functor;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  double operator()(typename Triangulation_2::Face_handle f) const
  {
    return min_angle_internal<Triangulation_2, No_finite_test_tag>(
        f, angle_functor);
  }
};

/******************************************************************************/

// Limit this function to use within this unit alone.
namespace
{
// Internal funtion to compute max_angle.
template <typename Triangulation_2, typename finite_test_tag>
static inline double max_angle_internal(
    const typename Triangulation_2::Face_handle& f,
    const Angle<Triangulation_2, finite_test_tag>& angle_functor)
{
  double max = 0.;
  for (unsigned short i = 0; i < 3; ++i)
  {
    double value = angle_functor(f, i);
    if (value > max)
      max = value;
  }
  return max;
}
}

/******************************************************************************/

template <typename Triangulation_2>
class Max_angle<Triangulation_2, Finite_test_tag>
{
  const Triangulation_2& tr;

  Angle<Triangulation_2, Finite_test_tag> angle_functor;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  Max_angle(const Triangulation_2& tr) : tr(tr), angle_functor(tr)
  {
  }

  double operator()(typename Triangulation_2::Face_handle f) const
  {
    return max_angle_internal<Triangulation_2,
                              CGAL::Properties::Finite_test_tag>(f,
                                                                 angle_functor);
  }
};

/******************************************************************************/

// Specialisation to disable finiteness tests.
template <typename Triangulation_2>
class Max_angle<Triangulation_2, No_finite_test_tag>
{
  // We cache the angle function object locally to save wasted constructions.
  Angle<Triangulation_2, No_finite_test_tag> angle_functor;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  double operator()(typename Triangulation_2::Face_handle f) const
  {
    return max_angle_internal<Triangulation_2, No_finite_test_tag>(
        f, angle_functor);
  }
};

/******************************************************************************/

// Limit this function to use within this unit alone.
namespace
{

// Internal funtion to compute min_edge_length.
template <typename Triangulation_2, typename finite_test_tag>
static inline double min_edge_length_internal(
    const typename Triangulation_2::Face_handle& f,
    const length<Triangulation_2, finite_test_tag>& length_functor)
{
  double min = std::numeric_limits<double>::infinity();

  for (unsigned short i = 0; i < 3; ++i)
  {
    double value = length_functor(f, i);
    if (value < min)
      min = value;
  }
  return min;
}
}

/******************************************************************************/

template <typename Triangulation_2>
class Min_edge_length<Triangulation_2, Finite_test_tag>
{
  const Triangulation_2& tr;

  Length<Triangulation_2, Finite_test_tag> length_functor;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  min_edge_length(const Triangulation_2& tr) : tr(tr), length_functor(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Face_handle f) const
  {
    return min_edge_length_internal<Triangulation_2, Finite_test_tag>(
        f, length_functor);
  }
};

/******************************************************************************/
// Specialisation to disable finiteness tests.

template <typename Triangulation_2>
class Min_edge_length<Triangulation_2, No_finite_test_tag>
{
  Length<Triangulation_2, No_finite_test_tag> length_functor;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  double operator()(typename Triangulation_2::Face_handle f) const
  {
    return min_edge_length_internal<Triangulation_2, No_finite_test_tag>(
        f, length_functor);
  }
};

/******************************************************************************/

// Limit this function to use within this unit alone.
namespace
{
// Internal funtion to compute max_edge_length.
template <typename Triangulation_2, typename finite_test_tag>
static inline double max_edge_length_internal(
    const typename Triangulation_2::Face_handle& f,
    const Length<Triangulation_2, finite_test_tag>& length_functor)
{
  double max = 0;
  for (int i = 0; i < 3; ++i)
  {
    double value = length_functor(f, i);
    if (value < max)
      max = value;
  }
  return max;
}
}

/******************************************************************************/

template <typename Triangulation_2>
class Max_edge_length<Triangulation_2, Finite_test_tag>
{
  const Triangulation_2& tr;

  Length<Triangulation_2, Finite_test_tag> length_functor;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  max_edge_length(const Triangulation_2& tr) : tr(tr), length_functor(tr)
  {
  }

  result_type operator()(typename Triangulation_2::Face_handle f) const
  {
    return max_edge_length_internal<Triangulation_2, Finite_test_tag>(
        f, length_functor);
  }
};

/******************************************************************************/
// Specialisation to disable finiteness tests.

template <typename Triangulation_2>
class Max_edge_length<Triangulation_2, No_finite_test_tag>
{
  Length<Triangulation_2, No_finite_test_tag> length_functor;

 public:
  typedef double result_type;
  typedef typename Triangulation_2::Face_handle argument_type;

  double operator()(typename Triangulation_2::Face_handle f) const
  {
    // return aspect_ratio_internal(f);
    return max_edge_length_internal<Triangulation_2, No_finite_test_tag>(
        f, length_functor);
  }
};

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Area<Triangulation_2, finite_test_tag> make_area(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Area<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Area<Triangulation_2, CGAL::Properties::No_finite_test_tag> make_area(
    const Triangulation_2& tr_2)
{
  return Area<Triangulation_2, CGAL::Properties::No_finite_test_tag>();
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Circumradius<Triangulation_2, finite_test_tag> make_circumradius(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Circumradius<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Circumradius<Triangulation_2, CGAL::Properties::Finite_test_tag>
    make_circumradius(const Triangulation_2& tr_2)
{
  return Circumradius<Triangulation_2, CGAL::Properties::Finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Aspect_ratio<Triangulation_2, finite_test_tag> make_aspect_ratio(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Aspect_ratio<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Aspect_ratio<Triangulation_2, CGAL::Properties::Finite_test_tag>
    make_aspect_ratio(const Triangulation_2& tr_2)
{
  return Aspect_ratio<Triangulation_2, CGAL::Properties::Finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Min_edge_length<Triangulation_2, finite_test_tag> make_min_edge_length(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Min_edge_length<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Min_edge_length<Triangulation_2, CGAL::Properties::Finite_test_tag>
    make_min_edge_length(const Triangulation_2& tr_2)
{
  return Min_edge_length<Triangulation_2, CGAL::Properties::Finite_test_tag>(
      tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Max_edge_length<Triangulation_2, finite_test_tag> make_max_edge_length(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Max_edge_length<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Max_edge_length<Triangulation_2, CGAL::Properties::Finite_test_tag>
    make_max_edge_length(const Triangulation_2& tr_2)
{
  return Max_edge_length<Triangulation_2, CGAL::Properties::Finite_test_tag>(
      tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Min_angle<Triangulation_2, finite_test_tag> make_min_angle(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Min_angle<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Min_angle<Triangulation_2, CGAL::Properties::Finite_test_tag>
    make_min_angle(const Triangulation_2& tr_2)
{
  return Min_angle<Triangulation_2, CGAL::Properties::Finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2, typename finite_test_tag>
Max_angle<Triangulation_2, finite_test_tag> make_max_angle(
    const Triangulation_2& tr_2,
    finite_test_tag)
{
  return Max_angle<Triangulation_2, finite_test_tag>(tr_2);
}

/******************************************************************************/

template <typename Triangulation_2>
Max_angle<Triangulation_2, CGAL::Properties::Finite_test_tag>
    make_max_angle(const Triangulation_2& tr_2)
{
  return Max_angle<Triangulation_2, CGAL::Properties::Finite_test_tag>(tr_2);
}

/******************************************************************************/

}  // namespace Triangulation_2
}  // namespace Properties
}  // namespace CGAL

/******************************************************************************/
#endif
/******************************************************************************/
