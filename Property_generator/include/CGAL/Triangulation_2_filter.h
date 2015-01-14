/******************************************************************************/
#ifndef TRIANGULATION_FILTER
#define TRIANGULATION_FILTER
/******************************************************************************/

#include <boost/iterator_adaptors.hpp>

namespace CGAL
{

// This class is only visible from within this file.
// Doxygen will not document it unless this namespace is named.
namespace internal
{

/*

 A filter iterator to work on CGAL iterators.

 We cannot use the default boost implementation because we apply functors
 directly to iterators, so we write our own by using boost's iterator
 adaptor class.

*/

template <typename Functor, typename Iterator, typename Primitive>
class my_filter_iterator
    : public boost::iterator_adaptor<
          my_filter_iterator<Functor, Iterator, Primitive>, /**<  Derived     */
          Iterator,                                         /**<  Base        */
          boost::use_default,                               /**<  Traversal   */
          boost::use_default,                               /**<  Reference   */
          boost::use_default                                /**<  Difference  */
          >
{

 private:
  struct enabler
  {
  };

  Iterator end;
  Functor f;

  /**************************************************************************/

  // public:
  //    my_filter_iterator()
  //        : my_filter_iterator::iterator_adaptor_() {}

  /**************************************************************************/
  // Construct the filter iterator.
 public:
  // Default constructor. Only valid if iterator type and functor type
  // are both default constructible.
  explicit my_filter_iterator()
  {
  }

  explicit my_filter_iterator(Iterator begin, Iterator end, Functor f)
      : my_filter_iterator::iterator_adaptor_(begin), end(end), f(f)
  {
  }

  // Only to be used for creating the 'end' iterator.
  explicit my_filter_iterator(Iterator end, Functor f)
      : my_filter_iterator::iterator_adaptor_(end), end(end), f(f)
  {
  }

  /**************************************************************************/
  // Since the compiler won't do auto conversions for us
  // for more than one level, we have to provide the "Primitive"
  // type information (e.g. Face_handle). This type isn't
  // available in the iterator (I can't find it). So we provide it
  // as a template.

  operator Primitive()
  {
    return this->base();
  }

  /**************************************************************************/

 private:
  friend class boost::iterator_core_access;

  /**************************************************************************/

  // This is called when the iterator is incremented. We skip elements
  // that do not conform, stopping if we hit end.
  void increment()
  {
    Iterator i = this->base_reference();
    ++i;
    // The filter.
    while (i != end && (!f(i)))
      ++i;
    this->base_reference() = i;
  }

  /**************************************************************************/

};  // my_filter_iterator

/******************************************************************************/

/// \addtogroup Property_generator_utilities
/// @{

/**
 * A class to filter 2D triangulations to bounded regions.
 *
 * The Triangulation_2_filter class is provided in order to iterate over
 * 2D triangulations whilst remaining within some bounded region. This is
 * acheived by iterating over all faces, edges or vertices and only visiting
 * those objects for which `on_bounded_side` evaluates to true. Thus, iterating
 * over the object will do one test for every object within the container,
 * regardless if is included within the output range.
 *
 * The class is templated on the region type, so that any region defining
 * `on_bounded_side` taking a `Point_2` may be used.
 */
template <typename Triangulation_2, typename Region>
class Triangulation_2_filter
{

 private:
  // Make 'my_filter_iterator' a friend, so that it can access
  // the special filter operators that we do not want the user to see.
  template <typename U, typename V, typename W>
  friend class my_filter_iterator;

  typedef typename Triangulation_2::Face Face;
  typedef typename Triangulation_2::Vertex Vertex;
  typedef typename Triangulation_2::Geom_traits::Point_2 Point_2;
  typedef typename Triangulation_2::Edge Edge;
  typedef typename Triangulation_2::Face_handle Face_handle;
  typedef typename Triangulation_2::Vertex_handle Vertex_handle;

 public:
/// Iterator over filtered faces.
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Faces_iterator;
#else
  typedef my_filter_iterator<Triangulation_2_filter<Triangulation_2, Region>,
                             typename Triangulation_2::Finite_faces_iterator,
                             Face_handle> Faces_iterator;
#endif

/// Iterator over filtered edges.
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Edges_iterator;
#else
  typedef my_filter_iterator<Triangulation_2_filter<Triangulation_2, Region>,
                             typename Triangulation_2::Finite_edges_iterator,
                             Edge> Edges_iterator;
#endif

/// Iterator over filtered vertices.
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Vertices_iterator;
#else
  typedef my_filter_iterator<Triangulation_2_filter<Triangulation_2, Region>,
                             typename Triangulation_2::Finite_vertices_iterator,
                             Vertex_handle> Vertices_iterator;
#endif

 private:
  const Triangulation_2& tr;
  Region region_;

 public:
  /**************************************************************************/

  /**
   * @{
   * \name Construction
   */

  /**************************************************************************/
  /// Construct a filter given a triangulation and a region.
  Triangulation_2_filter(const Triangulation_2& tr, Region region)
      : tr(tr), region_(region)
  {
  }

  /**/

  Triangulation_2_filter(const Triangulation_2& tr) : tr(tr)
  {
  }

  /**************************************************************************/

  /**
   *
   * @}
   *
   * @{
   * \name Iterators
   */

  /**************************************************************************/

  /// Start of the iterator range over the filtered edges.
  Edges_iterator edges_begin()
  {
    return Edges_iterator(
        tr.finite_edges_begin(), tr.finite_edges_end(), *this);
  }

  /**/

  /// End of the iterator range over the filtered edges.
  Edges_iterator edges_end()
  {
    return Edges_iterator(tr.finite_edges_end(), *this);
  }

  /**************************************************************************/

  /// Start of the iterator range over the filtered faces.
  Faces_iterator faces_begin() const
  {
    return Faces_iterator(
        tr.finite_faces_begin(), tr.finite_faces_end(), *this);
  }

  /**/

  /// End of the iterator range over the filtered faces.
  Faces_iterator faces_end() const
  {
    return Faces_iterator(tr.finite_faces_end(), *this);
  }

  /**************************************************************************/

  /// Start of the iterator range over the filtered vertices.
  Vertices_iterator vertices_begin() const
  {
    return Vertices_iterator(
        tr.finite_vertices_begin(), tr.finite_vertices_end(), *this);
  }

  /**/

  /// End of the iterator range over the filtered vertices.
  Vertices_iterator vertices_end() const
  {
    return Vertices_iterator(tr.finite_vertices_end(), *this);
  }

  /**************************************************************************/

  /**
   *
   * @}
   *
   * @{
   * \name Filter Operators
   * The filter operators return true if the given geometric object
   * is contained within the defined region.
   *
   */

  /**************************************************************************/
  bool operator()(const Face_handle& f) const
  {
    // This seems like the most sensible definition.
    for (int i = 0; i < 3; ++i)
    {
      if (!region_.has_on_bounded_side(f->vertex(i)->point()))
        return false;
    }

    return true;
  }

  /**************************************************************************/

  bool operator()(const Vertex_handle& v) const
  {
    return region_.has_on_bounded_side(v->point());
  }

  /**************************************************************************/

  bool operator()(const Edge& e) const
  {

    const Point_2& p1 = e.first->vertex((e.second + 1) % 3)->point();
    const Point_2& p2 = e.first->vertex((e.second + 2) % 3)->point();

    if (!region_.has_on_bounded_side(p1))
      return false;
    if (!region_.has_on_bounded_side(p2))
      return false;

    return true;
  }

  /**/

 private:
  // Deal with inconsistancy of edge iterator.
  bool operator()(const typename Triangulation_2::Finite_edges_iterator& e)
      const
  {
    return operator()(*e);
  }

 public:
  /**************************************************************************/

  /**
   * @}
   *
   * @{
   *
   * \name Miscellaneous
   *
   */

  /**************************************************************************/

  void set_region(Region r)
  {
    region_ = r;
  }

  Region region() const
  {
    return region_;
  }

  /**************************************************************************/
};

}  // namespace internal

/// @}

}  // namespace CGAL

/******************************************************************************/
#endif
/******************************************************************************/