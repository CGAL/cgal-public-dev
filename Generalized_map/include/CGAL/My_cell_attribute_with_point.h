#ifndef SRC_MY_CELL_ATTRIBUTE_WITH_POINT_H
#define SRC_MY_CELL_ATTRIBUTE_WITH_POINT_H

#include <CGAL/Cell_attribute.h>

namespace CGAL {

  /** @file Cell_attribute_with_point.h
   * Definition of cell attribute with point, with or without info.
   */

  /// Point associated with a cell.
  template < class Point >
  class My_point_for_cell
  {
  public:
    /// Contructor without parameter.
    My_point_for_cell()
    {}

    /// Contructor with a point in parameter.
    My_point_for_cell(const Point& apoint) : mpoint(apoint)
    {}

    /// Get the point associated with the cell.
    Point& point()
    { return mpoint; }

    /// Get the point associated with the cell.
    const Point& point() const
    { return mpoint; }

  protected:
    /// The point associated with the cell.
    Point mpoint;
  };

  /// Attribute associated with a point and an info.
  template < class LCC, class Info_=void, class Tag=Tag_true,
             class Functor_on_merge_=Null_functor,
             class Functor_on_split_=Null_functor >
  class My_cell_attribute_with_point :
    public Cell_attribute<LCC, Info_, Tag,
                          Functor_on_merge_, Functor_on_split_>,
    public Point_for_cell<typename LCC::Point>
  {
    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

  public:
    int index0;
    int prev;
    double distance;
    int used;
    int exist;
    int prevfinal;
    bool boundary;
    typedef My_cell_attribute_with_point<LCC, Info_, Tag, Functor_on_merge_,
                                      Functor_on_split_> Self;

    typedef Cell_attribute<LCC, Info_, Tag,
                           Functor_on_merge_, Functor_on_split_> Base1;
    typedef Point_for_cell<typename LCC::Point> Base2;

    typedef typename LCC::Point             Point;
    typedef typename LCC::Dart_handle       Dart_handle;
    typedef typename LCC::Dart_const_handle Dart_const_handle;

    typedef Info_                Info;
    typedef Functor_on_merge_    Functor_on_merge;
    typedef Functor_on_split_    Functor_on_split;

    using Base1::info;

    bool operator==(const Self& other) const
    { return Base1::operator==(other) && this->point()==other.point(); }

    bool operator!=(const Self& other) const
    { return !operator==(other); }

  protected:
    /// Default contructor.
    My_cell_attribute_with_point()
    {}

    /// Contructor with a point in parameter.
    My_cell_attribute_with_point(const Point& apoint) : Base2(apoint)
    {}

    /// Contructor with a point and an attribute in parameters.
    My_cell_attribute_with_point(const Point& apoint, const Info& ainfo) :
      Base1(ainfo),
      Base2(apoint)
    {}
  };

  /// Attribute associated with a point and without info.
  template < class LCC, class Tag,
             class Functor_on_merge_,
             class Functor_on_split_ >
  class My_cell_attribute_with_point<LCC, void, Tag,
                                  Functor_on_merge_, Functor_on_split_> :
    public Cell_attribute<LCC, void, Tag, Functor_on_merge_, Functor_on_split_>,
    public Point_for_cell<typename LCC::Point>
  {
    template <class, class, class, class>
    friend class Compact_container;

    template <class, class>
    friend class Concurrent_compact_container;

  public:

    int index0;
    int prev;

    double distance;
    typedef Cell_attribute<LCC, void, Tag,
                           Functor_on_merge_, Functor_on_split_> Base1;
    typedef Point_for_cell<typename LCC::Point> Base2;

    typedef void                            Info;
    typedef typename LCC::Point             Point;
    typedef typename LCC::Dart_handle       Dart_handle;
    typedef typename LCC::Dart_const_handle Dart_const_handle;

    typedef Functor_on_merge_ Functor_on_merge;
    typedef Functor_on_split_ Functor_on_split;

    bool operator==(const My_cell_attribute_with_point& other) const
    { return Base1::operator==(other) && this->point()==other.point(); }

    bool operator!=(const My_cell_attribute_with_point& other) const
    { return !operator==(other); }

    template<typename Cellattr>
    bool operator==(const Cellattr&) const
    { return false; }

  protected:
    /// Default contructor.
    My_cell_attribute_with_point()
    {
    	index0=0;
    	distance=0;
    	prev=-1;
    }

    /// Contructor with a point in parameter.
    My_cell_attribute_with_point(const Point& apoint) : Base2(apoint)
    {
    	index0=0;
    	distance=0;
    	prev=-1;

    }
  };
} // namespace CGAL

#endif // CGAL_CELL_ATTRIBUTE_WITH_POINT_H //
// EOF //
