#ifndef DCEL_TYPES_H
#define DCEL_TYPES_H
#include <CGAL/Arr_default_dcel.h>
#include <QColor>

class Face_with_color : public CGAL::Arr_face_base
{
  QColor    m_color;
  bool      m_visited;

public:
  Face_with_color() : CGAL::Arr_face_base(), m_color(), m_visited(false) { }

  QColor color() const { return m_color; }
  void set_color(const QColor& c) { m_color = c; }
  bool visited() const{ return m_visited; }
  void set_visited(bool b) { m_visited = b; }
};

template <class Traits>
class Dcel :
  public CGAL::Arr_dcel_base<CGAL::Arr_vertex_base<typename Traits::Point_2>,
                             CGAL::Arr_halfedge_base<
                               typename Traits::X_monotone_curve_2>,
                             Face_with_color>
{
public:
   /*! \struct
   * An auxiliary structure for rebinding the DCEL with a new traits class.
   */
  template <typename T>
  struct rebind
  {
    typedef Dcel<T> other;
  };

  // CREATION
  Dcel() {}
};

#endif // DCEL_TYPES_H
