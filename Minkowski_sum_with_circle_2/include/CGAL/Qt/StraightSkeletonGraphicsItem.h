#ifndef CGAL_QT_STRAIGHT_SKELETON_GRAPHICS_ITEM_H
#define CGAL_QT_STRAIGHT_SKELETON_GRAPHICS_ITEM_H

#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <CGAL/Bbox_2.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>
#include <boost/shared_ptr.hpp>

namespace CGAL {
namespace Qt {

template <typename S>
class StraightSkeletonGraphicsItem : public GraphicsItem
{
  typedef S Ss ;
  typedef boost::shared_ptr<Ss> Ss_ptr_2;  
  typedef typename Ss::Traits Traits;
  typedef typename Traits::Segment_2 Segment_2;
  typedef typename Traits::Point_2 Point_2;
public:
  StraightSkeletonGraphicsItem(const Ss_ptr_2& ss_ptr);

  void modelChanged();

public:
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  const QPen& verticesPen() const
  {
    return vertices_pen;
  }

  const QPen& edgesPen() const
  {
    return edges_pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    vertices_pen = pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    edges_pen = pen;
  }

  QRectF boundingRect() const;

  void updateBoundingBox();

protected:

  const Ss_ptr_2& m_ss_ptr;
  QPainter* m_painter;

  PainterOstream<Traits> painterostream;
  QRectF bounding_rect;

  QPen edges_pen;
  QPen vertices_pen;
};


template <typename S>
StraightSkeletonGraphicsItem<S>::StraightSkeletonGraphicsItem(const Ss_ptr_2& ss_ptr)
  :  m_ss_ptr(ss_ptr), painterostream(0)  
{
  if(!(m_ss_ptr) || (m_ss_ptr->size_of_halfedges() == 0) ){
    this->hide();
  }
  updateBoundingBox();
  setZValue(3);
}


template <typename S>
void 
StraightSkeletonGraphicsItem<S>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * widget)
{
  painter->setPen(this->edgesPen());
  painterostream = PainterOstream<Traits>(painter);
  
  if(m_ss_ptr)
  {
    typedef typename Ss::Halfedge_const_iterator Halfedge_const_iterator ;

    for ( Halfedge_const_iterator hit = m_ss_ptr->halfedges_begin(); 
      hit != m_ss_ptr->halfedges_end(); 
      ++hit )
    {
      if ( hit->is_bisector() ){
        Segment_2 seg(hit->opposite()->vertex()->point(), hit->vertex()->point());
        painterostream << seg;
      }
    }
  }
}

template <typename S>
void 
StraightSkeletonGraphicsItem<S>::updateBoundingBox()
{
  Converter<Traits> convert;

  struct local
  {
    static CGAL::Bbox_2 ConstructBbox_2(const Point_2& point)
    {
      double x = to_double(point.x());
      double y = to_double(point.y());
      bool is_x_neg = (x < 0);
      bool is_y_neg = (y < 0);
      return CGAL::Bbox_2((is_x_neg? x: 0), (is_y_neg? y: 0), (is_x_neg? 0: x), (is_y_neg? 0: y));
    }
  };

  if(m_ss_ptr)
  {
    CGAL::Bbox_2 skel_box(0, 0, 0, 0);
    typedef typename Ss::Vertex_const_iterator Vertex_const_iterator ;

    for ( Vertex_const_iterator vit = m_ss_ptr->vertices_begin(); 
      vit != m_ss_ptr->vertices_end(); 
      ++vit )
    {
      Point_2 vertex = vit->point();
      skel_box = skel_box + local::ConstructBbox_2(vertex);
    }
    
    bounding_rect = convert(skel_box);
  }
}


template <typename S>
QRectF 
StraightSkeletonGraphicsItem<S>::boundingRect() const
{
  return bounding_rect;
}

template <typename S>
void 
StraightSkeletonGraphicsItem<S>::modelChanged()
{
  if(!(m_ss_ptr) || (m_ss_ptr->size_of_halfedges() == 0) ){
    this->hide();
  } else if((m_ss_ptr->size_of_halfedges() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_STRAIGHT_SKELETON_GRAPHICS_ITEM_H
