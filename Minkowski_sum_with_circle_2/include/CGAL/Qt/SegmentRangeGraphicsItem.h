#ifndef CGAL_QT_SEGMENT_RANGE_GRAPHICS_ITEM_H
#define CGAL_QT_SEGMENT_RANGE_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/ArrPainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/utility.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename K>
class SegmentRangeGraphicsItem : public GraphicsItem
{
  typedef ArrPainterOstream<K> PainterOstreamType;
  typedef typename PainterOstreamType::TraitsHelper::Segment_2 Segment_2;

public:

  SegmentRangeGraphicsItem(const std::list< Segment_2 >* lines);
  
  void modelChanged();

public:
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, 
             const QStyleOptionGraphicsItem *option, 
             QWidget *widget);
  

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

protected:

  QPainter* m_painter;
  ArrPainterOstream<K> painterostream;

  QPen edges_pen;
  QPen vertices_pen;

  const std::list< Segment_2 > * const m_segments;
};

template <typename K>
SegmentRangeGraphicsItem<K>::SegmentRangeGraphicsItem(const std::list< Segment_2 >* lines)
  : m_segments(lines), 
    painterostream(0)
{
  this->hide();
  setZValue(13);
}

template <typename K>
QRectF 
SegmentRangeGraphicsItem<K>::boundingRect() const
{
  if(scene()){
    return CGAL::Qt::viewportsBbox(scene());
  }
  return QRectF();
}

template <typename K>
void 
SegmentRangeGraphicsItem<K>::paint(QPainter *painter, 
                           const QStyleOptionGraphicsItem *option,
                           QWidget * widget)
{
  painter->setPen(edges_pen);
  painterostream = PainterOstreamType(painter, boundingRect());
  for (typename std::list< Segment_2 >::const_iterator sit = m_segments->begin();
       sit != m_segments->end(); sit++) {
    painterostream << *sit;
  }
}

template <typename K>
void 
SegmentRangeGraphicsItem<K>::modelChanged()
{
  update();
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_SEGMENT_RANGE_GRAPHICS_ITEM_H

