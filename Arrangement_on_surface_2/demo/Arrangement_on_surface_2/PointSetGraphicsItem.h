#ifndef POINT_SET_GRAPHICS_ITEM_H
#define POINT_SET_GRAPHICS_ITEM_H
#include <QGraphicsItem>
#include <QRectF>
#include <QPainter>

class QStyleOptionGraphicsItem;
class QWidget;

class PointSetGraphicsItem : public QGraphicsItem
{

protected:
  std::vector< QPointF > m_points;

  QPen m_pen;

public:
  template < class InputIterator >
  PointSetGraphicsItem( InputIterator begin, InputIterator end ):
    QGraphicsItem( 0 ),
    m_pen( QPen( ::Qt::red, 3 ) )
  {
    std::copy( begin, end, std::back_inserter( m_points ) );
    setZValue( 5 );
  }

  virtual QRectF boundingRect( ) const;
  virtual void paint( QPainter *painter,
    const QStyleOptionGraphicsItem *option,
    QWidget *widget );
};

#endif // POINT_SET_GRAPHICS_ITEM_H
