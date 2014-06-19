#include "PointSetGraphicsItem.h"
#include "CGAL/Bbox_2.h"

QRectF PointSetGraphicsItem::boundingRect( ) const
{
  if (m_points.size() == 0)
  {
    return QRectF( );
  }

  double xmin, xmax, ymin, ymax;
  xmin = xmax = m_points[0].x();
  ymin = ymax = m_points[0].y();
  for (int i = 1; i < m_points.size(); ++i)
  {
    xmin = (xmin > m_points[i].x())? m_points[i].x() : xmin;
    xmax = (xmax < m_points[i].x())? m_points[i].x() : xmax;
    ymin = (ymin > m_points[i].y())? m_points[i].y() : ymin;
    ymax = (ymax < m_points[i].y())? m_points[i].y() : ymax;
  }
  return QRectF(xmin, ymin, xmax - xmin, ymax - ymin);
}

void PointSetGraphicsItem::paint( QPainter *painter,
  const QStyleOptionGraphicsItem *option,
  QWidget *widget )
{
  // draw scale-invariant points
  double scale = painter->worldTransform( ).m11( );
  double radius = m_pen.width( );
  painter->setBrush( m_pen.color( ) );
  radius /= scale;

  for (int i = 0; i < m_points.size(); ++i)
  {
    painter->drawEllipse( m_points[i], radius, radius );
  }
  painter->setBrush( QColor( ) );
}
