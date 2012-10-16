//    (c) 2009 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED
    
//  (used LineGraphicsItem as a crude base)

#ifndef CGAL_QT_ELLIPSE_BISECTOR_GRAPHICS_ITEM_H
#define CGAL_QT_ELLIPSE_BISECTOR_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/EllipseGraphicsItem.h>
#include <CGAL/Qt/utility.h>

#include <CGAL/Ellipse_bisector_2.h>

#include <QGraphicsScene>
#include <QGraphicsLineItem>
#include <QPainter>
#include <QPen>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <class ET>
class EllipseBisectorGraphicsItem : public GraphicsItem
{
  typedef Ellipse_bisector_2<ET> Ellipse_bisector;

  double oldpenw;

  void update_bbox(double x0, double y0, double x1, double y1, bool force = true);
  void init() { 
      this->hide();
      oldpenw = 0.0;
      setZValue(5);
  }  
    
public:
  EllipseBisectorGraphicsItem() { init(); }
  EllipseBisectorGraphicsItem(const Ellipse_bisector& bisec_, 
          const QPen pen_ = QPen()) {
      init(); 
      setEllipseBisector(bisec_);
      setPen(pen_);
      show();
  }

  void modelChanged();
  

public:
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  

  const QPen& Pen() const
  {
    return this->pen;
  }


  void setPen(const QPen& pen_)
  {
    this->pen = pen_;
    update_bbox(bbox.x(), bbox.y(), bbox.right(), bbox.bottom(), false);
  }

  void setEllipseBisector(const Ellipse_bisector& a);
//  void setEllipseBisector(const Ellipse_bisector& a, bool clip = true);

protected:

  QPen pen;
  QRectF bbox;

  QPolygonF poly;
};

template <class ET>
void EllipseBisectorGraphicsItem<ET>::
        update_bbox(double x0, double y0, double x1, double y1, bool force) {
  if (force || this->pen.widthF() != oldpenw) prepareGeometryChange();
  oldpenw = this->pen.widthF();
  bbox = QRectF(x0, y0, x1-x0, y1-y0);
  // BUG (on purpose): ignore pen width
}


template <class ET>
void 
EllipseBisectorGraphicsItem<ET>::
        setEllipseBisector(const Ellipse_bisector& bisec_)
{
  QVector<QPointF> vec;
  for (int i = 0; i < bisec_.xx.size(); i++) {
      vec.push_back(QPointF(bisec_.xx[i], bisec_.yy[i]));
  }
  poly = QPolygonF(vec);
  QRectF b = poly.boundingRect();
/*  if (clip) {
      EllipseGraphicsItem<ET> g1(bisec_.get_bisector_curve_apx().get_e1());
      EllipseGraphicsItem<ET> g2(bisec_.get_bisector_curve_apx().get_e2());
      b = b & (g1.boundingRect() | g2.boundingRect());
  }*/
  update_bbox(b.x(), b.y(), b.right(), b.bottom());
}

template <class ET>
QRectF 
EllipseBisectorGraphicsItem<ET>::boundingRect() const
{
  if (isVisible()) return bbox;
  return QRectF();
}

template <class ET>
void 
EllipseBisectorGraphicsItem<ET>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * widget)
{
  painter->setPen(this->Pen());
  
  painter->drawPolyline(poly); 
//  painter->setPen(QPen(::Qt::black, 0, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
//  painter->drawRect(bbox);
}


template <typename CK>
void 
EllipseBisectorGraphicsItem<CK>::modelChanged()
{
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ELLIPSE_BISECTOR_GRAPHICS_ITEM_H

