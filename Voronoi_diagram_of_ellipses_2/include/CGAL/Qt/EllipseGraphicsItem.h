//    (c) 2009 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED
    
//    (used LineGraphicsItem as a crude base)

#ifndef CGAL_QT_ELLIPSE_GRAPHICS_ITEM_H
#define CGAL_QT_ELLIPSE_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/utility.h>

#include <CGAL/Ellipse_2.h>

#include <QGraphicsScene>
#include <QGraphicsEllipseItem>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <class ET>
class EllipseGraphicsItem : public GraphicsItem
{
  typedef Ellipse_2<ET> Ellipse;

  bool is_circle;
  bool iptVisible;
  bool idVisible;

  double oldpenw;
  void update_bbox(double x0, double y0, double x1, double y1, bool force = true);

  void init() { 
      this->hide();
      oldpenw = 0.0;
      setZValue(10);
      is_circle = true;
      iptVisible = false;
      idVisible = false;
  }    

protected:

  QPen pen;
  QRectF bbox; // math box
  
  int eid;
  double xc, yc, a, b;
  double angle;

public:
  enum { Type = UserType + 1 };

  EllipseGraphicsItem() { init(); }
  EllipseGraphicsItem(const Ellipse& ellipse_) { 
      init(); setEllipse(ellipse_);
  }

  int type() const { return Type; }

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

  int get_id() const { return eid; }

  void setIptVisible(bool vis) { iptVisible = vis; }
  bool getIptVisible(bool vis) const { return iptVisible; }

  void setIdVisible(bool vis) { idVisible = vis; }
  bool getIdVisible(bool vis) const { return idVisible; }

  void setEllipse(const Ellipse& a);
};

template <class ET>
void EllipseGraphicsItem<ET>::
        update_bbox(double x0, double y0, double x1, double y1, bool force) {
  if (force || this->pen.widthF() != oldpenw) prepareGeometryChange();
  oldpenw = this->pen.widthF();
  bbox = QRectF(x0, y0, x1-x0, y1-y0);
  // BUG (on purpose): ignore i-point size, pen width
}


template <class ET>
void 
EllipseGraphicsItem<ET>::setEllipse(const Ellipse& ellipse_)
{
  xc = CGAL::to_double(ellipse_.x_center());
  yc = CGAL::to_double(ellipse_.y_center());
  a = CGAL::to_double(ellipse_.major_axis());
  b = CGAL::to_double(ellipse_.minor_axis());
  eid = ellipse_.get_id();
  double x0, x1, y0, y1;
  is_circle = ellipse_.is_circle();
  if (is_zero(ellipse_.rotation())) {
      angle = 0.0;
      x0 = xc - a;
      y0 = yc - b;
      x1 = xc + a;
      y1 = yc + b;
      // is_circle = true;
  /* } else if (is_zero(ellipse_.minor_axis())) {
      double w = CGAL::to_double(ellipse_.rotation());
      angle = atan(w)*360.0/M_PI;
      x0 = xc - a*cos(angle);
      y0 = yc - a*sin(angle);
      x1 = xc + a*cos(angle);
      y1 = yc + a*sin(angle);*/
  } else if (CGAL::abs(ellipse_.rotation()) == typename ET::QQ(1)) {
      angle = 90;
      x0 = xc - b;
      y0 = yc - a;
      x1 = xc + b;
      y1 = yc + a;
  } else {
      double w = CGAL::to_double(ellipse_.rotation());
      angle = atan(w)*360.0/M_PI;

      double t24, t16, t39, t38, t37, t23, t36, t35, t14, t21, t20, t9, t8, t12, 
             t34, t7, t33, t19, t32, t3, t31, t4, t30, t29, t28, t27, t26, t13, 
             t10, t2, t1;

      t24 = w*w;
      t16 = -1.0+t24;
      t39 = a*w;
      t38 = t16*a;
      t37 = 1.0+t24*t24;
      t23 = 1.0/(w*w);
      t36 = t23/4.0;
      t35 = -2.0*t39-yc;
      t14 = 2.0*t39;
      t21 = a*a;
      t20 = b*b;
      t9 = sqrt(t37*t20+(4.0*t21-2.0*t20)*t24);
      t8 = t14-t9;
      t12 = 1.0/t16;
      t34 = t8*t12;
      t7 = t14+t9;
      t33 = t7*t12;
      t19 = 1.0/(b*b);
      t32 = t19/(t16*t16);
      t3 = t8*t8;
      t31 = t3*t32;
      t4 = t7*t7;
      t30 = t4*t32;
      t29 = t19*t36;
      t28 = xc*t24+xc+t38;
      t27 = (t14-yc)*t32;
      t26 = t19*((t36+1.0/4.0)*xc+(-t23/4.0+1.0/4.0)*a);
      t13 = 1.0/(1.0+t24);
      t10 = sqrt(4.0*t20*t24+(-2.0*t24+t37)*t21);
      t2 = (t38-t10)*(t38-t10);
      t1 = (t38+t10)*(t38+t10);
      x0 = (-2.0*t10+t2*t26+t28)*t13/(1.0+t2*t29);
      x1 = (2.0*t10+t1*t26+t28)*t13/(1.0+t1*t29);
      y0 = -(-2.0*t33+t4*t27+(2.0*t33-yc-yc*t30)*t24+t35)*t13/(1.0+t30);
      y1 = -(-2.0*t34+t3*t27+(2.0*t34-yc-yc*t31)*t24+t35)*t13/(1.0+t31);
      if (x0 > x1) {
          double t = x0; x0 = x1; x1 = t;
          std::cerr << "x0 = " << x0 << " x1 = " << x1 << std::endl;
      }
      if (y0 > y1) {
          double t = y0; y0 = y1; y1 = t;
          std::cerr << "y0 = " << y0 << " y1 = " << y1 << std::endl;
      }
  }
/*  std::cerr << "x0 = " << x0 << " x1 = " << x1 << std::endl;
  std::cerr << "y0 = " << y0 << " y1 = " << y1 << std::endl;*/
  /*CGAL_assertion ( x0 <= x1 );
  CGAL_assertion ( y0 <= y1 );*/
  
  update_bbox(x0, y0, x1, y1);
}

template <class ET>
QRectF 
EllipseGraphicsItem<ET>::boundingRect() const
{
  if (isVisible()) return bbox;
  return QRectF();
}

template <class ET>
void 
EllipseGraphicsItem<ET>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * widget)
{
  painter->setPen(this->Pen());
  
  if (is_circle) {
      painter->drawEllipse(QPointF(xc, yc), a, a); 
  } else {
      painter->save();
      painter->translate(xc, yc);
      painter->rotate(angle);
      painter->drawEllipse(QPointF(0, 0), a, b);
      if (iptVisible) painter->drawEllipse(QPointF(-a,0),b/30.0,b/30.0);
      painter->restore();
      if (idVisible) {
          painter->save();
          painter->translate(xc, yc);
          painter->scale(1,-1);
          QRectF r2(-b/4,-b/4,b/2,b/2);
          QRectF r1 = painter->boundingRect(r2, ::Qt::AlignCenter, QString("%1").arg(eid));
          painter->scale(r2.width()/r1.width(),r2.height()/r1.height());
          painter->drawText(r1, ::Qt::AlignCenter, QString("%1").arg(eid));
          painter->restore();
      }
  }
//  painter->setPen(QPen(::Qt::black, 0, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
//  painter->drawRect(bbox);
}


template <typename CK>
void 
EllipseGraphicsItem<CK>::modelChanged()
{
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ELLIPSE_GRAPHICS_ITEM_H

