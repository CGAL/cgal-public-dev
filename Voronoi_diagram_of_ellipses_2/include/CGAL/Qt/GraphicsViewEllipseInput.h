//    (c) 2009 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED
    
//    (used GraphicsViewCircleInput as a crude base)

#ifndef CGAL_QT_GRAPHICS_VIEW_ELLIPSE_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_ELLIPSE_INPUT_H

#include <QGraphicsView>
#include <QRectF>
#include <QPointF>
#include <QLineF>
#include <QGraphicsItem>
#include <QGraphicsLineItem> 
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>

#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/GraphicsViewInput.h>

#include <CGAL/Qt/EllipseGraphicsItem.h> 

#include <CGAL/array.h>
#include <CGAL/simplest_rational_in_interval.h>

#include <CGAL/Ellipse_2.h>

namespace CGAL {
namespace Qt {

template <class ET>
class GraphicsViewEllipseInput : public GraphicsViewInput
{
public:
  GraphicsViewEllipseInput(QObject *parent, QGraphicsScene* s); 

  enum Ellipse_input_state { EMPTY, MAJOR_AXIS, MINOR_AXIS };

protected:
    
  virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);
  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  virtual void keyPressEvent(QKeyEvent *event);
  
  bool eventFilter(QObject *obj, QEvent *event);
  

  

private:
  void grabEllipse(QPointF qp, QPointF qq, QPointF qr, bool final);
        
  Ellipse_input_state state;
  EllipseGraphicsItem<ET> *ell;
  QGraphicsLineItem *qline;
  QPointF qp, qq, qr;
  QGraphicsScene *scene_;  
};


template <class ET>
GraphicsViewEllipseInput<ET>::GraphicsViewEllipseInput(QObject *parent, QGraphicsScene* s)
  : GraphicsViewInput(parent), 
    state(EMPTY), ell(new EllipseGraphicsItem<ET>()), 
    qline(new QGraphicsLineItem()), scene_(s)
{
  qline->hide();
  s->addItem(qline);
  s->addItem(ell);
}


template<class QQ> inline QQ rationalize(double x, double res) {
    return simplest_rational_in_interval<QQ>(x-res/2.0, x+res/2.0);
}

template <class ET>
void 
GraphicsViewEllipseInput<ET>::grabEllipse(QPointF qp, QPointF qq, QPointF qr,
        bool final = false)
{
    double res_x = fabs(qq.x() - qp.x())/1000.0;
    double res_y = fabs(qq.y() - qp.y())/1000.0;
    double res = sqrt(res_x*res_x + res_y+res_y);

    double a, b, w, xc, yc;
    w = tan(atan((qq.y()-qp.y())/(qq.x()-qp.x())) / 2.0);
    xc = (qp.x() + qq.x()) / 2.0;
    yc = (qp.y() + qq.y()) / 2.0;
    a = sqrt((qp.x()-xc)*(qp.x()-xc) + (qp.y()-yc)*(qp.y()-yc));

    b = fabs((qq.x()-qp.x())*(qp.y()-qr.y()) - (qp.x()-qr.x())*(qq.y()-qp.y()))/
       sqrt((qq.x()-qp.x())*(qq.x()-qp.x()) + (qq.y()-qp.y())*(qq.y()-qp.y()));

    typename ET::QQ ra, rb, rw, rxc, ryc;

    ra = rationalize<typename ET::QQ>(a,res);
    if (ra < typename ET::QQ(1,10000)) ra = typename ET::QQ(1,10000);
    rb = rationalize<typename ET::QQ>(b,res);
    if (rb < typename ET::QQ(1,10000)) rb = typename ET::QQ(1,10000);
    rw = rationalize<typename ET::QQ>(w,0.01);
    rxc = rationalize<typename ET::QQ>(xc,res);
    ryc = rationalize<typename ET::QQ>(yc,res);
    Ellipse_2<ET> e(ra, rb, rw, rxc, ryc);
    ell->setEllipse(e);
    //std::cerr << e << std::endl;
    //std::cerr << "a = " << a << " b = " << b << 
    //             " w = " << w << " xc = " << xc << " yc = " << yc << std::endl;
    if (final) {
        emit generate(CGAL::make_object(e));
    }
}

template <class ET>
void 
GraphicsViewEllipseInput<ET>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{ 
  if(event->modifiers()  & ::Qt::ShiftModifier){
    return;
  }
  if (event->button() == ::Qt::RightButton) {
      state = EMPTY;
      QRectF r(qline->boundingRect() | ell->boundingRect());
      qline->hide();
      ell->hide();
      scene_->update(r);
      return;
  }
  if (state == EMPTY) {
      qp = event->scenePos();
      qq = qp;
      qline->setLine(QLineF(qp, qq));
      qline->show();
      state = MAJOR_AXIS;
  } else if (state == MAJOR_AXIS) {
      qq = event->scenePos();
      qline->setLine(QLineF(qp, qq));
      qline->hide();
      state = MINOR_AXIS;
      qr = qq;
      grabEllipse(qp, qq, qr);
      ell->show();
  } else if (state == MINOR_AXIS) {
      qr = event->scenePos();
      state = EMPTY;
      grabEllipse(qp, qq, qr, true);
      ell->hide();
  }
}


template <class ET>
void 
GraphicsViewEllipseInput<ET>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if(state == EMPTY){
    return;
  } else if(state == MAJOR_AXIS) {
    qq = event->scenePos();
    qline->setLine(QLineF(qp, qq));
  } else if (state == MINOR_AXIS) {
      qr = event->scenePos();
      grabEllipse(qp, qq, qr);
  }
}

template <class ET>
void 
GraphicsViewEllipseInput<ET>::keyPressEvent ( QKeyEvent * event ) 
{
  if(event->key() == ::Qt::Key_Delete){
      state = EMPTY;
      QRectF r(qline->boundingRect() | ell->boundingRect());
      qline->hide();
      ell->hide();
      scene_->update(r);
  }
}


template <class ET>
bool 
GraphicsViewEllipseInput<ET>::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMousePress) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mousePressEvent(mouseEvent);
    return true;
  } else if (event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return true;
  } else if (event->type() == QEvent::KeyPress) {
    QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
    keyPressEvent(keyEvent);
    return true;
  } else{
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
} 

} // namespace Qt

} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_ELLIPSE_INPUT_H
