// Copyright (c) 2009  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://rozapoga@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Minkowski_sum_with_circle_2/demo/Minkowski_sum_with_circle_2/include/CircleRangeGraphicsItem.h $
// $Id: CircleRangeGraphicsItem.h 50309 2009-07-01 21:43:00Z eric $
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>
//                 Eric Berberich <ericb@post.tau.ac.il>

#ifndef CGAL_QT_CIRCLE_RANGE_GRAPHICS_ITEM_H
#define CGAL_QT_CIRCLE_RANGE_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/utility.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename K>
class CircleRangeGraphicsItem : public GraphicsItem
{
  typedef typename K::Circle_2 Circle_2;

public:
  // TODO remove list
  CircleRangeGraphicsItem(const std::list< Circle_2 >* circles);
  
  void modelChanged();

public:
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, 
             const QStyleOptionGraphicsItem *option, 
             QWidget *widget);
  

  const QPen& Pen() const
  {
    return this->pen;
  }


  void setPen(const QPen& pen_)
  {
    this->pen = pen_;
  }

protected:

  QPainter* m_painter;
  PainterOstream<K> painterostream;

  QPen pen;

  const std::list< Circle_2 > * const _m_circles;
};

template <typename K>
CircleRangeGraphicsItem<K>::CircleRangeGraphicsItem(
    const std::list< Circle_2 >* circles)
  : _m_circles(circles), 
    painterostream(0)
{
  this->hide();
  setZValue(13);
}

template <typename K>
QRectF 
CircleRangeGraphicsItem<K>::boundingRect() const
{
  if(scene()){
    return CGAL::Qt::viewportsBbox(scene());
  }
  return QRectF();
}

template <typename K>
void 
CircleRangeGraphicsItem<K>::paint(QPainter *painter, 
                           const QStyleOptionGraphicsItem *option,
                           QWidget * widget)
{
  painter->setPen(this->Pen());
  LOG_DEBUG << "painting " << _m_circles->size() << " circles." << std::endl;
  painterostream = PainterOstream<K>(painter, boundingRect());
  for (typename std::list< Circle_2 >::const_iterator cit = 
         _m_circles->begin();
       cit != _m_circles->end(); cit++) {
    painterostream << *cit;
  }
}

template <typename K>
void 
CircleRangeGraphicsItem<K>::modelChanged()
{
  update();
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_CIRCLE_RANGE_GRAPHICS_ITEM_H

