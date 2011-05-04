// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// $URL: svn+ssh://rozapoga@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Minkowski_sum_with_circle_2/demo/Minkowski_sum_with_circle_2/include/LineRangeGraphicsItem.h $
// $Id: LineRangeGraphicsItem.h 50283 2009-07-01 14:14:09Z eric $
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>
//                 Eric Berberich <ericb@post.tau.ac.il>

#ifndef CGAL_QT_LINE_RANGE_GRAPHICS_ITEM_H
#define CGAL_QT_LINE_RANGE_GRAPHICS_ITEM_H

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
class LineRangeGraphicsItem : public GraphicsItem
{
  typedef typename K::Line_2 Line_2;

public:
  // TODO remove list
  LineRangeGraphicsItem(const std::list< Line_2 >* lines);
  
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

  const std::list< Line_2 > * const _m_lines;
};

template <typename K>
LineRangeGraphicsItem<K>::LineRangeGraphicsItem(const std::list< Line_2 >* lines)
  : _m_lines(lines), 
    painterostream(0)
{
  this->hide();
  setZValue(13);
}

template <typename K>
QRectF 
LineRangeGraphicsItem<K>::boundingRect() const
{
  if(scene()){
    return CGAL::Qt::viewportsBbox(scene());
  }
  return QRectF();
}

template <typename K>
void 
LineRangeGraphicsItem<K>::paint(QPainter *painter, 
                           const QStyleOptionGraphicsItem *option,
                           QWidget * widget)
{
  painter->setPen(this->Pen());
  painterostream = PainterOstream<K>(painter, boundingRect());
  for (typename std::list< Line_2 >::const_iterator lit = _m_lines->begin();
       lit != _m_lines->end(); lit++) {
    painterostream << *lit;
  }
}

template <typename K>
void 
LineRangeGraphicsItem<K>::modelChanged()
{
  update();
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_LINE_RANGE_GRAPHICS_ITEM_H

