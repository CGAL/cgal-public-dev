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
// $URL: svn+ssh://rozapoga@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Minkowski_sum_with_circle_2/demo/Minkowski_sum_with_circle_2/include/SlopedPolygonGraphicsItem.h $
// $Id: SlopedPolygonGraphicsItem.h 50321 2009-07-02 10:18:26Z eric $
// 
//
// Author(s)     : Eric Berberich <ericb@post.tau.ac.il>

#ifndef CGAL_QT_SLOPED_POLYGON_GRAPHICS_ITEM_H
#define CGAL_QT_SLOPED_POLYGON_GRAPHICS_ITEM_H

#include <CGAL/Qt/PolygonGraphicsItem.h>

namespace CGAL {
namespace Qt {

template < typename P, typename CircleApproximation_2 >
class SlopedPolygonGraphicsItem : public PolygonGraphicsItem< P >
{
  typedef CircleApproximation_2 Circle_approximation_2;

  typedef PolygonGraphicsItem< P > Base;
  
  typedef typename P::Traits Traits;

public:

  SlopedPolygonGraphicsItem(const P* p_, const Circle_approximation_2* ca_, QColor& circle_color_);

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
private:

  const Circle_approximation_2 * const ca;
  QColor circle_color;
};


template <typename P, typename CircleApproximation_2>
SlopedPolygonGraphicsItem<P, CircleApproximation_2>::SlopedPolygonGraphicsItem(
    const P * p_, const Circle_approximation_2* ca_, QColor& circle_color_
)
  :  Base(p_), ca(ca_), circle_color(circle_color_)
{}

template <typename P, typename CircleApproximation_2>
void 
SlopedPolygonGraphicsItem<P, CircleApproximation_2>::paint(
    QPainter *painter, 
    const QStyleOptionGraphicsItem *option,
    QWidget * widget)
{
  painter->setPen(this->edgesPen());
  Base::painterostream = PainterOstream<Traits>(painter);
  
  QPen cpen(circle_color, 0, ::Qt::DashLine, 
            ::Qt::RoundCap, ::Qt::RoundJoin);
  
  if (this->drawEdges()) {
    for(typename P::Edge_const_iterator eit = Base::poly->edges_begin();
        eit != Base::poly->edges_end();
        ++eit) {
      if (ca->is_circle_slope(*eit)) {
        painter->setPen(cpen);
        Base::painterostream << *eit;
        painter->setPen(this->edgesPen());
      } else { 
        Base::painterostream << *eit;
      }
    }
  }
  
  if (this->drawVertices()) {
    Converter<Traits> convert;
    
    painter->setPen(this->verticesPen());
    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
    for(typename P::Vertex_iterator it = Base::poly->vertices_begin();
        it != Base::poly->vertices_end();
        it++){
      QPointF point = matrix.map(convert(*it));
      painter->drawPoint(point);
    }
  }
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_SLOPED_POLYGON_GRAPHICS_ITEM_H
