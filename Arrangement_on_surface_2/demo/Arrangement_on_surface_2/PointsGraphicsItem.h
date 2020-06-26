// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef POINTS_GRAPHICS_ITEM_H
#define POINTS_GRAPHICS_ITEM_H

#include <vector>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/number_utils.h>
#include <QPen>

class QPainter;
class QPen;

/**
   Add a set of points to the QGraphicsScene.
*/
class PointsGraphicsItem: public CGAL::Qt::GraphicsItem
{
public:
  PointsGraphicsItem( );

  void paint(
    QPainter* painter, const QStyleOptionGraphicsItem* option = nullptr,
    QWidget* widget = nullptr) override;

  QRectF boundingRect( ) const override;				//!< virtual function for the bounding box

  /** Template type
     *  adds the points to the vector
     */
  template < class Point >
  inline void insert( const Point& point );

  void clear( );

  void setColor( QColor c );				//!< sets the color of the curve.
  QColor getColor( ) const;					//!< returns the color of the curve

  void setPointRadius( double d );			//!< sets the user defined radius of the curve
  double getPointRadius( ) const;			//!< returns the radius of the curve

public Q_SLOTS:
  virtual void modelChanged( );

protected:
  std::vector< QPointF > points;  		/*!< vector of points of the curve */
  double pointRadius;					/*!< area of the curve draw */
  QColor color;                       	/*!< QColor object for the curve */

}; // class PointsGraphicsItem

template <class Point>
inline void PointsGraphicsItem::insert(const Point& point)
{
  this->prepareGeometryChange();

  double x = CGAL::to_double(point.x());
  double y = CGAL::to_double(point.y());
  this->points.push_back(QPointF(x, y));
}

// TODO: clean this up
#include "ArrangementTypes.h"
template <>
inline void PointsGraphicsItem::insert<Bezier_point>(const Bezier_point& point)
{
  this->prepareGeometryChange();
  auto xy = point.approximate();
  this->points.push_back(QPointF(xy.first, xy.second));
}

#endif // POINTS_GRAPHICS_ITEM_H
