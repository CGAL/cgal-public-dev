// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>


#include "ArrangementGraphicsItem.h"

namespace CGAL {
namespace Qt {

#if 0

QRectF ArrangementGraphicsItemBase::getViewportRect( ) const
{
  QRectF clipRect;
  if ( this->scene == NULL || this->scene->views( ).size( ) == 0 )
  {
    return clipRect;
  }

  QGraphicsView* view = this->scene->views( ).first( );
  QPointF p1 = view->mapToScene( 0, 0 );
  QPointF p2 = view->mapToScene( view->width( ), view->height( ) );
  clipRect = QRectF( p1, p2 );

  return clipRect;
}

#endif

template < typename Arr_, class ArrTraits >
ArrangementGraphicsItem< Arr_, ArrTraits >::
ArrangementGraphicsItem( Arrangement* arr_ ):
  arr( arr_ ),
  painterostream( 0 )
{
  if ( this->arr->number_of_vertices( ) == 0 )
  {
    this->hide( );
  }

  this->updateBoundingBox( );
  this->setZValue( 3 );
}

template < typename Arr_, typename ArrTraits >
QRectF
ArrangementGraphicsItem< Arr_, ArrTraits >::
boundingRect( ) const
{
  QRectF rect = this->convert( this->bb );
  return rect;
}

template < typename Arr_, typename ArrTraits >
void
ArrangementGraphicsItem< Arr_, ArrTraits >::
paint(QPainter* painter,
      const QStyleOptionGraphicsItem* /* option */,
      QWidget*  /*widget*/)
{
  this->paint( painter, ArrTraits( ) );
}

template < typename Arr_, typename ArrTraits >
template < typename TTraits >
void ArrangementGraphicsItem< Arr_, ArrTraits >::
paint(QPainter* painter, TTraits /* traits */)
{
  std::cout<<"In paint ArrTraits"<<std::endl;

  painter->setPen( this->verticesPen );

  this->painterostream =
    ArrangementPainterOstream< Traits >( painter, this->boundingRect( ) );
  this->painterostream.setScene( this->scene );

  QRectF rect = this->boundingRect( );
  std::cout<<"Curve boundingRect rect\n";
  std::cout<<"left, right, bottom, top:\n";
  std::cout<<rect.left()<<", "<<rect.right()<<", "<<rect.bottom()<<", "<<rect.top()<<std::endl;

  for ( Vertex_iterator it = this->arr->vertices_begin( );
        it != this->arr->vertices_end( ); ++it )
  {
    Point_2 p = it->point( );
    Kernel_point_2 pt( p.x( ), p.y( ) );
    this->painterostream << pt;
  }

  painter->setPen( this->edgesPen );
  for ( Edge_iterator it = this->arr->edges_begin( );
        it != this->arr->edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );

    Bbox_2 bbox = curve.bbox();
    std::cout<<"Curve bounding box\n";
    std::cout<<"xmin, xmax, ymin, ymax:\n";
    std::cout<<bbox.xmin()<<", "<<bbox.xmax()<<", "<<bbox.ymin()<<", "<<bbox.ymax()<<std::endl;
    this->painterostream << curve;
  }
}

template < typename Arr_, typename ArrTraits >
template < typename CircularKernel >
void ArrangementGraphicsItem< Arr_, ArrTraits >::
paint(QPainter* painter,
      CGAL::Arr_circular_arc_traits_2< CircularKernel > /* traits */)
{
  std::cout<<"In paint Arr_circular_arc_traits_2"<<std::endl;
  typedef Kernel_point_2 Non_arc_point_2;
  typedef typename Traits::Point_2 Arc_point_2;

  painter->setPen( this->verticesPen );
  this->painterostream =
    ArrangementPainterOstream< Traits >( painter, this->boundingRect( ) );
  this->painterostream.setScene( this->scene );

  for ( Vertex_iterator it = this->arr->vertices_begin( );
        it != this->arr->vertices_end( ); ++it )
  {
    Arc_point_2 pt = it->point( );
    Non_arc_point_2 pt2(CGAL::to_double(pt.x( )), CGAL::to_double(pt.y()) );
    this->painterostream << pt2;
  }
  
  painter->setPen( this->edgesPen );
  for ( Edge_iterator it = this->arr->edges_begin( );
        it != this->arr->edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );
    this->painterostream << curve;
  }
}


template < typename Arr_, typename ArrTraits >
template < typename Coefficient_ >
void ArrangementGraphicsItem< Arr_, ArrTraits >::
paint(QPainter* painter,
      CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > /* traits */)
{
  painter->setPen( this->verticesPen );
  QRectF clipRect = this->boundingRect( );
  if ( std::isinf(clipRect.left( )) ||
       std::isinf(clipRect.right( )) ||
       std::isinf(clipRect.top( )) ||
       std::isinf(clipRect.bottom( )) )
  {
    clipRect = this->viewportRect( );
  }
  this->painterostream =
    ArrangementPainterOstream< Traits >( painter, clipRect );
  this->painterostream.setScene( this->scene );
  for ( Vertex_iterator it = this->arr->vertices_begin( );
        it != this->arr->vertices_end( ); ++it )
  {
    Point_2 p = it->point( );
    //std::pair< double, double > approx = p.to_double( );
    //Kernel_point_2 pt( approx.first, approx.second );
    //this->painterostream << pt;
    this->painterostream << p;
  }
  painter->setPen( this->edgesPen );
  for ( Edge_iterator it = this->arr->edges_begin( );
        it != this->arr->edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );
    this->painterostream << curve;
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template < typename Arr_, typename ArrTraits >
void ArrangementGraphicsItem< Arr_, ArrTraits >::updateBoundingBox( )
{
  this->updateBoundingBox( ArrTraits( ) );
}

template < typename Arr_, typename ArrTraits >
template < typename TTraits >
void ArrangementGraphicsItem< Arr_, ArrTraits >::
updateBoundingBox(TTraits /* traits */)
{
  this->prepareGeometryChange( );
  if ( this->arr->number_of_vertices( ) == 0 )
  {
    this->bb = Bbox_2( 0, 0, 0, 0 );
    this->bb_initialized = false;
    return;
  }
  else
  {
    this->bb = this->arr->vertices_begin( )->point( ).bbox( );
    this->bb_initialized = true;
  }

  for ( Curve_iterator it = this->arr->curves_begin( );
        it != this->arr->curves_end( );
        ++it )
  {
    if ( this->curveBboxMap.count( it ) == 0 )
    {
      this->curveBboxMap[ it ] = it->bbox( );
    }
    this->bb = this->bb + this->curveBboxMap[ it ];
  }
}

template < typename Arr_, typename ArrTraits >
template < typename Kernel_ >
void
ArrangementGraphicsItem< Arr_, ArrTraits >::
updateBoundingBox(CGAL::Arr_linear_traits_2< Kernel_ > /* traits */)
{
  std::cout<<"In updateBoundingBox Arr_linear_traits_2\n";
  this->prepareGeometryChange( );
  QRectF clipRect = this->viewportRect( );
  this->convert = Converter<Kernel>( clipRect );

  std::cout<<"left, right, bottom, top:\n";
  std::cout<<clipRect.left()<<", "<<clipRect.right()<<", "<<clipRect.bottom()<<", "<<clipRect.top()<<std::endl;

  if ( ! clipRect.isValid( ) /*|| this->arr->number_of_vertices( ) == 0*/ )
  {
    this->bb = Bbox_2( 0, 0, 0, 0 );
    this->bb_initialized = false;
    return;
  }
  else
  {
    this->bb = this->convert( clipRect ).bbox( );
    this->bb_initialized = true;
  }

  for ( Curve_iterator it = this->arr->curves_begin( );
        it != this->arr->curves_end( );
        ++it )
  {
    if ( it->is_segment( ) )
    {
      this->bb = this->bb + it->segment( ).bbox( );
    }
    else if ( it->is_ray( ) )
    {
      QLineF qclippedRay = this->convert( it->ray( ) );
      Segment_2 clippedRay = this->convert( qclippedRay );
      this->bb = this->bb + clippedRay.bbox( );
    }
    else // ( it->is_line( ) )
    {
      QLineF qclippedLine = this->convert( it->line( ) );
      Segment_2 clippedLine = this->convert( qclippedLine );
      this->bb = this->bb + clippedLine.bbox( );
    }
  }
}

template < typename Arr_, typename ArrTraits >
template < typename Coefficient_ >
void ArrangementGraphicsItem< Arr_, ArrTraits >::
updateBoundingBox(CGAL::Arr_algebraic_segment_traits_2<Coefficient_> traits)
{
  this->prepareGeometryChange( );
  if ( this->arr->number_of_vertices( ) == 0 )
  {
    this->bb = Bbox_2( 0, 0, 0, 0 );
    this->bb_initialized = false;
    return;
  }
  else
  {
    //std::pair< double, double > approx =
    //  this->arr->vertices_begin( )->point( ).to_double( );
    //this->bb = CGAL::Bbox_2( approx.first, approx.second,
    //                         approx.first, approx.second );
    this->bb = CGAL::Bbox_2( 0, 0, 0, 0 );
    this->bb_initialized = true;
  }
#if 0
  typename Traits::Make_x_monotone_2 make_x_monotone_2 =
    traits.make_x_monotone_2_object( );
  for ( Curve_iterator it = this->arr->curves_begin( );
        it != this->arr->curves_end( );
        ++it )
  {
    std::vector< CGAL::Object > cvs;
    make_x_monotone_2( *it, std::back_inserter( cvs ) );
    for ( unsigned int i = 0 ; i < cvs.size( ); ++i )
    {
      X_monotone_curve_2 cv;
      CGAL::assign( cv, cvs[ i ] );
      this->bb = this->bb + cv.bbox( );
    }
  }
#endif

  int curve_cnt = 0;

  for ( Edge_iterator it = this->arr->edges_begin( );
        it != this->arr->edges_end( ); ++it )
  {
    X_monotone_curve_2 curve = it->curve( );
    this->bb = this->bb + curve.bbox( );
    std::cout<<"In updateBoundingBox for"<<std::endl;
    std::cout<<curve.bbox( ).xmin()<<"\t";
    std::cout<<curve.bbox( ).xmax()<<"\t";
    std::cout<<curve.bbox( ).ymin()<<"\t";
    std::cout<<curve.bbox( ).ymax()<<"\t"<<std::endl;
    curve_cnt++;

  }

  std::cout<<"curve_cnt\t"<<curve_cnt<<std::endl;
}

template < typename Arr_, typename ArrTraits >
void ArrangementGraphicsItem< Arr_, ArrTraits >::modelChanged( )
{
  std::cout<<"In ArrangementGraphicsItem modelChanged"<<std::endl;
  if ( this->arr->is_empty( ) )
  {
    this->hide( );
  }
  else
  {
    this->show( );
  }

  std::cout<<"In ArrangementGraphicsItem modelChanged after if"<<std::endl;
  this->updateBoundingBox( );
  this->update( );
}

} // namespace Qt
} // namespace CGAL

