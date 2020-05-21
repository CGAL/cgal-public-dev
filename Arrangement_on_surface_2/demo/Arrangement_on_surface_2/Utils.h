// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>

#ifndef CGAL_ARRANGEMENTS_DEMO_UTILS_H
#define CGAL_ARRANGEMENTS_DEMO_UTILS_H

#include <CGAL/iterator.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_circular_arc_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

#include "ArrangementTypes.h"

class QGraphicsScene;

class QGraphicsSceneMixin
{
public:
  /*! Costructor */
  QGraphicsSceneMixin( ) : scene( 0 ) { }

  /*! Destructor (virtual) */
  virtual ~QGraphicsSceneMixin() {}

  virtual void setScene( QGraphicsScene* scene_ ) { this->scene = scene_; }

  virtual QGraphicsScene* getScene( ) const { return this->scene; }

  virtual QRectF viewportRect( ) const
  {
    // std::cout<<"In QGraphicsSceneMixin viewportRect\n";
    QRectF res;
    if ( this->scene == NULL )
    {
      return res;
    }

    QList< QGraphicsView* > views = this->scene->views( );
    if ( views.size( ) == 0 )
    {
      // std::cout<<"Return: views.size( ) == 0\n";
      return res;
    }
    // assumes the first view is the right one
    QGraphicsView* viewport = views.first( );
    QPointF p1 = viewport->mapToScene( 0, 0 );
    QPointF p2 = viewport->mapToScene(viewport->width(), viewport->height());

    double xmin = (std::min)(p1.x(), p2.x());
    double xmax = (std::max)(p1.x(), p2.x());
    double ymin = (std::min)(p1.y(), p2.y());
    double ymax = (std::max)(p1.y(), p2.y());

    res = QRectF( QPointF( xmin, ymin ), QPointF( xmax, ymax ) );

    // std::cout<<"Return: OKAY\n";
    return res;
  }

  QPoint fromScene( QPointF p, bool* ok = 0 )
  {
    QPoint res;
    if ( this->scene == NULL )
    {
      if ( ok ) { *ok = false; }
      return res;
    }
    QList< QGraphicsView* > views = this->scene->views( );
    if ( views.size( ) == 0 )
    {
      if ( ok ) { *ok = false; }
      return res;
    }

    // assumes the first view is the right one
    QGraphicsView* viewport = views.first( );

    if ( ok ) { *ok = true; }
    res = viewport->mapFromScene( p );
    return res;
  }

  QPointF toScene( QPoint p, bool* ok = 0 )
  {
    QPointF res;
    if ( this->scene == NULL )
    {
      if ( ok ) { *ok = false; }
      return res;
    }
    QList< QGraphicsView* > views = this->scene->views( );
    if ( views.size( ) == 0 )
    {
      if ( ok ) { *ok = false; }
      return res;
    }

    // assumes the first view is the right one
    QGraphicsView* viewport = views.first( );

    if ( ok ) { *ok = true; }
    res = viewport->mapToScene( p );
    return res;
  }

  int fromScene( double d, bool* ok = 0 )
  {
    QPointF p( d, 0 );
    QPoint pp = this->fromScene( p, ok );
    return pp.x( );
  }

  double toScene( int i, bool* ok = 0 )
  {
    QPoint p( i, 0 );
    QPointF pp = this->toScene( p, ok );
    return pp.x( );
  }

protected: // fields
  QGraphicsScene* scene;
};

BOOST_MPL_HAS_XXX_TRAIT_DEF( Approximate_2 )

template <typename Arr_, bool b = has_Approximate_2< Arr_ >::value>
struct Supports_landmarks
{
  typedef CGAL::Boolean_tag< b > Tag;
  struct LandmarksType { };
};

template <typename Arr_>
struct Supports_landmarks< Arr_, true >
{
  typedef CGAL::Tag_true Tag;
  typedef CGAL::Arr_landmarks_point_location< Arr_ > LandmarksType;
};

/**
   Support for new ArrTraits should specify types:

   * Kernel - a not-necessarily-exact kernel to represent the arrangement
   graphically. We'll use the Point_2 type provided by this kernel for
   computing distances
   * Point_2 - the point type used in the particular arrangement
   * CoordinateType - the coordinate type used by the point type
   */
template <typename ArrTraits>
class ArrTraitsAdaptor
{ };

template <typename Kernel_>
class ArrTraitsAdaptor< CGAL::Arr_segment_traits_2< Kernel_ > >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2< Kernel > ArrTraits;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename Kernel_>
class ArrTraitsAdaptor< CGAL::Arr_linear_traits_2< Kernel_ > >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_linear_traits_2< Kernel > ArrTraits;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename SegmentTraits>
class ArrTraitsAdaptor< CGAL::Arr_polyline_traits_2< SegmentTraits > >
{
public:
  typedef CGAL::Arr_polyline_traits_2< SegmentTraits > ArrTraits;
  typedef typename SegmentTraits::Kernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename CircularKernel >
class ArrTraitsAdaptor< CGAL::Arr_circular_arc_traits_2< CircularKernel > >
{
public:
  typedef CGAL::Arr_circular_arc_traits_2< CircularKernel > ArrTraits;
  typedef CircularKernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::Root_of_2 CoordinateType;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits >
class ArrTraitsAdaptor< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel,
                                                  NtTraits > >
{
public:
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > ArrTraits;
  typedef AlgKernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class ArrTraitsAdaptor< CGAL::Arr_Bezier_curve_traits_2< RatKernel, AlgKernel,
                                                         NtTraits > >
{
public:
  typedef CGAL::Arr_Bezier_curve_traits_2< RatKernel, AlgKernel, NtTraits > ArrTraits;
  typedef AlgKernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;
  typedef typename Kernel::FT CoordinateType;
};

template <typename Coefficient_>
class ArrTraitsAdaptor< CGAL::Arr_algebraic_segment_traits_2< Coefficient_ > >
{
public:
  typedef Coefficient_ Coefficient;
  typedef typename CGAL::Arr_algebraic_segment_traits_2<Coefficient>
                                                        ArrTraits;
  typedef typename ArrTraits::Point_2                   Point_2; // CKvA_2
  typedef typename ArrTraits::Algebraic_real_1          CoordinateType;
  typedef CGAL::Cartesian< typename ArrTraits::Bound >  Kernel;
  //typedef typename ArrTraits::CKvA_2                  Kernel;
};

template <typename ArrTraits>
class Arr_compute_y_at_x_2 : public QGraphicsSceneMixin
{
public:
  typedef ArrTraits Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename ArrTraitsAdaptor< Traits >::CoordinateType CoordinateType;
  // typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef std::pair< typename Traits::Point_2, Multiplicity >
                                                        IntersectionResult;

  /*! Constructor */
  Arr_compute_y_at_x_2( ) :
    intersectCurves( this->traits.intersect_2_object( ) )
  { }

  /*! Destructor (virtual) */
  virtual ~Arr_compute_y_at_x_2() {}

  CoordinateType operator() ( const X_monotone_curve_2& curve,
                              const CoordinateType& x )
  {
    typename Traits::Left_side_category category;
    return this->operator()( curve, x, this->traits, category );
  }

  double approx( const X_monotone_curve_2& curve, const CoordinateType& x )
  {
    return CGAL::to_double( (*this)( curve, x ) );
  }

protected:
  template <typename TTraits >
  CoordinateType operator() ( const X_monotone_curve_2& curve,
                              const CoordinateType& x, TTraits traits_,
                              CGAL::Arr_oblivious_side_tag )
  {
    typedef typename TTraits::Construct_x_monotone_curve_2
      Construct_x_monotone_curve_2;
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2 =
      traits_.construct_x_monotone_curve_2_object( );
    CoordinateType res( 0 );
    CGAL::Bbox_2 clipRect = curve.bbox( );
    Point_2 p1c1( x, CoordinateType( clipRect.ymin( ) - 1 ) ); // clicked point
    // upper bounding box
    Point_2 p2c1( x, CoordinateType( clipRect.ymax( ) + 1 ) );

    const X_monotone_curve_2 verticalLine =
      construct_x_monotone_curve_2( p1c1, p2c1 );
    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );

    this->intersectCurves( curve, verticalLine, oi );

    IntersectionResult pair;
    if ( CGAL::assign( pair, o ) )
    {
      Point_2 pt = pair.first;
      res = pt.y( );
    }
    return res;
  }

  template <typename TTraits >
  CoordinateType operator() ( const X_monotone_curve_2& curve,
                              const CoordinateType& x, TTraits traits_,
                              CGAL::Arr_open_side_tag )
  {
    typename TTraits::Construct_x_monotone_curve_2
      construct_x_monotone_curve_2 =
      traits_.construct_x_monotone_curve_2_object( );
    CoordinateType res( 0 );
    // QRectF clipRect = this->viewportRect( );
    Line_2 line = curve.supporting_line( );
    // FIXME: get a better bounding box for an unbounded segment
    Point_2 p1c1( x, CoordinateType( -10000000 ) ); // clicked point
    Point_2 p2c1( x, CoordinateType(  10000000 ) ); // upper bounding box

    const X_monotone_curve_2 verticalLine =
      construct_x_monotone_curve_2( p1c1, p2c1 );
    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );

    this->intersectCurves( curve, verticalLine, oi );

    IntersectionResult pair;
    if ( CGAL::assign( pair, o ) )
    {
      Point_2 pt = pair.first;
      res = pt.y( );
    }
    return res;
  }

protected:
  Traits traits;
  Intersect_2 intersectCurves;
};

template <typename CircularKernel>
class Arr_compute_y_at_x_2< CGAL::Arr_circular_arc_traits_2<CircularKernel> > :
  public QGraphicsSceneMixin
{
public:
  typedef CGAL::Arr_circular_arc_traits_2< CircularKernel > Traits;
  typedef CircularKernel                                Kernel;
  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::Root_of_2                    Root_of_2;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Traits::Point_2                      Arc_point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Kernel::Line_arc_2                   Line_arc_2;
  // Circular_arc_2
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef std::pair< typename Traits::Point_2, Multiplicity >
                                                        IntersectionResult;

  /*! Constructor */
  Arr_compute_y_at_x_2( ) :
    intersectCurves( this->traits.intersect_2_object( ) )
  { }

  /*! Destructor (virtual) */
  virtual ~Arr_compute_y_at_x_2() {}

  Root_of_2 operator() ( const X_monotone_curve_2& curve, const FT& x )
  {
    Root_of_2 res( 0 );
    CGAL::Bbox_2 clipRect = curve.bbox( );
    Point_2 p1c1( x, FT( clipRect.ymin( ) - 1 ) ); // clicked point
    Point_2 p2c1( x, FT( clipRect.ymax( ) + 1 ) ); // upper bounding box
    Line_arc_2 verticalLine( Segment_2( p1c1, p2c1 ) );

    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );

    this->intersectCurves( curve, verticalLine, oi );

    IntersectionResult pair;
    if ( CGAL::assign( pair, o ) )
    {
      Arc_point_2 pt = pair.first;
      res = pt.y( );
    }
    return res;
  }

  double approx( const X_monotone_curve_2& curve, const FT& x )
  {
    return CGAL::to_double( (*this)( curve, x ) );
  }

  // FIXME: inexact projection
  Root_of_2 operator() ( const X_monotone_curve_2& curve, const Root_of_2& x )
  {
    FT approx( CGAL::to_double( x ) );
    return this->operator()( curve, approx );
  }

  double approx( const X_monotone_curve_2& curve, const Root_of_2& x )
  {
    return CGAL::to_double( (*this)( curve, x ) );
  }

protected:
  Traits traits;
  Intersect_2 intersectCurves;
};

template <typename Coefficient_>
class Arr_compute_y_at_x_2< CGAL::Arr_algebraic_segment_traits_2<
                              Coefficient_ > > : public QGraphicsSceneMixin
{
public:
  typedef Coefficient_ Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2< Coefficient > Traits;
  typedef typename Traits::Algebraic_real_1             CoordinateType;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

  CoordinateType operator() ( const X_monotone_curve_2& curve,
                              const CoordinateType& x )
  {
    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );
    Intersect_2 intersect = traits.intersect_2_object( );
    X_monotone_curve_2 c2 = this->makeVerticalLine( x );
    intersect( curve, c2, oi );
    std::pair< Point_2, Multiplicity > res;
    if ( CGAL::assign( res, o ) ) // TODO: handle failure case
    {
      Point_2 p = res.first;
      // std::cout << "approx y: " << p.to_double( ).second << std::endl;
      CoordinateType coord = p.y( );
      return coord;
    }
    else
    {
      std::cout << "Warning: vertical projection failed" << std::endl;
      return CoordinateType( 0 );
    }
  }

  double approx( const X_monotone_curve_2& curve, const CoordinateType& x )
  {
    CGAL::Object o;
    CGAL::Oneset_iterator< CGAL::Object > oi( o );
    Intersect_2 intersect = traits.intersect_2_object( );
    X_monotone_curve_2 c2 = this->makeVerticalLine( x );
    intersect( curve, c2, oi );
    std::pair< Point_2, Multiplicity > res;
    if ( CGAL::assign( res, o ) ) // TODO: handle failure case
    {
      Point_2 p = res.first;
      std::pair< double, double > tmp = p.to_double();
      return tmp.second;
    }
    else
    {
      return 0;
    }
  }

protected:
  X_monotone_curve_2 makeVerticalLine( const CoordinateType& x )
  {
    typename Traits::Construct_point_2 constructPoint =
      traits.construct_point_2_object( );
    typename Traits::Construct_x_monotone_segment_2 constructSegment =
      traits.construct_x_monotone_segment_2_object( );

    std::vector< X_monotone_curve_2 > curves;
    Point_2 p1 = constructPoint( x, CoordinateType( -1000000 ) );
    Point_2 p2 = constructPoint( x, CoordinateType( +1000000 ) );
    constructSegment( p1, p2, std::back_inserter( curves ) );
    return curves[ 0 ]; // by construction, there is one curve in curves
  }
  Traits traits;
};

template < class ArrTraits >
class Compute_squared_distance_2_base : public QGraphicsSceneMixin
{
public:
  typedef CGAL::Cartesian< double > InexactKernel;

public:
  // ctors
  Compute_squared_distance_2_base( ) { }

public: // methods

  template < class T1, class T2 >
  double operator() ( const T1& t1, const T2& t2 )
  {
    return this->squaredDistance( t1, t2 );
  }

protected: // fields
  typename Kernel::Compute_squared_distance_2 squared_distance;
  InexactKernel::Compute_squared_distance_2 squaredDistance;
};

template <typename ArrTraits >
class Compute_squared_distance_2 :
  public Compute_squared_distance_2_base< ArrTraits >
{ };

template <typename Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_segment_traits_2< Kernel_ > > :
  public Compute_squared_distance_2_base<CGAL::Arr_segment_traits_2<Kernel_> >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2< Kernel > Traits;
  typedef Compute_squared_distance_2_base< Traits > Superclass;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_linear_traits_2< Kernel_ > > :
  public Compute_squared_distance_2_base<CGAL::Arr_linear_traits_2<Kernel_> >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_linear_traits_2< Kernel > Traits;
  typedef Compute_squared_distance_2_base< Traits > Superclass;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Ray_2 Ray_2;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_polyline_traits_2< Kernel_ > > :
  public Compute_squared_distance_2_base<CGAL::Arr_polyline_traits_2<Kernel_> >
{
public:
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_polyline_traits_2< Kernel > Traits;
  typedef Compute_squared_distance_2_base< Traits > Superclass;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Curve_2::Subcurve_const_iterator Seg_const_it;

  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename CircularKernel >
class Compute_squared_distance_2< CGAL::Arr_circular_arc_traits_2<
                                    CircularKernel > > :
  public Compute_squared_distance_2_base< CGAL::Arr_circular_arc_traits_2<
                                            CircularKernel > >
{
public: // typedefs
  typedef CircularKernel Kernel;
  typedef CGAL::Arr_circular_arc_traits_2< CircularKernel > Traits;
  typedef Compute_squared_distance_2_base< Traits > Superclass;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Point_2 Arc_point_2;
  typedef typename Kernel::Point_2 Non_arc_point_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Root_of_2 Root_of_2;

public: // methods
  double operator() ( const Non_arc_point_2& p, const X_monotone_curve_2& c );
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Compute_squared_distance_2< CGAL::Arr_conic_traits_2< RatKernel,
                                                            AlgKernel,
                                                            NtTraits > > :
  public Compute_squared_distance_2_base< CGAL::Arr_conic_traits_2< RatKernel,
                                                                    AlgKernel,
                                                                    NtTraits > >
{
public:
  typedef AlgKernel                                     Kernel;
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
  typedef Compute_squared_distance_2_base< Traits >     Superclass;
  // _Conic_point_2< AlgKernel > : public AlgKernel::Point_2
  typedef typename Traits::Point_2                      Conic_point_2;
  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public: // methods
  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Compute_squared_distance_2< CGAL::Arr_Bezier_curve_traits_2< RatKernel,
                                                            AlgKernel,
                                                            NtTraits > > :
  public Compute_squared_distance_2_base< CGAL::Arr_Bezier_curve_traits_2< RatKernel,
                                                                    AlgKernel,
                                                                    NtTraits > > {
public:
  typedef AlgKernel                                     Kernel;
  typedef CGAL::Arr_Bezier_curve_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
  typedef Compute_squared_distance_2_base< Traits >     Superclass;
  // _Conic_point_2< AlgKernel > : public AlgKernel::Point_2
  typedef typename Traits::Point_2                      Conic_point_2;
  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public: // methods
  double operator() ( const Point_2& p, const X_monotone_curve_2& c ) const;
};
template <typename Coefficient_ >
class Compute_squared_distance_2< CGAL::Arr_algebraic_segment_traits_2<
                                    Coefficient_ > > :
  public Compute_squared_distance_2_base< CGAL::Arr_algebraic_segment_traits_2<
                                            Coefficient_ > >
{
public:
  typedef Coefficient_                                  Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2<Coefficient>
                                                        Traits;
  typedef typename Traits::Bound                        FT;
  typedef typename ArrTraitsAdaptor<Traits>::Kernel     Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::CKvA_2                       CKvA_2;
  typedef std::pair< double, double >                   Coord_2;
  typedef std::vector< Coord_2 >                        Coord_vec_2;
  typedef typename X_monotone_curve_2::Curve_analysis_2 Curve;
  typedef typename Curve::Polynomial_traits_2           Polynomial_traits_2;
  typedef typename Curve::Polynomial_2                  Polynomial_2;
  typedef typename Polynomial_traits_2::Multivariate_content
                                                        Multivariate_content;
  typedef typename Polynomial_traits_2::Substitute      Substitute;
  typedef typename Polynomial_traits_2::
	  Construct_innermost_coefficient_const_iterator_range
		  ConstructInnerCoeffIter;

public:
  double operator()(const Point_2& p, const X_monotone_curve_2& c) const;

};

// TODO: Make Construct_x_monotone_subcurve_2 more generic
template <typename ArrTraits>
class Construct_x_monotone_subcurve_2
{
public:
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel  Kernel;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraits::Split_2                   Split_2;
  typedef typename ArrTraits::Intersect_2               Intersect_2;
  typedef typename ArrTraits::Multiplicity              Multiplicity;
  typedef typename ArrTraits::Construct_min_vertex_2    Construct_min_vertex_2;
  typedef typename ArrTraits::Construct_max_vertex_2    Construct_max_vertex_2;
  typedef typename ArrTraits::Compare_x_2               Compare_x_2;
  typedef typename Kernel::FT                           FT;
  typedef typename ArrTraitsAdaptor< ArrTraits >::CoordinateType
                                                        CoordinateType;
  typedef typename ArrTraits::Point_2                   Point_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  //typedef typename Kernel::Line_2 Line_2;
  //typedef typename Kernel::Compute_y_at_x_2 Compute_y_at_x_2;

  Construct_x_monotone_subcurve_2( ):
    intersect_2( this->traits.intersect_2_object( ) ),
    split_2( this->traits.split_2_object( ) ),
    compare_x_2( this->traits.compare_x_2_object( ) ),
    construct_min_vertex_2( this->traits.construct_min_vertex_2_object( ) ),
    construct_max_vertex_2( this->traits.construct_max_vertex_2_object( ) )
  { }

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.

    We assume pLeft and pRight don't lie on the curve and always do a vertical
    projection.
  */
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const Point_2& pLeft, const Point_2& pRight )
  {
    Point_2 pMin = this->construct_min_vertex_2( curve );
    Point_2 pMax = this->construct_max_vertex_2( curve );
    X_monotone_curve_2 subcurve;
    X_monotone_curve_2 unusedTrimmings;
    X_monotone_curve_2 finalSubcurve;
    if ( this->compare_x_2( pLeft, pMin ) == CGAL::LARGER )
    {
      CoordinateType y1 = this->compute_y_at_x( curve, pLeft.x( ) );
      Point_2 splitPoint( pLeft.x( ), y1 );
      this->split_2( curve, splitPoint, unusedTrimmings, subcurve );
    }
    else
    {
      subcurve = curve;
    }

    if ( this->compare_x_2( pRight, pMax ) == CGAL::SMALLER )
    {
      CoordinateType y2 = this->compute_y_at_x( subcurve, pRight.x( ) );
      Point_2 splitPoint( pRight.x( ), y2 );
      this->split_2( subcurve, splitPoint, finalSubcurve, unusedTrimmings );
    }
    else
    {
      finalSubcurve = subcurve;
    }

    return finalSubcurve;
  }

protected:
  ArrTraits traits;
  Intersect_2 intersect_2;
  Split_2 split_2;
  Compare_x_2 compare_x_2;
  Arr_compute_y_at_x_2< ArrTraits > compute_y_at_x;
  Construct_min_vertex_2 construct_min_vertex_2;
  Construct_max_vertex_2 construct_max_vertex_2;
}; // class Construct_x_monotone_subcurve_2

template <typename CircularKernel>
class Construct_x_monotone_subcurve_2< CGAL::Arr_circular_arc_traits_2<
                                         CircularKernel > >
{
public:
  typedef CGAL::Arr_circular_arc_traits_2<CircularKernel> ArrTraits;
  typedef typename ArrTraits::Intersect_2               Intersect_2;
  typedef typename ArrTraits::Split_2                   Split_2;
  typedef typename ArrTraits::Compare_x_2               Compare_x_2;
  typedef typename ArrTraits::Construct_min_vertex_2    Construct_min_vertex_2;
  typedef typename ArrTraits::Construct_max_vertex_2    Construct_max_vertex_2;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename CircularKernel::Point_2              Non_arc_point_2;
  typedef typename ArrTraits::Point_2                   Arc_point_2;
  typedef typename CircularKernel::FT                   FT;
  typedef typename CircularKernel::Root_of_2            Root_of_2;
  typedef typename CircularKernel::Root_for_circles_2_2 Root_for_circles_2_2;

public:
  Construct_x_monotone_subcurve_2( ):
    intersect_2( this->traits.intersect_2_object( ) ),
    split_2( this->traits.split_2_object( ) ),
    compare_x_2( this->traits.compare_x_2_object( ) ),
    construct_min_vertex_2( this->traits.construct_min_vertex_2_object( ) ),
    construct_max_vertex_2( this->traits.construct_max_vertex_2_object( ) )
  { }

  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const Arc_point_2& pLeft,
                                  const Arc_point_2& pRight )
  {
    Arc_point_2 pMin = this->construct_min_vertex_2( curve );
    Arc_point_2 pMax = this->construct_max_vertex_2( curve );
    X_monotone_curve_2 subcurve;
    X_monotone_curve_2 unusedTrimmings;
    X_monotone_curve_2 finalSubcurve;
    if ( this->compare_x_2( pLeft, pMin ) == CGAL::LARGER )
    {
      Arr_compute_y_at_x_2< ArrTraits > compute_y_at_x;
      FT x_approx( CGAL::to_double( pLeft.x( ) ) );
      Root_of_2 y1 = compute_y_at_x( curve, x_approx );
      Root_for_circles_2_2 intersectionPoint( x_approx, y1 );
      Arc_point_2 splitPoint( intersectionPoint );
      this->split_2( curve, splitPoint, unusedTrimmings, subcurve );
    }
    else
    {
      subcurve = curve;
    }

    if ( this->compare_x_2( pRight, pMax ) == CGAL::SMALLER )
    {
      Arr_compute_y_at_x_2< ArrTraits > compute_y_at_x;
      FT x_approx( CGAL::to_double( pRight.x( ) ) );
      Root_of_2 y2 = compute_y_at_x( subcurve, x_approx );
      Root_for_circles_2_2 intersectionPoint( x_approx, y2 );
      Arc_point_2 splitPoint( intersectionPoint );
      this->split_2( subcurve, splitPoint, finalSubcurve, unusedTrimmings );
    }
    else
    {
      finalSubcurve = subcurve;
    }

    return finalSubcurve;
  }

protected:
  ArrTraits traits;
  Intersect_2 intersect_2;
  Split_2 split_2;
  Compare_x_2 compare_x_2;
  Construct_min_vertex_2 construct_min_vertex_2;
  Construct_max_vertex_2 construct_max_vertex_2;
};

/*
 * This specialization for conic traits makes use of X_monotone_curve_2::trim,
 * which is not necessarily available.
 */
template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Construct_x_monotone_subcurve_2< CGAL::Arr_conic_traits_2< RatKernel,
                                                                 AlgKernel,
                                                                 NtTraits > >
{
public:
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > ArrTraits;
  typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename AlgKernel::Point_2 Point_2;

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.
  */
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const Point_2& pLeft, const Point_2& pRight )
  {
    // find the points on the curve
    Point_2 left = curve.point_at_x( pLeft );
    Point_2 right = curve.point_at_x( pRight );

    // make sure the points are oriented in the direction that the curve is
    // going
    AlgKernel ker;
    if (! (((curve.is_directed_right( )) &&
            ker.compare_xy_2_object() ( left, right ) == CGAL::SMALLER) ||
           ((! curve.is_directed_right( )) &&
            ker.compare_xy_2_object() ( left, right ) == CGAL::LARGER)))
    {
      Point_2 tmp = left;
      left = right;
      right = tmp;
    }

    X_monotone_curve_2 res = curve.trim( left, right );
    return res;
  }
}; // class Construct_x_monotone_subcurve_2 for Arr_conic_traits_2

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Construct_x_monotone_subcurve_2< CGAL::Arr_Bezier_curve_traits_2< RatKernel,
                                                                 AlgKernel,
                                                                 NtTraits > >
{
public:
  typedef CGAL::Arr_Bezier_curve_traits_2< RatKernel, AlgKernel, NtTraits > ArrTraits;
  typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename AlgKernel::Point_2 Point_2;

  /*
    Return the subcurve of curve bracketed by pLeft and pRight.
  */
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const Point_2& pLeft, const Point_2& pRight )
  {
    // find the points on the curve
    Point_2 left = curve.point_at_x( pLeft );
    Point_2 right = curve.point_at_x( pRight );

    // make sure the points are oriented in the direction that the curve is
    // going
    AlgKernel ker;
    if (! (((curve.is_directed_right( )) &&
            ker.compare_xy_2_object() ( left, right ) == CGAL::SMALLER) ||
           ((! curve.is_directed_right( )) &&
            ker.compare_xy_2_object() ( left, right ) == CGAL::LARGER)))
    {
      Point_2 tmp = left;
      left = right;
      right = tmp;
    }

    X_monotone_curve_2 res = curve.trim( left, right );
    return res;
  }
}; // class Construct_x_monotone_subcurve_2 for Arr_conic_traits_2

template <typename Kernel_ >
class Construct_x_monotone_subcurve_2< CGAL::Arr_linear_traits_2< Kernel_ > >
{
public: // typedefs
  typedef CGAL::Arr_linear_traits_2< Kernel_ > ArrTraits;
  typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
  typedef Kernel_ Kernel;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Ray_2 Ray_2;

public: // methods
  // curve can be unbounded. if curve is unbounded to the left,
  // pLeft is a point on the left edge of viewport.
  X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve,
                                  const Point_2& pLeft, const Point_2& pRight )
  {

    Segment_2 subsegment;

    if ( curve.is_segment( ) )
    {
      subsegment = this->constructSubsegment( curve.segment( ), pLeft, pRight );
    }
    else if (curve.is_ray( ))
    {
      Ray_2 ray = curve.ray();
      Point_2 rightP;
      Point_2 leftP;

      if (ray.source() == pRight)
      {
        rightP = pRight;
      }
      else
      {
        double right_y = arr_compute_y_at_x_2.approx(curve, CGAL::to_double(pRight.x()));
        rightP = Point_2(CGAL::to_double(pRight.x()), right_y);
      }

      if (ray.source() == pLeft)
      {
        leftP = pLeft;
      }
      else
      {
        double left_y = arr_compute_y_at_x_2.approx(curve, CGAL::to_double(pLeft.x()));
        leftP = Point_2(CGAL::to_double(pLeft.x()), left_y);
      }

      subsegment = Segment_2(leftP, rightP);
    }
    else if (curve.is_line( ))
    {
      double left_y = arr_compute_y_at_x_2.approx(curve, CGAL::to_double(pLeft.x()));
      double right_y = arr_compute_y_at_x_2.approx(curve, CGAL::to_double(pRight.x()));

      Point_2 leftP(CGAL::to_double(pLeft.x()), left_y);
      Point_2 rightP(CGAL::to_double(pRight.x()), right_y);

      subsegment = Segment_2(leftP, rightP);
    }
    else
    {
      return curve;
    }

    return X_monotone_curve_2( subsegment );
  }


protected:
  Construct_x_monotone_subcurve_2< CGAL::Arr_segment_traits_2< Kernel_ > >
    constructSubsegment;
  Arr_compute_y_at_x_2<ArrTraits> arr_compute_y_at_x_2;
};

template <typename Coefficient_ >
class Construct_x_monotone_subcurve_2< CGAL::Arr_algebraic_segment_traits_2<
                                         Coefficient_ > >
{
public: // typedefs
  typedef Coefficient_ Coefficient;
  typedef CGAL::Arr_algebraic_segment_traits_2< Coefficient > ArrTraits;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel Kernel;
  typedef typename ArrTraits::Point_2                   Point_2;
  //typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;

public: // methods
  // curve can be unbounded. if curve is unbounded to the left, pLeft is a
  // point on the left edge of viewport.
  X_monotone_curve_2 operator()(const X_monotone_curve_2& curve,
                                const Point_2& /* pLeft */,
                                const Point_2& /* pRight */)
  {
    // TODO: trim the algebraic curve
    return curve;
  }
};

// FIXME: return Traits::Point_2 instead of Kernel::Point_2
template <typename ArrTraits >
class SnapStrategy : public QGraphicsSceneMixin
{
public:
  //typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel Kernel;
  typedef typename ArrTraits::Point_2 Point_2;

  virtual Point_2 snapPoint( QGraphicsSceneMouseEvent* event ) = 0;

protected:
  SnapStrategy( QGraphicsScene* scene_ );
}; // class SnapStrategy

template <typename ArrTraits >
SnapStrategy< ArrTraits >::SnapStrategy( QGraphicsScene* scene_ )
{
  this->scene = scene_;
}

template < typename ArrTraits>
class SnapToGridStrategy : public SnapStrategy<ArrTraits>
{
public:
  typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel Kernel;
  typedef typename ArrTraits::Point_2                   Point_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;
  typedef SnapStrategy< ArrTraits >                     Superclass;

  /*! Constructors */
  SnapToGridStrategy( ) :
    Superclass( NULL ),
    gridSize( 50 )
  { }

  SnapToGridStrategy( QGraphicsScene* scene ) :
    Superclass( scene ),
    gridSize( 50 )
  { }

  /*! Destructors (virtual) */
  ~SnapToGridStrategy() {}

  Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
  {
    return this->snapPoint( event, ArrTraits( ) );
  }

  template <typename TTraits>
  Point_2 snapPoint(QGraphicsSceneMouseEvent* event, TTraits /* traits */)
  {
    QPointF clickedPoint = event->scenePos( );
    QRectF viewportRect = this->viewportRect( );
    if ( viewportRect == QRectF( ) )
    { // fallback case; we usually shouldn't end up here
      Kernel_point_2 res = this->convert( event->scenePos( ) );
      return Point_2( CGAL::to_double(res.x( )), CGAL::to_double(res.y()) );
    }

    qreal d( this->gridSize / 2.0 );
    int left = int( viewportRect.left( ) ) -
      (int( viewportRect.left( ) ) % this->gridSize);
    int right = int( viewportRect.right( ) ) +
      (this->gridSize - int( viewportRect.right( ) ) % this->gridSize);
    int x = int(clickedPoint.x( ));
    int y = int(clickedPoint.y( ));
    for ( int i = left - this->gridSize; i <= right; i += this->gridSize )
    {
      if ( i - d <= clickedPoint.x( ) && clickedPoint.x( ) <= i + d )
      {
        x = i;
        break;
      }
    }
    int top = int( viewportRect.top( ) ) -
      (int( viewportRect.top( ) ) % this->gridSize);
    int bottom = int( viewportRect.bottom( ) ) +
      (this->gridSize - int( viewportRect.bottom( ) ) % this->gridSize);
    for ( int i = top - this->gridSize; i <= bottom; i += this->gridSize )
    {
      if ( i - d <= clickedPoint.y( ) && clickedPoint.y( ) <= i + d )
      {
        y = i;
        break;
      }
    }
    //return this->convert( QPointF( x, y ) );
    Point_2 res( x, y );
    return res;
  }

  template <typename CircularKernel>
  Point_2 snapPoint(QGraphicsSceneMouseEvent* event,
                    CGAL::Arr_circular_arc_traits_2<CircularKernel>
                    /* traits */)
  {
    QPointF clickedPoint = event->scenePos( );
    QRectF viewportRect = this->viewportRect( );
    if ( viewportRect == QRectF( ) )
    {
      Kernel_point_2 res = this->convert( event->scenePos( ) );
      return Point_2( res );
    }

    qreal d( this->gridSize / 2.0 );
    int left = int( viewportRect.left( ) ) -
      (int( viewportRect.left( ) ) % this->gridSize);
    int right = int( viewportRect.right( ) ) +
      (this->gridSize - int( viewportRect.right( ) ) % this->gridSize);
    int x = int(clickedPoint.x( ));
    int y = int(clickedPoint.y( ));
    for ( int i = left - this->gridSize; i <= right; i += this->gridSize )
    {
      if ( i - d <= clickedPoint.x( ) && clickedPoint.x( ) <= i + d )
      {
        x = i;
        break;
      }
    }
    int top = int( viewportRect.top( ) ) -
      (int( viewportRect.top( ) ) % this->gridSize);
    int bottom = int( viewportRect.bottom( ) ) +
      (this->gridSize - int( viewportRect.bottom( ) ) % this->gridSize);
    for ( int i = top - this->gridSize; i <= bottom; i += this->gridSize )
    {
      if ( i - d <= clickedPoint.y( ) && clickedPoint.y( ) <= i + d )
      {
        y = i;
        break;
      }
    }
    //return this->convert( QPointF( x, y ) );
    Kernel_point_2 res( x, y );
    return Point_2( res );
  }

  void setGridSize( int size )
  {
    this->gridSize = size;
  }

protected:
  int gridSize;
  CGAL::Qt::Converter< Kernel > convert;
}; // class SnapToGridStrategy

template <typename Arr_>
class SnapToArrangementVertexStrategy:
  public SnapStrategy< typename Arr_::Geometry_traits_2 >
{
public:
  typedef Arr_                                          Arrangement;
  typedef typename Arrangement::Geometry_traits_2       Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef SnapStrategy< Traits >                        Superclass;
  typedef typename Arrangement::Vertex_iterator         Vertex_iterator;
  typedef typename Kernel::Compute_squared_distance_2
    Compute_squared_distance_2;
  typedef typename Kernel::FT                           FT;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Kernel::Point_2                      Kernel_point_2;

  SnapToArrangementVertexStrategy( ):
    Superclass( NULL ),
    arrangement( NULL )
  { }

  SnapToArrangementVertexStrategy( Arrangement* arr, QGraphicsScene* scene_ ):
    Superclass( scene_ ),
    arrangement( arr )
  { }

  Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
  {
    Kernel_point_2 clickedPoint = this->convert( event->scenePos( ) );
    return this->snapPoint( clickedPoint, Traits( ) );
  }

  template <typename TTraits>
  Point_2 snapPoint(const Kernel_point_2& clickedPoint, TTraits /* traits */)
  {
    Point_2 initialPoint( CGAL::to_double(clickedPoint.x()),
                          CGAL::to_double(clickedPoint.y()) );
    Point_2 closestPoint( CGAL::to_double(clickedPoint.x()),
                          CGAL::to_double(clickedPoint.y()) );
    bool first = true;
    FT minDist( 0 );
    QRectF viewportRect = this->viewportRect( );
    if ( viewportRect == QRectF( ) )
    {
      return initialPoint;
    }

    FT maxDist( ( viewportRect.right( ) - viewportRect.left( ) ) / 4.0 );
    for ( Vertex_iterator vit = this->arrangement->vertices_begin( );
          vit != this->arrangement->vertices_end( ); ++vit )
    {
      Point_2 point = vit->point( );
      Kernel_point_2 thisPoint( CGAL::to_double(point.x()),
                                CGAL::to_double(point.y()) );
      FT dist = this->compute_squared_distance_2( clickedPoint, thisPoint );
      if ( first || ( dist < minDist ) )
      {
        first = false;
        minDist = dist;
        closestPoint = point;
      }
    }
    if ( ! first && minDist < maxDist )
    {
      return closestPoint;
    }
    else
    {
      return initialPoint;
    }
  }

  template <typename CircularKernel>
  Point_2 snapPoint(const Kernel_point_2& clickedPoint,
                    CGAL::Arr_circular_arc_traits_2<CircularKernel>
                    /* traits */)
  {
    typedef Kernel_point_2 Non_arc_point_2;
    typedef typename CircularKernel::Circular_arc_point_2 Arc_point_2;

    Non_arc_point_2 closestKernelPoint = clickedPoint;
    Arc_point_2 closestPoint( closestKernelPoint );
    bool first = true;
    FT minDist( 0 );
    QRectF viewportRect = this->viewportRect( );
    if ( viewportRect == QRectF( ) )
    {
      return clickedPoint;
    }

    FT maxDist( ( viewportRect.right( ) - viewportRect.left( ) ) / 4.0 );
    for ( Vertex_iterator vit = this->arrangement->vertices_begin( );
          vit != this->arrangement->vertices_end( ); ++vit )
    {
      Arc_point_2 point = vit->point( );
      Non_arc_point_2 point2( CGAL::to_double(point.x( )),
                              CGAL::to_double(point.y()) );
      FT dist = this->compute_squared_distance_2( clickedPoint, point2 );
      if ( first || ( dist < minDist ) )
      {
        first = false;
        minDist = dist;
        //closestPoint = point2;
        closestPoint = point;
      }
    }
    if ( ! first && minDist < maxDist )
    {
      return closestPoint;
    }
    else
    {
      return clickedPoint;
    }
  }

  void setArrangement( Arrangement* arr )
  {
    this->arrangement = arr;
  }

protected:
  Arrangement* arrangement;
  Compute_squared_distance_2 compute_squared_distance_2;
  CGAL::Qt::Converter< Kernel > convert;
}; // class SnapToArrangementVertexStrategy

/**
   Converts between Kernel points and Arrangement points.

   The conversion is not necessarily exact.
*/
template <typename ArrTraits>
class Arr_construct_point_2
{
  typedef typename ArrTraits::Point_2                            Point_2;
  typedef typename ArrTraitsAdaptor< ArrTraits >::CoordinateType CoordinateType;
  typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel         Kernel;
  typedef typename Kernel::Point_2                               Kernel_point_2;

public:
  Point_2 operator()( const Kernel_point_2& pt );

  template <typename T >
  Point_2 operator()( const T& x, const T& y );

protected:
  template <typename T, typename TTraits >
  Point_2 operator()(const T& x, const T& y, TTraits /* traits */);

  template <typename T, typename CircularKernel >
  Point_2 operator()(const T& x, const T& y,
                     CGAL::Arr_circular_arc_traits_2<CircularKernel>
                     /* traits */);
};

class Find_nearest_edge_base : public QGraphicsSceneMixin
{
public:
  /*! Destructor (virtual) */
  virtual ~Find_nearest_edge_base() {}
};

template <typename Arr_, typename ArrTraits = typename Arr_::Geometry_traits_2>
class Find_nearest_edge : public Find_nearest_edge_base
{
public: // typedefs
  typedef Arr_ Arrangement;
  //typedef typename Arrangement::Geometry_traits_2 ArrTraits;
  typedef Compute_squared_distance_2< ArrTraits > Point_curve_distance;
  typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
  typedef CGAL::Arr_walk_along_line_point_location< Arrangement >
                                                        Point_location_strategy;
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel  Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Arrangement::Face_const_handle       Face_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Point_curve_distance::FT             FT;
  typedef typename Arrangement::Hole_const_iterator     Hole_const_iterator;
  typedef typename Arrangement::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;

public:
  /*! constructor */
  Find_nearest_edge( Arrangement* arr_ ) :
    Find_nearest_edge_base( ),
    arr( arr_ ),
    pointLocationStrategy( Point_location_strategy( *arr_ ) )
  { }

  /*! Destructor (virtual) */
  virtual ~Find_nearest_edge() {}

public: // member methods
  Halfedge_const_handle operator()( const Point_2& queryPt );

  virtual void setScene( QGraphicsScene* scene_ )
  {
    this->pointCurveDistance.setScene( scene_ );
    Find_nearest_edge_base::setScene( scene_ );
  }

protected: // member methods
  Face_const_handle getFace( const CGAL::Object& obj );

protected: // member fields
  Arrangement* arr;
  Point_curve_distance pointCurveDistance;
  Point_location_strategy pointLocationStrategy;
  Arr_construct_point_2< ArrTraits > toArrPoint;

}; // class Find_nearest_edge

template <typename Arr_, typename Coefficient_>
class Find_nearest_edge<Arr_, CGAL::Arr_algebraic_segment_traits_2<
                                Coefficient_> >: public Find_nearest_edge_base
{
public: // typedefs
  typedef Arr_ Arrangement;
  typedef typename CGAL::Arr_algebraic_segment_traits_2<
                                Coefficient_>   ArrTraits;

  typedef Compute_squared_distance_2< ArrTraits > Point_curve_distance;
  typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
  typedef CGAL::Arr_walk_along_line_point_location< Arrangement >
                                                        Point_location_strategy;
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel  Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Arrangement::Face_const_handle       Face_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Point_curve_distance::FT             FT;
  typedef typename Arrangement::Hole_const_iterator     Hole_const_iterator;
  typedef typename Arrangement::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;

public:
  /*! constructor */
  Find_nearest_edge( Arrangement* arr_ ) :
    Find_nearest_edge_base( ),
    arr( arr_ ),
    pointLocationStrategy( Point_location_strategy( *arr_ ) )
  { }

  /*! Destructor (virtual) */
  virtual ~Find_nearest_edge() {}

public: // member methods
  Halfedge_const_handle operator()( const Point_2& queryPt );

  virtual void setScene( QGraphicsScene* scene_ )
  {
    this->pointCurveDistance.setScene( scene_ );
    Find_nearest_edge_base::setScene( scene_ );
  }

protected: // member methods
  Face_const_handle getFace( const CGAL::Object& obj );

protected: // member fields
  Arrangement* arr;
  Point_curve_distance pointCurveDistance;
  Point_location_strategy pointLocationStrategy;
  Arr_construct_point_2<ArrTraits> toArrPoint;

}; // class Find_nearest_edge

template <typename Arr_, typename Kernel_ >
class Find_nearest_edge< Arr_, CGAL::Arr_linear_traits_2< Kernel_ > > :
  public Find_nearest_edge_base
{
public: // typedefs
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Geometry_traits_2       ArrTraits;
  typedef Compute_squared_distance_2<ArrTraits>         Point_curve_distance;
  typedef typename ArrTraits::X_monotone_curve_2        X_monotone_curve_2;
  typedef CGAL::Arr_walk_along_line_point_location<Arrangement>
                                                        Point_location_strategy;
  typedef typename ArrTraitsAdaptor<ArrTraits>::Kernel  Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Arrangement::Face_const_handle       Face_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Point_curve_distance::FT             FT;
  typedef typename Arrangement::Hole_const_iterator     Hole_const_iterator;
  typedef typename Arrangement::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;

public: // constructors
  Find_nearest_edge( Arrangement* arr_ ) :
    Find_nearest_edge_base( ),
    arr( arr_ ),
    pointLocationStrategy( Point_location_strategy( *arr_ ) )
  { }

public: // member methods
  Halfedge_const_handle operator()( const Point_2& queryPt );

protected: // member methods
  Face_const_handle getFace( const CGAL::Object& obj );

protected: // member fields
  Arrangement* arr;
  Point_curve_distance pointCurveDistance;
  Point_location_strategy pointLocationStrategy;
}; // class Find_nearest_edge

#endif // CGAL_ARRANGEMENTS_DEMO_UTILS_H
