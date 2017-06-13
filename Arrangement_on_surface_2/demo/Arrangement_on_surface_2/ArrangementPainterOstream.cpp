
#include "ArrangementPainterOstream.h"

namespace CGAL {
namespace Qt {

template < typename Kernel_ >
ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> >& 
ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> > ::
operator<<( const X_monotone_curve_2& curve )
{
    std::cout<<"In operator<< Arr_segment_traits_2 X_monotone_curve_2"<<std::endl;
    const Point_2& p1 = curve.source( );
    const Point_2& p2 = curve.target( );
    Segment_2 seg( p1, p2 );

    // skip segments outside our view
    QRectF seg_bb = this->convert( seg.bbox( ) );
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.intersects( seg_bb ) 
         && (!seg.is_horizontal() && !seg.is_vertical()))
    {
      return *this;
    }

    this->painterOstream << seg;
    return *this;
}

template < typename Kernel_ >
ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> >& 
ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> > ::
operator<<( const Point_2& p )
{
    std::cout<<"In operator<< Arr_segment_traits_2 Point_2"<<std::endl;

    QPointF qpt = this->convert( p );
    // clip the point if possible
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.contains( qpt ) )
    {
      return *this;
    }

    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    radius /= this->scale;

    this->qp->drawEllipse( qpt, radius, radius );

    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );
    return *this;
}

} // namespace Qt
} // namespace CGAL
