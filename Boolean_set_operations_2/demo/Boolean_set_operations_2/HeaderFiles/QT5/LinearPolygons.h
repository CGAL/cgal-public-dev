#ifndef CGAL_QT_LINEAR_POLYGONS_H
#define CGAL_QT_LINEAR_POLYGONS_H

#include <CGAL/Qt/Converter.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <QT5/PiecewiseSetGraphicsItem.h>
#include <QT5/BoundaryPiecesGraphicsItem.h>
#include "Typedefs.h"

namespace CGAL {

namespace Qt {
struct Linear_X_monotone_bbox
{
  template<class X_monotone_linear_segment_2>
  CGAL::Bbox_2 operator()( X_monotone_linear_segment_2 const& aC ) const 
  {
    return aC.bbox();
  }
} ;
    
struct Linear_bbox
{
  template<class Linear_segment_2>
  CGAL::Bbox_2 operator()( Linear_segment_2 const& aC ) const 
  {
    double x_min = to_double(aC.source().x());
    double x_max = to_double(aC.target().x());   
    double y_min = to_double(aC.source().y()); 
    double y_max = to_double(aC.target().y());
    if(x_min > x_max)
    {
      std::swap(x_min, x_max);
      std::swap(y_min, y_max);
    }

    if(y_min > y_max)
      std::swap(y_min, y_max);

    return Bbox_2(x_min, y_min, x_max, y_max);
  }
} ;


struct Draw_linear_X_monotone_curve
{
  template<class X_monotone_linear_segment_2, class Path>
  void operator()( X_monotone_linear_segment_2 const& curve, Path& aPath, int aIdx ) const 
  {
    //typedef Simple_cartesian<double> Linear_kernel ;
    
    //commenting it gives errors
    typedef Linear_kernel::Point_2 Linear_point ;
    
    typedef CGAL::Qt::Converter<Linear_kernel> Converter ;
    
    Converter convert ;
    
    Linear_point ps( CGAL::to_double(curve.source().x()), CGAL::to_double(curve.source().y()) ) ;
    Linear_point pt( CGAL::to_double(curve.target().x()), CGAL::to_double(curve.target().y()) ) ;
      
    if( aIdx == 0 ) 
       aPath.moveTo( convert( ps ) ) ;
    else 
        aPath.lineTo( convert( ps ) ) ;

      aPath.lineTo( convert( pt ) ) ;
   }
} ;

    
struct Draw_linear_curve
{
  template<class Linear_segment_2, class Path>
  void operator()( Linear_segment_2 const& curve, Path& aPath, int aIdx ) const 
  {
    //typedef Simple_cartesian<double> Linear_kernel ;
    
    //commenting it gives errors
    typedef Linear_kernel::Point_2 Linear_point ;
      
    typedef Qt::Converter<Linear_kernel> Converter ;
    
    Converter convert ;
    
    Linear_point ps(CGAL::to_double(curve.source().x()), CGAL::to_double(curve.source().y()) ) ;
    Linear_point pt( CGAL::to_double(curve.target().x()), CGAL::to_double(curve.target().y()) ) ;
     
    if ( aIdx == 0 ) 
       aPath.moveTo( convert( ps ) ) ;
    else 
        aPath.lineTo( convert( ps ) ) ;

      aPath.lineTo( convert( pt ) ) ;
  }
} ;


template<class Linear_boundary_pieces>
class Linear_boundary_pieces_graphics_item : public Boundary_pieces_graphics_item<Linear_boundary_pieces,Draw_linear_curve,Linear_bbox>
{
  typedef Boundary_pieces_graphics_item<Linear_boundary_pieces,Draw_linear_curve,Linear_bbox> Base ;
  
public :

  Linear_boundary_pieces_graphics_item( Linear_boundary_pieces* aPieces ) : Base(aPieces) {}
} ;
    
template<class Linear_boundary>
class Linear_boundary_graphics_item : public Piecewise_boundary_graphics_item<Linear_boundary,Draw_linear_X_monotone_curve,Linear_X_monotone_bbox>
{
  typedef Piecewise_boundary_graphics_item<Linear_boundary,Draw_linear_X_monotone_curve,Linear_X_monotone_bbox> Base ;
  
public :

  Linear_boundary_graphics_item( Linear_boundary* aBoundary ) : Base(aBoundary) {}
} ;

template<class Linear_region>
class Linear_region_graphics_item : public Piecewise_region_graphics_item<Linear_region,Draw_linear_X_monotone_curve,Linear_X_monotone_bbox>
{

  typedef Piecewise_region_graphics_item<Linear_region,Draw_linear_X_monotone_curve,Linear_X_monotone_bbox> Base ;
  
public:

  Linear_region_graphics_item(Linear_region* aRegion ) : Base(aRegion) {}  
} ;
//**********************************************************************
template<class Linear_set>
class Linear_set_graphics_item : public Piecewise_set_graphics_item<Linear_set,Draw_linear_X_monotone_curve,Linear_X_monotone_bbox>
{

  typedef Piecewise_set_graphics_item<Linear_set,Draw_linear_X_monotone_curve,Linear_X_monotone_bbox> Base ;
  
public:

  Linear_set_graphics_item(Linear_set* aSet) : Base(aSet) {}
} ;

    
    
    
}//namespace Qt

}//namespace CGAL

#endif//CGAL_QT_LINEAR_POLYGONS_H
