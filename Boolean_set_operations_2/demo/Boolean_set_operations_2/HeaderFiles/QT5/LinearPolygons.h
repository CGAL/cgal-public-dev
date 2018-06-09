#ifndef CGAL_QT_LINEAR_POLYGONS_H
#define CGAL_QT_LINEAR_POLYGONS_H

#include <CGAL/Qt/Converter.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <QT5/PiecewiseSetGraphicsItem.h>
#include <QT5/BoundaryPiecesGraphicsItem.h>
#include "Typedefs.h"

namespace CGAL {

<<<<<<< HEAD:Boolean_set_operations_2/demo/Boolean_set_operations_2/HeaderFiles/QT5/LinearPolygons.h
namespace Qt {
=======
template<class Linear_polygons_>    
class Iterator_and_polygons
{
public:    
    typedef Linear_traits::Curve_const_iterator Curve_const_iterator;
    typedef Linear_polygons_                    Linear_polygons;
  typedef typename Linear_polygons_::Base Base;
} ;
    
>>>>>>> af6cb4b3ce9eab9443661804548debe57e5c5f44:Boolean_set_operations_2/demo/Boolean_set_operations_2/include/QT5/Linear_polygons.h
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
/*//take care it might be circular
    if(aC.is_circular())
    {
      typedef typename Circle_segment_2::Circle_2 Circle_2 ;

      const Circle_2& circ = aC.supporting_circle();

      return circ.bbox();
    }
  */  
    return Bbox_2(x_min, y_min, x_max, y_max);
  }
} ;


struct Draw_linear_X_monotone_curve
{
  template<class X_monotone_linear_segment_2, class Path>
  void operator()( X_monotone_linear_segment_2 const& curve, Path& aPath, int aIdx ) const 
  {
    typedef Simple_cartesian<double> Linear_kernel ;
    
    typedef Point_2<Linear_kernel> Linear_point ;
    
    typedef CGAL::Qt::Converter<Linear_kernel> Converter ;
    
    Converter convert ;
    
    Linear_point lS=new Linear_point( CGAL::to_double(curve.source().x()), CGAL::to_double(curve.source().y()) ) ;
    Linear_point lT( CGAL::to_double(curve.target().x()), CGAL::to_double(curve.target().y()) ) ;
      
    if( aIdx == 0 ) 
       aPath.moveTo( convert( lS ) ) ;
    else 
        aPath.lineTo( convert( lS ) ) ;

      aPath.lineTo( convert( lT ) ) ;
   }
} ;

    
struct Draw_linear_curve
{
  template<class Linear_segment_2, class Path>
  void operator()( Linear_segment_2 const& curve, Path& aPath, int aIdx ) const 
  {
    typedef Simple_cartesian<double> Linear_kernel ;
    
    typedef Point_2<Linear_kernel> Linear_point ;
      
    typedef Qt::Converter<Linear_kernel> Converter ;
    
    Converter convert ;
    
    Linear_point lS=new Linear_point( CGAL::to_double(curve.source().x()), CGAL::to_double(curve.source().y()) ) ;
    Linear_point lT( CGAL::to_double(curve.target().x()), CGAL::to_double(curve.target().y()) ) ;
     
    if ( aIdx == 0 ) 
       aPath.moveTo( convert( lS ) ) ;
    else 
        aPath.lineTo( convert( lS ) ) ;

      aPath.lineTo( convert( lT ) ) ;
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