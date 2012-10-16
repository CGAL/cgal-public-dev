//    (c) 2009 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED
    
// TODO   (might be possible to use already existing TriangulationGraphicsItem )

#ifndef CGAL_QT_DELAUNAY_ELLIPSES_GRAPHICS_ITEM_H
#define CGAL_QT_DELAUNAY_ELLIPSES_GRAPHICS_ITEM_H

//#include <CGAL/Qt/GraphicsItem.h>
#include <QGraphicsItemGroup>
#include <QPen>

#include <CGAL/Qt/utility.h>

namespace CGAL {
namespace Qt {

template <class DT>
class DelaunayEllipsesGraphicsItem : public QGraphicsItemGroup
{
  typedef typename DT::Geom_traits::Ellipse_traits   ET;

public:
  DelaunayEllipsesGraphicsItem(const DT& dt, const QPen& pen = QPen()) {

      for(typename DT::Finite_edges_iterator eit = dt.finite_edges_begin();
          eit != dt.finite_edges_end(); eit++){

          typename DT::Geom_traits::Site_2 p1 = 
                  eit->first->vertex(eit->first->cw(eit->second))->site();
          typename DT::Geom_traits::Site_2 p2 = 
                  eit->first->vertex(eit->first->ccw(eit->second))->site();
          QGraphicsLineItem *e = 
                new QGraphicsLineItem(to_double(p1.x_center()), 
                                      to_double(p1.y_center()), 
                                      to_double(p2.x_center()), 
                                      to_double(p2.y_center()));
          e->setPen(pen);
          e->setZValue(5);
          addToGroup(e);
      }
//      setZValue(5);
//      show();
  }

//  void modelChanged();
};


/*template <typename DT>
void 
DelaunayEllipsesGraphicsItem<DT>::modelChanged()
{
  update();
}*/


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_DELAUNAY_ELLIPSES_GRAPHICS_ITEM_H

