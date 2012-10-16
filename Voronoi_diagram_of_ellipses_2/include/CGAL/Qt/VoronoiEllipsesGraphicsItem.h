//    (c) 2008-2009 National and Kapodistrian University of Athens
//    (c) 2009-2011 INRIA Nancy
//    (c) 2011-2012 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED
    
//  (could not inherit from VoronoiEllipsesGraphicsItem, since it has
//    private members)

#ifndef CGAL_QT_VORONOI_ELLIPSES_GRAPHICS_ITEM_H
#define CGAL_QT_VORONOI_ELLIPSES_GRAPHICS_ITEM_H

#include <CGAL/Qt/EllipseBisectorGraphicsItem.h>

//#include <CGAL/Qt/GraphicsItem.h>
#include <QGraphicsItemGroup>
#include <QPen>

#include <CGAL/Qt/EllipseBisectorGraphicsItem.h>
#include <CGAL/Qt/utility.h>

//#include <omp.h>

namespace CGAL {
namespace Qt {

template <class DT>
class VoronoiEllipsesGraphicsItem : public QGraphicsItemGroup
{
  typedef typename DT::Geom_traits::Ellipse_traits   ET;

public:
//  enum { Type = UserType + 2 };

  VoronoiEllipsesGraphicsItem(const DT& dt, const QPen& pen = QPen()) {
      // workaround for iterators, as they are supported only in OpenMP 3.0
      std::vector<typename DT::Edge> elist;

      for(typename DT::Finite_edges_iterator eit = dt.finite_edges_begin();
          eit != dt.finite_edges_end(); eit++){

          elist.push_back(*eit);
      }

      int cur = 0;
      int m = elist.size();
      #pragma omp parallel for
      for(int i = 0; i < m; i++) {

          EllipseBisectorGraphicsItem<ET> *b = 
                new EllipseBisectorGraphicsItem<ET>(dt.dual(elist[i]), pen);

          #pragma omp critical(cgal_qt_bisector_drawing)
          {
            addToGroup(b);
            cur++;
            std::cerr << "Drawing bisectors... " << cur*100.0/m << "%       \r";
//            std::cerr << omp_get_thread_num()  << " : " << cur << "\n";
          }
      }
      if (m > 0) std::cerr << std::endl;
//      setZValue(5);
//      show();
  }

//  void modelChanged();
};


/*template <typename DT>
void 
VoronoiEllipsesGraphicsItem<DT>::modelChanged()
{
  update();
}*/


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_VORONOI_ELLIPSES_GRAPHICS_ITEM_H

