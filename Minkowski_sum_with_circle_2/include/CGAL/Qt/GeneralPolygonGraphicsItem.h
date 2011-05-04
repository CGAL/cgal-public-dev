#ifndef CGAL_QT_GENERAL_POLYGON_GRAPHICS_ITEM_H
#define CGAL_QT_GENERAL_POLYGON_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
//#include <CGAL/Circular_arc_2.h>
//#include <CGAL/Segment_2.h>
#include <CGAL/apply_to_range.h>
#include <CGAL/Qt/ArrPainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename P>
class GeneralPolygonGraphicsItem : public GraphicsItem
{
  typedef typename P::Traits_2 Traits;

  typedef ArrPainterOstream<Traits> PainterOstreamType;
  typedef typename PainterOstreamType::ConverterType ConverterType;

//  typedef typename CGAL::Segment_2<Traits> Segment_2;
//  typedef typename CGAL::Circular_arc_2<Traits> Circular_arc_2;
public:
  GeneralPolygonGraphicsItem(const P* p_);

  void modelChanged();

public:
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  

  const QPen& verticesPen() const
  {
    return vertices_pen;
  }

  const QPen& edgesPen() const
  {
    return edges_pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    vertices_pen = pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    edges_pen = pen;
  }

  bool drawVertices() const
  {
    return draw_vertices;
  }

  void setDrawVertices(const bool b)
  {
    draw_vertices = b;
    update();
  }

  bool drawEdges() const
  {
    return draw_edges;
  }

  void setDrawEdges(const bool b)
  {
    draw_edges = b;
    update();
  }

protected:
  void updateBoundingBox();

  const P * const poly;
  QPainter* m_painter;
  PainterOstreamType painterostream;

  typename P::Point_2 p;
  QRectF bounding_rect;

  QPen vertices_pen;
  QPen edges_pen;
  bool draw_edges;
  bool draw_vertices;
};


template <typename P>
GeneralPolygonGraphicsItem<P>::GeneralPolygonGraphicsItem(const P * p_)
  :  poly(p_), painterostream(0),
     draw_edges(true), draw_vertices(true)   
{
  setVerticesPen(QPen(::Qt::red, 3.));
  if(poly->size() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(3);
}

template <typename P>
QRectF 
GeneralPolygonGraphicsItem<P>::boundingRect() const
{
  return bounding_rect;
}




template <typename P>
void 
GeneralPolygonGraphicsItem<P>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * widget)
{
  painter->setPen(this->edgesPen());
  painterostream = PainterOstreamType(painter);
  if(drawEdges()) {
    for(typename P::Curve_const_iterator eit = poly->curves_begin();
        eit != poly->curves_end();
        ++eit)
    {
  	  painterostream << *eit;
    }
  }


 /* if(drawVertices()) {
    ConverterType convert;

    painter->setPen(verticesPen());
    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
    for(typename P::Vertex_iterator it = poly->vertices_begin();
        it != poly->vertices_end();
        it++){
      QPointF point = matrix.map(convert(*it));
      painter->drawPoint(point);
    }
  }
  */
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename P>
void 
GeneralPolygonGraphicsItem<P>::updateBoundingBox()
{
  ConverterType convert;
  prepareGeometryChange();
  if(poly->size() == 0){
    return;
  }
  bounding_rect = convert(poly->bbox());
}


template <typename P>
void 
GeneralPolygonGraphicsItem<P>::modelChanged()
{
  if((poly->size() == 0) ){
    this->hide();
  } else if((poly->size() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GENERAL_POLYGON_GRAPHICS_ITEM_H
