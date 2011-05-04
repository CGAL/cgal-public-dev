#ifndef CGAL_QT_GENERAL_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H
#define CGAL_QT_GENERAL_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_with_holes_2.h>
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
class GeneralPolygonWithHolesGraphicsItem : public GraphicsItem
{
  typedef typename P::General_polygon_2::Traits_2 Traits;

  typedef ArrPainterOstream<Traits> PainterOstreamType;
  typedef typename PainterOstreamType::ConverterType ConverterType;

public:
  GeneralPolygonWithHolesGraphicsItem(const P* p_);

  void modelChanged();

public:
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  

  const QBrush& brush() const
  {
    return brush_;
  }

  
  void setBrush(const QBrush& b)
  {
    brush_ = b;
  }

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

  typename P::General_polygon_2::Point_2 p;
  QRectF bounding_rect;

  QBrush brush_;
  QPen vertices_pen;
  QPen edges_pen;
  bool draw_vertices;
  bool draw_edges;
};


template <typename P>
GeneralPolygonWithHolesGraphicsItem<P>::GeneralPolygonWithHolesGraphicsItem(const P * p_)
  :  poly(p_), painterostream(0),
     draw_edges(true), draw_vertices(true)   
{
  setVerticesPen(QPen(::Qt::red, 3.));
  if(poly->outer_boundary().size() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(0);
}

template <typename P>
QRectF 
GeneralPolygonWithHolesGraphicsItem<P>::boundingRect() const
{
  return bounding_rect;
}


template <typename P>
void 
GeneralPolygonWithHolesGraphicsItem<P>::paint(QPainter *painter, 
				       const QStyleOptionGraphicsItem *option,
				       QWidget * widget)
{
  ConverterType convert;

  painterostream = PainterOstreamType(painter);
  if(drawEdges()) {
	painter->setPen(this->edgesPen());

    if(!poly->is_unbounded())
    {
      typename P::General_polygon_2 boundary =  poly->outer_boundary();
      for(typename P::General_polygon_2::Curve_const_iterator eit = boundary.curves_begin();
          eit != boundary.curves_end();
          ++eit)
      {
        painterostream << *eit;
      }
    }

	for(typename P::Hole_const_iterator hit = poly->holes_begin();
      hit != poly->holes_end();
      ++hit){
		for(typename P::General_polygon_2::Curve_const_iterator eit = hit->curves_begin();
			eit != hit->curves_end();
			++eit)
		{
		  painterostream << *eit;
		}
	}
  }
/*
  if(drawVertices()) {

    painter->setPen(verticesPen());
    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
    for(typename P::General_polygon_2::Vertex_iterator it = poly->outer_boundary().vertices_begin();
        it != poly->outer_boundary().vertices_end();
        ++it){
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
GeneralPolygonWithHolesGraphicsItem<P>::updateBoundingBox()
{
  ConverterType convert;
  prepareGeometryChange();
  if(poly->is_unbounded()){
    return;
  }
  bounding_rect = convert(poly->outer_boundary().bbox());
}


template <typename P>
void 
GeneralPolygonWithHolesGraphicsItem<P>::modelChanged()
{
  if((poly->is_unbounded()) ){
    this->hide();
  } else if((!poly->is_unbounded()) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GENERAL_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H
