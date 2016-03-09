// This class connects a triangulation with a scene via
// event filters on to allow manipulation
// including moving points, panning and displaying extremely large triangulations.
// This is dealt with by only allowing the vewier to view a small window of an
// over-large triangulation at one time.

/******************************************************************************/
#ifndef TRIANGULATION_VIEWER_H
#define TRIANGULATION_VIEWER_H
/******************************************************************************/

//Qt.
#include <QtGui>
#include <QMainWindow>
#include <QLineF>
#include <QRectF>
#include <QGraphicsPolygonItem>

// CGAL.
#include <CGAL/Random.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/Qt/Converter.h>

// Other.
#include <boost/format.hpp>
#include <vector>

// Local.
#include "quad_tree.h"

/******************************************************************************/

template <typename T>
class Triangulation_Viewer : public CGAL::Qt::GraphicsItem
{

    typedef typename T::Face                            Face;
    typedef typename T::Face_handle                     Face_handle;
    typedef typename T::Edge                            Edge;
    typedef typename T::Vertex                          Vertex;
    typedef typename T::Vertex_handle                   Vertex_handle;
    typedef typename T::Point                           Point;
    typedef typename T::Geom_traits                     Gt;
    typedef typename T::Segment                         Segment;
    typedef typename T::Edge_circulator                 Edge_circulator;

    CGAL::Qt::Converter<Gt>                             c;

public:

    /**************************************************************************/

    void setValid(bool v)
    {
        valid = v;
    }

    /**************************************************************************/

    void modelChanged() {

        if (!valid)
            return;

        prepareGeometryChange();

        clear_highlighted_faces();
        clear_highlighted_edges();
        clear_highlighted_vertices();

        // Free the old tree.
        if (quadtree != NULL)
            delete quadtree;

        // Rebuild the tree.
        init_structure();

        return;
    }

    /**************************************************************************/

    QRectF viewingRect() const
    {
        return view_rect;
    }

    /**************************************************************************/

    // Return the rectangle we are currently viewing.
    QRectF boundingRect() const
    {
        return bounding_rect;
    }

    /**************************************************************************/

    void paint(       QPainter*                 painter,
                const QStyleOptionGraphicsItem* option,
                      QWidget*                  widget   )
    {
        if (!valid)
            return;

        // The view rectangle.
        Point p = Point(view_rect.x(), view_rect.y());
        float w = view_rect.width();
        float h = view_rect.height();

        // Query for the correct points.
        std::list<Vertex_handle> *li;

        // Sample the current view rectangle from the quad tree.
        // TODO: improve the way the list is propogated.
        li = quadtree->sample_from_rectangle(p, w, h,  density/(w*h));

        // Draw the quadtree for debugging if needs be.
        //	quadtree->draw(scene);

        // Store the number of points we sampled.
        last_number_of_points_selected = li->size();

        QPen edge_pen      = QPen( Qt::black, 1., Qt::SolidLine );
        QPen edge_highlight = QPen( Qt::blue, 1., Qt::SolidLine );
        QPen vertex_pen    = QPen( Qt::red,   5.                );
        QPen highlight_pen = QPen( Qt::blue,  5.                );

        QMatrix matrix = painter->matrix();
        painter->resetMatrix();

        // Draw highlighted items
        painter->setBrush(QBrush(Qt::yellow, Qt::SolidPattern));
        painter->setPen( Qt::NoPen );


        typename std::vector<Face_handle>::iterator fit;
        for (fit=highlighted_faces.begin(); fit!=highlighted_faces.end(); ++fit)
        {
            QPointF points[3];
            points[0] = matrix.map( c((*fit)->vertex(0)->point()) );
            points[1] = matrix.map( c((*fit)->vertex(1)->point()) );
            points[2] = matrix.map( c((*fit)->vertex(2)->point()) );

            painter->drawConvexPolygon(points, 3);
        }

        // Draw the edges.
        painter->setPen(edge_pen);
        typename std::list<Vertex_handle>::iterator it;
        // Optimise the conditoinal outside of the loop.
        if (selection_rect == QRectF())
        {
            for (it = li->begin(); it != li->end(); ++it)
            {
                Edge_circulator eit = triangulation->incident_edges(*it),
                                        done(eit);
                if (eit!=0)
                {
                    do
                    {
                        if (!triangulation->is_infinite(*eit))
                        {
                            typename T::Segment s =triangulation->segment(*eit);
                            Point a = s.point(0);
                            Point b = s.point(1);
                            painter->drawLine( matrix.map(c(a)),
                                               matrix.map(c(b)) );
                        }
                    } while ( ++eit != done);
                }
            }
        } else {
            for (it = li->begin(); it != li->end(); ++it)
            {
                Edge_circulator eit = triangulation->incident_edges(*it),
                                          done(eit);
                if (eit!=0)
                {
                    do
                    {
                        if (!triangulation->is_infinite(*eit))
                        {
                            typename T::Segment s =triangulation->segment(*eit);
                            Point a = s.point(0);
                            Point b = s.point(1);

                            if ( selection_rect.contains(c(a))
                                     && selection_rect.contains(c(b)) )
                            {
                                painter->setPen(edge_highlight);
                                painter->drawLine( matrix.map(c(a)),
                                                   matrix.map(c(b))  );
                                painter->setPen(edge_pen);
                                continue;
                            }

                            painter->drawLine( matrix.map(c(a)),
                                               matrix.map(c(b)) );
                        }
                    } while ( ++eit != done);
                }
            }
        }

        // Draw highlighted edges.
        painter->setPen(QPen( Qt::blue, 3., Qt::SolidLine ));
        typename std::vector<Segment>::iterator eit;
        for (eit=highlighted_edges.begin(); eit!=highlighted_edges.end(); ++eit)
        {
            Point a = (*eit).point(0);
            Point b = (*eit).point(1);
            painter->drawLine(matrix.map(c(a)), matrix.map(c(b)));
        }

        // Draw the vertices.
        painter->setPen(vertex_pen);

        // If we have defined a selection rectangle.
        // We keep this conditional on the outside for efficiency.
        if (selection_rect == QRectF())
        {
            for (it = li->begin(); it != li->end(); ++it)
                painter->drawPoint(matrix.map(c((*it)->point())));

        } else {
            for (it = li->begin(); it != li->end(); ++it)
            {
                // Is the point within the selection rectangle?
                QPointF p = matrix.map(c((*it)->point()));

                if (selection_rect.contains(c((*it)->point())))
                {
                    painter->setPen(highlight_pen);
                    painter->drawPoint(p);
                    painter->setPen(vertex_pen);
                    continue;
                }

                painter->drawPoint(p);
            }
        }

    }

    /**************************************************************************/

    int signed_ceil(double v)
    {
        return ceil(abs(v)) * ((v>0) ? 1:-1);
    }

    /**************************************************************************/

    void init_structure()
    {
        if (!valid)
            return;

        // View approx. this number of points at one time.
        // TODO: Remove this magic number.
        density         = 7000;

        // Just draw all the points of small triangulations.
        if (density >= triangulation->number_of_vertices())
            density = std::numeric_limits<double>::infinity();

        // Construct a quad tree for fast pruning and sampling of the dataset.

        // Shuffle the points.
        // TODO: More elegant approach!
        std::vector<Vertex_handle> temp_list;

        // Insert all of the verticies from the triangulation into the quadtree.
        typename T::Finite_vertices_iterator v;
        for (v  = triangulation->finite_vertices_begin();
             v != triangulation->finite_vertices_end();   ++v)
         {
             temp_list.push_back(v);
         }

        std::random_shuffle( temp_list.begin(), temp_list.end() );

        quadtree = new Quad_Tree<T>(&temp_list,  1000);

        Point bl       = quadtree->get_bottom_left();
        view_rect      = QRectF( bl.x(),
                                 bl.y(),
                                 quadtree->get_w(),
                                 quadtree->get_h()  );

        bounding_rect  = QRectF(view_rect);
        selection_rect = QRectF();
        panning        = false;
    }

    /**************************************************************************/

    Triangulation_Viewer(T* triangulation)
    {
        this->valid=true;

        // Set the local pointers.
        this->triangulation = triangulation;
        this->quadtree      = NULL;

        // Initialise the data structure.
        init_structure();
    }

    /**************************************************************************/

    // This is how we intercept user events.
    // We allow the user to move the scene as usual,
    // but we each time such a movement occurs, we
    // have to redraw the view.
    bool eventFilter(QObject *obj, QEvent *event)
    {
        if (!valid)
            return false;

        // We are only interested in events sent to the graphicsview.
        QGraphicsView* v = qobject_cast<QGraphicsView*>(obj);
        if(v == NULL) {
            QWidget* viewport = qobject_cast<QWidget*>(obj);
            if(viewport == NULL) {
              return false;
            }
            v = qobject_cast<QGraphicsView*>(viewport->parent());
            if(v == NULL) {
              return false;
            }
        }

        switch (event->type())
        {
            case QEvent::MouseButtonRelease:
            {
                if (panning)
                {
                    panning = false;
                    // Don't 'gobble' release events,
                    // we want the application to accept these.
                    return false;
                }
                return false;
                break;
            }

            case QEvent::MouseButtonPress: {

                QMouseEvent* mouseEvent = static_cast<QMouseEvent*>(event);

                if( mouseEvent->button() == ::Qt::LeftButton )
                {
                    QMouseEvent* mouseEvent = static_cast<QMouseEvent*>(event);
                    start_point = v->mapToScene( mouseEvent->pos() );
                    panning     = true;
                    return false;
                }
                return false;
                break;
            }

            case QEvent::MouseMove: {
                if (panning)
                {
                    QMouseEvent* mouseEvent = static_cast<QMouseEvent*>(event);

                    // How much the scene has moved by.
                    QPointF offset = start_point -
                                         v->mapToScene( mouseEvent->pos() );

                    // Move the current center.
                    current_center = current_center + offset;
                    v->centerOn( current_center );

                    // Set a new start point.
                    start_point = v->mapToScene( mouseEvent->pos() );

                    // Update what we can see.
                    update_view_rectangle(v, false);
                    update(boundingRect() );

                    return true;
                    break;
                }
                return false;
                break;
            }

            case QEvent::Wheel: {
                QWheelEvent *wheelEvent = static_cast<QWheelEvent*>(event);
                if(wheelEvent->orientation() != ::Qt::Vertical) {
                    // This will stop the graphicsview using its built
                    // in scroll routines.
                    return true;
                }
                // The amount to zoom by.
                double zoom_ratio = 240.0;

                if( (wheelEvent->modifiers()  & ::Qt::ShiftModifier)
                  || (wheelEvent->modifiers() & ::Qt::ControlModifier) ) {
                    zoom_ratio = 120.0;
                }

                double scale_factor = pow( (double)2,
                                           wheelEvent->delta() / zoom_ratio );

                // Scale the view and recenter.
                v->scale(scale_factor, scale_factor);
                v->centerOn( current_center );

                // The view rectangle has changed.
                update_view_rectangle(v, false);
                update( boundingRect() );

                return true;
                break;
            }

            default:
                return false;
        } // switch (event->type())

        return false;
    }

    /**************************************************************************/

    void update_view_rectangle(QGraphicsView *v, bool recalculate_center)
    {
        view_rect = v->mapToScene(v->viewport()->geometry()).boundingRect();

        // The problem is, we can't get an accurate value for the center using
        // the built in Qt routines (that I have found). So when we know
        // exactly what the center should be, we set the center value ourselves.
        // Otherwise, we compute it using the following.
        if (recalculate_center)
            current_center =
                v->mapToScene(v->viewport()->rect()).boundingRect().center();
    }

    /**************************************************************************/

    void add_highlighted_face(Face_handle f)
    {
        highlighted_faces.push_back(f);
    }

    /**************************************************************************/

    void clear_highlighted_faces()
    {
        highlighted_faces.clear();
    }

    /**************************************************************************/

    void add_highlighted_vertex(Vertex_handle v)
    {
        highlighted_vertices.push_back(v);
    }

    /**************************************************************************/

    void clear_highlighted_vertices()
    {
        highlighted_vertices.clear();
    }

    /**************************************************************************/

    void add_highlighted_edge(Segment e)
    {
        highlighted_edges.push_back(e);
    }

    /**************************************************************************/

    void clear_highlighted_edges()
    {
        highlighted_edges.clear();
    }

    /**************************************************************************/

    int number_of_points_in_view()
    {
        return last_number_of_points_selected;
    }

    /**************************************************************************/

    QRectF get_view_rect()
    {
        return view_rect;
    }

    /**************************************************************************/

    double get_density()
    {
        return density;
    }

    /**************************************************************************/

    void set_density(int val)
    {
        // TODO: deal with non unit-area triangulations.
        density = (double)val;

        if (val==0)
            density = std::numeric_limits<double>::infinity();

        update( boundingRect() );
    }

    /**************************************************************************/

    bool sampling_enabled()
    {
        return density != std::numeric_limits<double>::infinity();
    }

    /**************************************************************************/

    // TODO: Currently naming conventions are a mess... fix them.
    void setSelection(QRectF r)
    {
        selection_rect = r;
    }

    /**************************************************************************/

    void clearSelection()
    {
        selection_rect = QRectF();
    }

    /**************************************************************************/
    // Return a list of vertices that are in the selection window.

    std::vector<Vertex_handle>* get_selected_vertices()
    {
        // iterate over all the vertices in the triangulation and return those
        // that are in the selection window.
        // TODO: Should this be part of a different class? perhaps
        // we should subclass the triangulation to do this?
        // Leave it here for now since it is convenient.

        std::vector<Vertex_handle>* output = new std::vector<Vertex_handle>();


        // TODO: Enable more complex selections.
        // Insert all of the verticies from the triangulation into the quadtree.
        typename T::Finite_vertices_iterator v;
        for (v  = triangulation->finite_vertices_begin();
             v != triangulation->finite_vertices_end();   ++v)
         {
             if (selection_rect.contains(c(v->point())) )
                 output->push_back(v);
         }
         return output;
    }

    /**************************************************************************/

    std::vector<Face_handle>* get_selected_faces()
    {
        std::vector<Face_handle>* output = new std::vector<Face_handle>();

        typename T::Finite_faces_iterator f;
        for (f  = triangulation->finite_faces_begin();
             f != triangulation->finite_faces_end();   ++f)
         {
             if (    selection_rect.contains(c(f->vertex(0)->point()))
                  && selection_rect.contains(c(f->vertex(1)->point()))
                  && selection_rect.contains(c(f->vertex(2)->point()))	)
                 output->push_back(f);
         }
         return output;
    }

    /**************************************************************************/

    std::vector<Edge>* get_selected_edges()
    {
        std::vector<Edge>* output = new std::vector<Edge>();

        typename T::Finite_edges_iterator e;
        for (e  = triangulation->finite_edges_begin();
             e != triangulation->finite_edges_end();   ++e)
         {
             const Segment s  = triangulation->segment(e);
             const Point & p1 = s.point(0);
             const Point & p2 = s.point(1);
             if (    selection_rect.contains(c(p1))
                  && selection_rect.contains(c(p2))	)
                 output->push_back(*e);
         }
         return output;
    }

    /**************************************************************************/

private:

    T*                         triangulation;
    Quad_Tree<T>*              quadtree;
    QPointF 			       start_point;
    bool                       panning;

    // Store the number of points we painted most recently.
    int                        last_number_of_points_selected;

    std::vector<Vertex_handle> highlighted_vertices;
    std::vector<Face_handle>   highlighted_faces;
    std::vector<Segment>       highlighted_edges;

    // Density we use when drawing triangulation.
    double                     density;

    // The rectangle that we can currently see.
    QRectF                     view_rect;

    // The rectangle bounding the whole item.
    QRectF				       bounding_rect;

    // Rectangle defining selected items.
    // We change colour of items in this rectangle.
    // to indicate user selection.
    QRectF                     selection_rect;

    // Try to keep this updated as much as possible
    // without resorting to looking at the viewport
    // transfor (which introduces inaccuracies).
    QPointF                   current_center;

    bool 			          valid;

};

/******************************************************************************/
#endif
/******************************************************************************/
