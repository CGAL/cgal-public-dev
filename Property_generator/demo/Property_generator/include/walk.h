/*******************************************************************************
* Written by Ross Hemsley for INRIA.fr.
* This is a class heirachy to perform different walks on triangulations,
* providing methods to create a QGraphicsItem so that they can be drawn
* directly to a QGraphicsScene.
*
* We also consider the number of triangles and properties of each walk.
*
*******************************************************************************/
#ifndef WALK_H
#define WALK_H
/******************************************************************************/

// Qt.
#include <QtGui>
#include <QMainWindow>
#include <QLineF>
#include <QRectF>
#include <QGraphicsPolygonItem>

// CGAL.
#include <CGAL/Random.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/point_generators_2.h>

// Other.
#include <boost/format.hpp>
#include <vector>

/******************************************************************************/

// Use this to stop the walks keeping track of the faces and pivots,
// to gain more accurate timing information.
// #define DOING_TIMING

/*******************************************************************************
* Abstract class to contain different walking strategies
*******************************************************************************/
// TODO: Refactor this: it needed be a class heirarchy, just
// add a paramter for choosing the walk type.

template <typename T>
class Walk
{
public:
    typedef typename T::Face                            Face;
    typedef typename T::Face_handle                     Face_handle;
    typedef typename T::Line_face_circulator            Lfc;
    typedef typename T::Geom_traits                     Gt;
    typedef typename T::Point                           Point;

public:
    long                         index_count;

    // Base class used for initialisation purposes.
                                 Walk(T* dt);
    Face_handle 				 do_walk( Point p0,
                                          Face_handle f = Face_handle() );

    void 					     set_dt( T* dt );
    void 						 clear();

    // Create a graphics item for drawing this triangulation.
    QGraphicsItemGroup*          getGraphics( QPen       pen   = QPen(),
                                              QBrush     brush = QBrush() ) const;
    int                          getNumTrianglesVisited()      const;
    int                          getNumOrientationsPerformed() const; 
	
	int                          get_face_index(Face_handle f) const;

    // This provides an interface to properties of the faces.
    std::vector<Face_handle>*    get_faces();
    Point                        circumcenter( const Face_handle f ) const;
    bool                         is_infinite ( const Face_handle f ) const;
    Point                        destination()                       const;

    // Static helper functions to draw 2D faces to QQrahpicsItems.
    static QGraphicsPolygonItem* drawTriangle( Face_handle f,
                                               QPen        pen   = QPen(),
                                               QBrush      brush = QBrush() );

    static QGraphicsItemGroup*   drawArrow( const QPoint p,
                                            const QPoint q,
                                            int   size,
                                            QPen  pen = QPen() );

	typename std::vector<Face_handle>::const_iterator begin() const;	
    typename std::vector<Face_handle>::const_iterator end()   const;
	

protected:
    // Pointer to the triangluation this walk is on.
    T*                           dt;
    // The destination point and start face.
    Point						 p0;
    Face_handle					 f0;

    // This allows subclasses to add faces to the current walk.
    // Doing this enables the base-class functions to work.
    void                         addToWalk( Face_handle f );

    CGAL::Orientation            orientation( Point p, Point q, Point r );

private:
    // List of faces this walk intersects.
    std::vector<Face_handle>          faces;
    int 						 o_count;

};

/*******************************************************************************
* Straight walk strategy
*******************************************************************************/

template <typename T>
class StraightWalk : public Walk<T>
{
public:
    typedef typename T::Face                            Face;
    typedef typename T::Point                           Point;
    typedef typename T::Face_handle                     Face_handle;
    typedef typename T::Line_face_circulator            Lfc;
    typedef typename T::Geom_traits                     Gt;

private:
    using Walk<T>::addToWalk;

public:

    StraightWalk(T* dt) : Walk<T>(dt)
    {
    }

    Face_handle do_walk(Point p,Face_handle f=Face_handle())
    {
        this->clear();

        this->p0 = p;
        this->f0 = f;

        // Create a circulator describing the walk.
        CGAL::Qt::Converter<Gt> c;

        if (f==Face_handle())
            f=this->dt->infinite_face();

        Point x = f->vertex(0)->point();

        // Use CGAL built-in line walk.
        Lfc lfc = this->dt->line_walk (x,p), done(lfc);

        // Take all the items from the circulator and add them to a list.
        if (lfc != 0)
        {
            do {
                Face_handle f = lfc;
                addToWalk(f);
            } while (++lfc != done);
        }

        return lfc;
    }
};

/*******************************************************************************
* Pivot Walk strategy
*******************************************************************************/

template <typename T>
class PivotWalk : public Walk<T>
{
public:
    typedef typename T::Face                            Face;
    typedef typename T::Point                           Point;
    typedef typename T::Face_handle                     Face_handle;
    typedef typename T::Vertex_handle                   Vertex_handle;
    typedef typename T::Geom_traits                     Gt;

private:

    // We store the pivot points for this walk so we can draw them later.
    #ifndef DOING_TIMING
        std::vector<Point> pivots;
    #endif

    // Make sure we use the right orientation predicate.
    using Walk<T>::orientation;
    using Walk<T>::addToWalk;

public:

    /**************************************************************************/

    PivotWalk(T* dt) : Walk<T>(dt) {}

    /**************************************************************************/

    PivotWalk( Point        p,
               T*           dt,
               Face_handle  f      = Face_handle(),
               long         rseed  = time(NULL)                 ) : Walk<T>(dt)
    {
        do_walk(p,f,rseed);
    }

    /**************************************************************************/

    Face_handle do_walk(  Point       p,
                          Face_handle f      = Face_handle(),
                          long        rseed  = time(NULL)      )
    {

        // Count the number of index tests... TODO: remove this debugging
        // information.
        this->index_count = 0;

        // Clear the walk status.
        #ifndef DOING_TIMING
            pivots.clear();
        #endif
        this->clear();
        this->p0 = p;
        this->f0 = f;

        // We use these to switch the direction without using comparisons.
        CGAL::Orientation left_turn;
        int cw_offset, ccw_offset;

        // Create a binary random number generator.
        CGAL::Random random(rseed);

        // Initialise the first face properly.
        if(f == Face_handle())
        {
            // Index of a non-infinite face.
            this->index_count++;
            // TODO: make this neater.
            int i=this->dt->infinite_face()->index(this->dt->infinite_vertex());
            f    =this->dt->infinite_face()->neighbor(i);
        } else if(this->dt->is_infinite(f)){
            this->index_count++;
            f = f->neighbor(f->index(this->dt->infinite_vertex()));
        }

        // This is where we store the current face.
        Face_handle c    = f;
        Face_handle prev = c;
        addToWalk(c);

        // **     FIND FIRST FACE      ** //
        bool found = false;
        for (int i=0; i<3; i++)
        {
            const Point & p0 = c->vertex(i)->point();
            const Point & p1 = c->vertex(c->cw(i))->point();

            // If we have found a face that can see the point.
            if ( orientation(p0,p1,p) == CGAL::POSITIVE )
            {
                found = true;
                c     = c->neighbor(c->ccw(i));
                addToWalk(c);
                break;
            }
        }
        // Point was in the first face.
        if (!found)
            return c;
        // ** END OF FIND FIRST FACE ** //

        // Loop over the pivots.
        while(1)
        {
            // If we have walked out of the hull, we are done.
            if (this->dt->is_infinite(c))
                return c;

            // Choose which direction we are going to go around this pivot
            // at random.
            bool clockwise = random.get_bool();

            if (clockwise)
            {
                cw_offset  = 2;
                ccw_offset = 1;
                left_turn  = CGAL::LEFT_TURN;
            } else {
                cw_offset  = 1;
                ccw_offset = 2;
                left_turn  = CGAL::RIGHT_TURN;
            }

            // 'i' is always the index of the previously visited triangle
            // relative to the current one.
            this->index_count++;
            int i = c->index(prev);

            // This is the pivot we are going to use, it is opposite
            // the last face.
            const Point & p_pivot = c->vertex(i)->point();

            // Maintain a list of pivots for drawing and statistics.

            #ifndef DOING_TIMING
                pivots.push_back(p_pivot);
            #endif

            // Skip forwards by two orientations in accordance with the
            // 'skipping' strategy.
            prev = c;
            c    = c->neighbor((i+ccw_offset)%3);
            addToWalk(c);
            this->index_count++;
            i    = c->index(prev);
            prev = c;
            c    = c->neighbor((i+cw_offset)%3);
            addToWalk(c);

            // If we end up in an infinite cell, we have to keep walking
            // until we  are back inside the convex hull, otherwise
            // the orientations might be invalid.
            while (this->dt->is_infinite(c))
            {
                this->index_count++;
                i    = c->index(prev);
                prev = c;
                c    = c->neighbor((i+cw_offset)%3);
                addToWalk(c);
            }

            // Re-calculate index of the previous face relative to current.
            this->index_count++;
            i    = c->index(prev);

            // Now decide which direction we are going to go based on the
            // result of the following orientation.
            // To do the test, we use the pivot, and the remaining point
            // which is the point opposite from the previous face.
            const Point & p1 = c->vertex(i)->point();

            // If it's a left turn then we need to start going back the other
            // way, otherwise we can continue going in the same direction.
            if (orientation(p_pivot, p1, p) != left_turn)
            {

                // This is an idea for an optimisation that removes a few
                // index tests. It uses the fact that we already
                // know some of the values in some cases.
                if (orientation( p_pivot,
                                 c->vertex((i+cw_offset)%3)->point(),
                                 p ) != left_turn)
                {
                    Face_handle tmp = c;
                    c    = prev;
                    prev = tmp;

                    // Reverse direction without comparisons.
                    // TODO: probably shouldn't assume this is valid code.
                    clockwise  = !clockwise;
                    cw_offset  = cw_offset  *2 % 3;
                    ccw_offset = ccw_offset *2 % 3;
                    left_turn  = CGAL::Orientation((int)left_turn*-1);

                } else{
                    // Check the point isn't contained in this face:
                    if (orientation( p1,
                                     c->vertex((i+cw_offset)%3)->point(),
                                     p) == left_turn)
                    {
                        // New pivot.
                        prev = c;
                        c    = c->neighbor((i+ccw_offset)%3); // (right)
                        addToWalk(c);
                        continue;
                    }
                    else
                    {
                        // Point was located.
                        return c;
                    }
                }

                // Previous code before optimisation.
                /*
                prev = c->neighbor((i+cw_offset)%3);

                clockwise = !clockwise;
                if (clockwise)
                {
                    cw_offset  = 2;
                    ccw_offset = 1;
                    left_turn  = CGAL::LEFT_TURN;
                } else {
                    cw_offset  = 1;
                    ccw_offset = 2;
                    left_turn  = CGAL::RIGHT_TURN;
                }
                */

            } else {
                prev = c;
                c = c->neighbor((i+cw_offset)%3);
                addToWalk(c);
            }

            // We now have a direction, a previous, and a current triangle.
            // Start walking in this direction until we find the sink.
            // We only need one orientation per triangle we walk through at
            // this point.
            while(1)
            {
                // Check that we haven't walked outside of the convex hull.
                if (this->dt->is_infinite(c))
                    return c;

                this->index_count++;
                int i = c->index(prev);

                // The orientation for the next 'spoke' around this pivot.
                const Point & p1 = c->vertex(i)->point();
                if (orientation(p_pivot, p1, p) != left_turn)
                {
                    // If the test 'failed' then we have either arrived, or
                    // found the sink. Do an extra test to check to see if we
                    // are at the sink.
                    int   p2_index   = (i+cw_offset)%3;
                    const Point & p2 = c->vertex( p2_index )->point();
                    if ( orientation(p1,p2,p) == left_turn)
                    {
                        // NEW PIVOT //

                        // Go through this edge to the next pivot.
                        prev = c;
                        c    = c->neighbor((i+ccw_offset)%3);
                        addToWalk(c);
                        break;
                    } else {
                        // POINT LOCATED //
                        return c;
                    }
                } else {
                    // Continue going in the same direction.
                    prev = c;
                    c    = c->neighbor((i+cw_offset)%3);
                    addToWalk(c);
                }
            }
        }

        // Error condition.
        return Face_handle();
    }

    /**************************************************************************/

    int getNumPivots()
    {
        #ifndef DOING_TIMING
            return pivots.size();
        #else
            return 0;
        #endif
    }

    /**************************************************************************/

    // Overload get graphics to draw the pivots.
    QGraphicsItemGroup* getGraphics( QPen pen=QPen(), QBrush brush=QBrush() ) const
    {

        CGAL::Qt::Converter<Gt> c;

        // Invoke the base-class drawing method to get
        // the triangles involved.
        QGraphicsItemGroup* g = Walk<T>::getGraphics(pen,brush);

        // The drawing style for the pivots.
        QPen   e_pen(Qt::blue);
        QBrush e_brush(Qt::blue);

        #ifndef DOING_TIMING

        // Iterate over pivots in this walk.
        typename std::vector<typename T::Point>::const_iterator i;
        for (i = pivots.begin(); i != pivots.end(); ++i)
        {

            QGraphicsEllipseItem *e;
            e = new QGraphicsEllipseItem(QRectF(QPointF(-4,-4), QSizeF(8,8)));

            // Don't resize this point, but do move it!
            e->setFlag(QGraphicsItem::ItemIgnoresTransformations);
            e->setPos(c(*i));
            e->setBrush(e_brush);
            e->setPen(e_pen);
            g->addToGroup(e);

        }
        #endif

        return g;
    }

    /**************************************************************************/

};

/*******************************************************************************
* Visibility walk strategy
*******************************************************************************/

template <typename T>
class VisibilityWalk : public Walk<T>
{
public:
	
    typedef typename T::Face                            Face;
    typedef typename T::Point                           Point;
    typedef typename T::Face_handle                     Face_handle;
    typedef typename T::Geom_traits                     Gt;

private:
    using Walk<T>::orientation;
    using Walk<T>::addToWalk;

public:
    /**************************************************************************/

    VisibilityWalk(T* dt) : Walk<T>(dt) {}

    /**************************************************************************/

    VisibilityWalk(Point           p,
                   T*              dt,
                   Face_handle     f      = Face_handle(),
                   long            rseed  = time(NULL)          ) : Walk<T>(dt)
    {
        do_walk(p,f,rseed);
    }

    /**************************************************************************/

    Face_handle do_walk( Point       p,
                         Face_handle f=Face_handle(),
                         long        rseed=time(NULL)  )
    {
        // Initialise the base class.
        this->clear();
        this->p0 = p;
        this->f0 = f;
        this->index_count=0;


        // The user did not provide a face handle. So just use the infinite
        // face.

        if(f == Face_handle()){
            // Index of a non-infinite face.
            this->index_count++;
            int i=this->dt->infinite_face()->index(this->dt->infinite_vertex());
            f    =this->dt->infinite_face()->neighbor(i);
        } else if(this->dt->is_infinite(f)){
            this->index_count++;

            f = f->neighbor(f->index(this->dt->infinite_vertex()));
        }

        // This is where we store the current face.
        Face_handle    c = f;
        Face_handle prev = c;

        addToWalk(c);

        // Create a binary random number generator.
        CGAL::Random random(rseed);

        // **     FIND FIRST FACE      ** //
        bool found = false;
        for (int i=0; i<3; i++)
        {
            const Point & p0 = c->vertex(i       )->point();
            const Point & p1 = c->vertex(c->cw(i))->point();

            // If we have found a face that can see the point.
            if ( orientation(p0,p1,p) == CGAL::POSITIVE )
            {
                found=true;
                c = c->neighbor(c->ccw(i));
                break;
            }
        }
        // The point could have been in the first face.
        if (!found)
            return c;
        // ** END OF FIND FIRST FACE ** //

        // Loop until we find our destination point.
        while(1) //for (int i=0; i<10000; i++)
        {
            addToWalk(c);

            // This will stop us getting stuck in loops.
            if (this->dt->is_infinite(c))
                return c;

            this->index_count++;
            int i = c->index(prev);

            const Point & p0 = c->vertex( i                )->point();
            const Point & p1 = c->vertex( this->dt->cw(i)  )->point();
            const Point & p2 = c->vertex( this->dt->ccw(i) )->point();

            int left_first   = random.get_bool();

            // We randomise the order in which we test the
            // faces we are walking through
            if (left_first)
            {

                if ( orientation(p0,p1,p) == CGAL::POSITIVE ) {
                    prev = c;
                    c = c->neighbor( this->dt->ccw(i) );
                    continue;
                }

                if ( orientation(p2,p0,p) == CGAL::POSITIVE ) {
                    prev = c;
                    c = c->neighbor( this->dt->cw(i) );
                    continue;
                }

            } else {

                if ( orientation(p2,p0,p) == CGAL::POSITIVE ) {
                    prev = c;
                    c = c->neighbor( this->dt->cw(i) );
                    continue;
                }

                if ( orientation(p0,p1,p) == CGAL::POSITIVE ) {
                    prev = c;
                    c = c->neighbor( this->dt->ccw(i) );
                    continue;
                }
            }

            // We must have located the point.
            return c;
        }
    }

    /**************************************************************************/
};

/*******************************************************************************
* Walk base-class functions
*
* These provide basic functionality that is common to all of the walk types
* so that we do not have to re-implement this functionality multiple times.
*
*******************************************************************************/

// Overload the orientation predicate to enable statistics gathering.
template <typename T>
CGAL::Orientation Walk<T>::orientation(Point p, Point q, Point r)
{
    o_count++;
    return CGAL::orientation(p,q,r);
}

/******************************************************************************/

// Add a face to the walk. We gather statistics about the walk here too.
template <typename T>
inline void Walk<T>::addToWalk(Face_handle f)
{
    // Don't keep track of faces if we are doing timing.
    #ifndef DOING_TIMING
        faces.push_back(f);
    #endif
}

/******************************************************************************/

// Return the number of faces visited in this walk.
template <typename T>
int Walk<T>::getNumTrianglesVisited() const
{
    return faces.size();
}

/******************************************************************************/

// Create a graphics item representing this walk.
template <typename T>
QGraphicsItemGroup* Walk<T>::getGraphics( QPen pen, QBrush brush ) const
{
    // This GraphicsItem Group will store the triangles from the walk.
    QGraphicsItemGroup* g = new QGraphicsItemGroup();

    CGAL::Qt::Converter<Gt> c;

    // If we are drawing as a graph or a set of faces.
    if (FALSE)
    {
        // Get the circumcenters
        typename std::vector<typename T::Face_handle>::const_iterator f;

        Point last;

        bool first=true;
        for (f=faces.begin(); f != faces.end(); ++f)
        {
            if (dt->is_infinite(*f)) continue;

            Point center  = centroid( (*f)->vertex(0)->point(),
                                      (*f)->vertex(1)->point(),
                                      (*f)->vertex(2)->point()  );

            if (!first)
                g->addToGroup(drawArrow(c(last).toPoint(), c(center).toPoint(),
                                                 7, QPen(QColor("#FF0000"),6)));

            last = center;

            first=false;
        }


    } else {
        // Iterate over faces in this walk.
        typename std::vector<typename T::Face_handle>::const_iterator i;
        for (i = faces.begin(); i != faces.end(); ++i)
        {
            // Draw this triangle in the walk.
            if (! dt->is_infinite( *i ) )
            {
                QGraphicsPolygonItem *tr = drawTriangle(*i,pen,brush);
                g->addToGroup(tr);
            }
        }
    }

    return g;
}

/******************************************************************************/

template <typename T>
Walk<T>::Walk(T *dt)
{
    this->dt = dt;
    o_count   = 0;
    index_count =0;
}

/******************************************************************************/

template <typename T>
void Walk<T>::clear()
{
    o_count   = 0;
    #ifndef DOING_TIMING
        faces.clear();
    #endif
}

/******************************************************************************/

template <typename T>
int Walk<T>::getNumOrientationsPerformed() const
{
    return o_count;
}

/******************************************************************************/

template <typename T>
typename T::Point Walk<T>::destination() const
{
    return p0;
}

template <typename T>
std::vector<typename T::Face_handle>* Walk<T>::get_faces()
{
    return &faces;
}

/******************************************************************************/

template <typename T>
typename std::vector<typename T::Face_handle>::const_iterator Walk<T>::begin() const
{
	return faces.begin();
}	

/******************************************************************************/

template <typename T>
typename std::vector<typename T::Face_handle>::const_iterator Walk<T>::end() const
{
	return faces.end();
}	

/******************************************************************************/

template <typename T>
typename T::Point Walk<T>::circumcenter( Face_handle f) const
{
    return dt->circumcenter(f);
}

/******************************************************************************/

template <typename T>
bool Walk<T>::is_infinite( Face_handle f) const
{
    return dt->is_infinite(f);
}

/******************************************************************************/

template <typename T>
void Walk<T>::set_dt(T* dt)
{
    this->dt = dt;
    clear();
}

/******************************************************************************/

// Helper-function to create a triangle graphics item.
// Note that this is publically accessible and static.
template <typename T>
QGraphicsPolygonItem* Walk<T>::drawTriangle( Face_handle f,
                                             QPen        pen,
                                             QBrush      brush )
{
    // Helper to convert between different point types.
    CGAL::Qt::Converter<Gt> c;

    // We store this triangle as a polygonItem.
    QGraphicsPolygonItem *polygonItem = 0;

    // Convert a face into a polygon for plotting.
    QVector<QPointF>      polygon;

    polygon << c(f->vertex(0)->point())
            << c(f->vertex(1)->point())
            << c(f->vertex(2)->point());

    polygonItem = new QGraphicsPolygonItem(QPolygonF(polygon));

    // The "look" of the triangle.
    polygonItem->setPen(pen);
    polygonItem->setBrush(brush);

    return polygonItem;
}

/******************************************************************************/

template <typename T>
int Walk<T>::get_face_index( Face_handle f ) const
{
    return std::find(faces.begin(), faces.end(), f) - faces.begin();	
}

/******************************************************************************/

template <typename T>
QGraphicsItemGroup* Walk<T>::drawArrow( const QPoint p,
                                        const QPoint q,
                                              int    size,
                                              QPen    pen   )
{
    // The line of the arrow.
    QLineF line(p, q);

    int half_size = size/2;

    // The center point of the line.
    QPointF c  = line.pointAt(0.5);

    // Unit vectors running normal and parallel to the line.
    QLineF  u = line.normalVector().unitVector();
    QLineF  v = line.unitVector();


    // Point on line which is at front of arrow head.
    QPointF c1 = c + half_size * QPointF( v.dx(), v.dy() );
    // Point on the line which is at the back of the arrow head.
    QPointF c2 = c - half_size * QPointF( v.dx(), v.dy() );

    // Bottom and top of the arrow head.
    QPointF arrow_bottom = c2 +  size * QPointF(u.dx(), u.dy());
    QPointF arrow_top    = c2 -  size * QPointF(u.dx(), u.dy());


    // Create the line items for the arrow head.
    QLineF a1(c1, arrow_bottom);
    QLineF a2(c1, arrow_top   );


    // Create the graphics item for this arrow.
    QGraphicsItemGroup* g = new QGraphicsItemGroup();

    QGraphicsLineItem* l1 = new QGraphicsLineItem(line);
    QGraphicsLineItem* l2 = new QGraphicsLineItem(a1);
    QGraphicsLineItem* l3 = new QGraphicsLineItem(a2);

    l1->setPen(pen);
    l2->setPen(pen);
    l3->setPen(pen);

    g->addToGroup(l1);
    g->addToGroup(l2);
    g->addToGroup(l3);

    return g;
}

/******************************************************************************/

#endif

/******************************************************************************/
