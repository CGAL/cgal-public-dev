/******************************************************************************
* Written by Ross Hemsley for INRIA.fr.
* This is a class heirachy to perform different walks on triangulations,
* providing methods to create a QGraphicsItem so that they can be drawn
* directly to a QGraphicsScene.
*
*
******************************************************************************/
#ifndef VERTEX_WALK_H
#define VERTEX_WALK_H
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

template <typename T>
class Vertex_Walk
{
    typedef typename T::Vertex_circulator                   Vertex_circulator;
    typedef typename T::Vertex_handle                       Vertex_handle;
    typedef typename T::Face                                Face;
    typedef typename T::Face_handle                         Face_handle;
    typedef typename T::Geom_traits::Vector_2               Vector_2;
    typedef typename T::Line_face_circulator                Lfc;
    typedef typename T::Geom_traits                         Gt;
    typedef typename T::Point                               Point;
    // Just to cut down on line length later on.
    typedef typename std::pair<Vertex_handle,Vertex_handle> Vertex_pair;
    typedef typename std::vector<Vertex_pair>               Vertex_pair_vector;


public:
    double distance_per_step;
    double radius;
    double progress;
    double angle;
    double step_length;
    double path_length;
    int    substeps;
    int    path_steps;
    int    neighbour_count;
    int    step_count;

public:

    /**************************************************************************/

    enum Strategy {STANDARD, RESTRICTED, MONOTONE};

    /**************************************************************************/

    Vertex_Walk(T* dt)
    {
        this->dt = dt;
    }


private:
    /**************************************************************************/
    // We consider the line (p1,p3), and then grow the smallest circle
    // that has p0 on its border and whose diamter completely intersects
    // the line.        
    //                               
    //             , - ~ ~ ~ - ,
    //         , '               ' ,  
    //       ,                       o   p2  
    //      ,                        |,      
    //     ,                         | ,                    
    //     o-------------+-----------|-,------------------ o  p3
    //     , p1                        ,
    //      ,                         ,
    //       ,                       ,
    //         ,                  , '
    //           ' - , _ _ _ ,  '
    //
    // 
    //
    
    double calculuate_disc_radius(Point p1, Point p2, Point p3)
    {

        Vector_2  v1;
        Vector_2  v2;
            
        // The angle p3,p1,p2
        v1 = (p3-p1);
        v2 = (p2-p1);
        v1 = v1 / std::sqrt(v1*v1);
        v2 = v2 / std::sqrt(v2*v2);

        // If the angle is >= a right angle,
        // we would need an infinite radius.
        if (std::acos( v1*v2 ) > boost::math::constants::pi<double>()/2 )
            return std::numeric_limits<double>::infinity();

        // The distance from p2 to p1.
        double x = sqrt((p2-p1).squared_length());

        return x/(2*(v1*v2));
    }

    /**************************************************************************/
    // Get the point that we used to add a given neighbour into the
    // neighbours list.

    Vertex_handle get_predecessor(
        Vertex_handle        v,
        Vertex_pair_vector&  l
        )
    {
        typename Vertex_pair_vector::const_iterator i;
        
        for (i=l.begin(); i!= l.end(); ++i)
            if ((*i).first == v)
                return (*i).second;
        
        // We should have found the point, if we didn't then
        // something has broken.
        std::cerr << "Something strange happened..." << std::endl;
        return NULL;
    }

    /**************************************************************************/
    // Clear the given list and then add all the neighbours of u to it.

    // Store a pair containing the vertex that added each neighbour
    // so that we can find a path in the Delaunay triangulation
    // at the end of each step.
    void add_neighbours(
        Vertex_pair_vector&  l,
        Vertex_handle        u, 
        bool                 keep
        )
    {
        // Empty the list to start.
        if (!keep)
            l.clear();
        
        int x =0;

        Vertex_circulator v = dt->incident_vertices(u), done(u);
        if (v != 0)
        {            
            do
            {
                if (x>1000) return;
                ++x;
                if (dt->is_infinite(v))
                {
                    l.clear();
                    return;
                }
                
                // Is the piont already in the list?
                bool contained=false;
                typename Vertex_pair_vector::const_iterator i;
                for (i= l.begin(); i!= l.end(); ++i)
                {
                    if ((*i).first == v)
                    {   
                        contained=true;
                        break;
                    }
                }
                
                // If the point wasn't already in the list, add it,
                // storing the point whose neighbour it was
                // so we can efficiently construct a path later.
                if(! contained )
                    l.push_back(std::make_pair(v,u));

            } while ( ++v != done );
        }
    }

    /**************************************************************************/
    // The restricted vertex walk walking strategy.

    Face_handle restricted_walk()
    {
        CGAL::Qt::Converter<Gt> c;

        // Keep track of these quantities for testing.
        step_count        = 0;
        distance_per_step = 0;
        radius            = 0;
        progress          = 0;
        angle             = 0;
        step_length       = 0;
        path_length       = 0;
        substeps          = 0;
        path_steps        = 0;
        neighbour_count   = 0;

        Point A = f0->vertex(0)->point();
        steps.push_back(A);

        // Use this to store the current neighbours at each step.
        Vertex_pair_vector neighbours;
        add_neighbours(neighbours, f0->vertex(0), false);
                
        // Current radius stores the last radius during a substep
        // so that we can ensure it always grows.        
        double current_radius = 0;    
        
        // Unit vector describing the direction.
        Vector_2 direction = (p0-A);
        direction          = direction/sqrt(direction.squared_length());
        Point last_point   = A;

        while(true)
        {            
            double min_radius = std::numeric_limits<double>::infinity();
            Vertex_handle min_radius_vertex;
            
            typename Vertex_pair_vector::const_iterator i;
            for (i = neighbours.begin(); i != neighbours.end(); ++i)
            {
                Vector_2 u = ((*i).first->point()-A);
                
                // Calculate this radius.
                double this_radius = u.squared_length()/(2*u*direction);
                                
                if (    this_radius > current_radius 
                     && this_radius < min_radius     )
                {
                    min_radius        = this_radius;
                    min_radius_vertex = (*i).first;
                }                
            }

            // Termination condition.
            if (sqrt((A-p0).squared_length()) < 2*min_radius)
            {
                vertices.push_back(min_radius_vertex);
                return Face_handle();
            }
                    
            Vector_2 u = min_radius_vertex->point()-A;

            // Is the point in the decision cone?
            // Avoid inverse cos using a little re-arranging.
            if (u*direction/sqrt(u.squared_length()) >= sqrt( 2+sqrt(2) )/2 )
            {
                ////////////////////////////////////////////////////
                // Compute the path in the Delaunay triangulation //

                // This is the vertex that we found inside the current search
                // cone.
                Vertex_handle z = min_radius_vertex;
                Vertex_handle w;

                std::list<Vertex_handle> path;

                // Keep going until we end up back at A.
                while (z->point() != A)
                {
                    // Go back to the point that was a neighbour of 
                    // the current point and rempeat.

                    // Keep track of the number of path steps that we do.
                    // note that this is different to the number of substeps.

                    w = get_predecessor(z, neighbours);
                    path_length += 
                        sqrt((z->point()-w->point()).squared_length());
                    z=w;
                    
                    // Count the number of points, not the number
                    // of edges.
                    if (z->point() != A)
                        path_steps++;

                    path.push_back(z);
                }

                // Keep track of the whole path for displaying later.
                path.reverse();
                typename std::list<Vertex_handle>::const_iterator j;
                for (j=path.begin(); j!=path.end(); ++j)
                    vertices.push_back(*j);

                ////////////////////////////////////////////////////

                step_count++;
                neighbour_count += neighbours.size();
                progress +=
                        sqrt( (A-p0).squared_length()) -
                        sqrt( (min_radius_vertex->point()-p0).squared_length());

                radius += calculuate_disc_radius(
                              A,
                              min_radius_vertex->point(),
                              p0);

                step_length += sqrt(
                    (last_point - min_radius_vertex->point()).squared_length());

                steps.push_back(min_radius_vertex->point());
                current_radius = 0;
                A              = min_radius_vertex->point();
                last_point     = A;
                direction      = (p0-A);
                direction      = direction/sqrt(direction.squared_length());
                add_neighbours(neighbours, min_radius_vertex, false);
            } else {
                // If the closest point is not in the decision cone,
                // do a new substep.
                substeps++;
                step_length += sqrt(
                    (last_point - min_radius_vertex->point()).squared_length());
                add_neighbours(neighbours, min_radius_vertex, true);
                current_radius = min_radius;
                last_point     = min_radius_vertex->point();
            }
        }
        
        std::cout << "Didn't terminate correctly" << std::endl;
        return Face_handle();        
    }

    /**************************************************************************/
    // The monotone walk walking strategy.

    Face_handle monotone_walk()
    {
        CGAL::Qt::Converter<Gt> c;
        
        // Keep track of these quantities for testing.
               step_count        = 0;
        double distance_per_step = 0;
        double radius            = 0;
        double progress          = 0;
        double angle             = 0;
        int    substeps          = 0;
        double step_length       = 0;
        int    neighbour_count   = 0;

        // The list of vertices that we have yet to visit,
        // and that we are currently visiting.
        std::vector<std::pair<Vertex_handle,Vertex_handle> > neighbours;

        // Arbitrary start point given start face.
        Point A = f0->vertex(0)->point();

        // This is a unit vector describing our direction.
        Vector_2 direction = this->p0 - A;
        direction = direction/sqrt(direction.squared_length());
        
        // For drawring purposes.
        steps.push_back(A);
        vertices.push_back(f0->vertex(0));

        add_neighbours(neighbours, f0->vertex(0), false);
                
        double current_radius = 0;    
                
        // Kill the walk when this grows larger than say... 30.        
        int substep_count     = 0;        
        int max_substep_count = 30;        
        double substep_length = 0;
        
        Point last_point = A;
        
        // This will gradually grow until we are in the decision cone.
        while(true) //for (int k=0; k<100; k++)
        {
            double min_radius = std::numeric_limits<double>::infinity();
            Vertex_handle min_radius_vertex;
            
            typename Vertex_pair_vector::const_iterator i;
            for (i = neighbours.begin(); i != neighbours.end(); ++i)
            {
                Vector_2 u = ((*i).first->point()-A);
                
                // Calculate this radius.
                double this_radius = u.squared_length()/(2*u*direction);
                                
                if (this_radius > current_radius && this_radius < min_radius)
                {
                    min_radius        = this_radius;
                    min_radius_vertex = (*i).first;
                }                
            }
            
            // Should never happen, keep for now for debugging.
            if (min_radius == std::numeric_limits<double>::infinity())
            {
                qDebug() << "Something weird happened "  << endl;                
                break;
            }
                    
            Vector_2 u = min_radius_vertex->point()-A;
            double alpha = std::acos(u*direction/sqrt(u.squared_length())); 
                        
            if( alpha <= boost::math::constants::pi<double>()/8 )
            {
                ////////////////////////////////////////////////////
                // Compute the path in the Delaunay triangulation // 
                
                // This is the vertex that we found inside the current search
                // cone. 
                Vertex_handle z = min_radius_vertex;
                
                std::list<Vertex_handle> path;
                
                // Keep going until we end up back at A.
                while (z->point() != A)
                {
                    // Go back to the point that was a neighbour of 
                    // the current point and rempeat. 
                    
                    z = get_predecessor(z, neighbours);
                    path.push_back(z);
                }
                
                typename std::list<Vertex_handle>::const_iterator j;
                
                path.reverse();
                
                for (j=path.begin(); j!=path.end(); ++j)
                    vertices.push_back(*j);

                ////////////////////////////////////////////////////
                
                substep_length += 
                sqrt((last_point-min_radius_vertex->point()).squared_length());

                
                // TODO: How are we going to make sure this is unique?
                neighbour_count += neighbours.size();
                
                substeps        += substep_count;
                step_length     += substep_length;
                
                // If the closest point is in the decision cone, do a new step.
                // Statistics gathering for an individual step //
                step_count ++;
                distance_per_step += sqrt(u.squared_length());
                radius            += u.squared_length()/(2*u*direction);
                progress          += u*direction;
                angle             += alpha;
                // 
                
               // vertices.push_back(min_radius_vertex);
                steps.push_back(min_radius_vertex->point());                                            
                substep_count  = 0;
                substep_length = 0;
                current_radius = 0;
                A              = min_radius_vertex->point();            
                last_point     = A;
                add_neighbours(neighbours, min_radius_vertex, false);
            } else {
                // If the closest point is not in the decision cone, do a new substep.

                substep_length += sqrt( ( last_point - min_radius_vertex->point() ).squared_length() );
                last_point      = min_radius_vertex->point();
                substep_count++;
               // vertices.push_back(min_radius_vertex);
                add_neighbours(neighbours, min_radius_vertex, true);
                current_radius = min_radius;
            }
            
            // Just to deal with border effects: TODO: make less messy.
            if( substep_count > max_substep_count )    
                break;
            
            // If we hit the infinite vertex.
            if (neighbours.size() == 0)
                break;    
        }            

        // average_step_length       = step_length / step_count;
        // average_neighbour_count   = neighbour_count   /(double)step_count;
        // average_distance_per_step = distance_per_step / step_count;
        // average_x_progress        = progress / step_count;
        // average_angle             = angle    / step_count;
        // average_radius            = radius   / step_count;
        // average_substeps          = substeps / (double)step_count;
        
        return Face_handle();
    }

    /**************************************************************************/

    Face_handle standard_walk()
    {

        // Choose any point on the current face to start the walk from.
        Vertex_handle current    = f0->vertex(0);
        double current_distance  = (current->point()-p0).squared_length();
        
        vertices.push_back(current);
        
        // Continue the walk until we reach the destination.
        while(true)
        {
        
            Vertex_circulator v = dt->incident_vertices(current), done(v);
            double min          = std::numeric_limits<double>::infinity();
            Vertex_handle best  = Vertex_handle();
        
            // Loop through all of the neighbouring vertices
            // until we find the closest neighbour.
            // TODO: Try to visit fewer vertices?
            if (v != 0)
            {
                do
                {
                    double distance = (v->point() - p0).squared_length();
                    if (distance < min)
                    {
                        min  = distance;
                        best = v;
                    }
                } while ( ++v != done );
            }
        
            // We are at the closest vertex.
            if (min >= current_distance)
                break;
        
            // Continue walking.
            current          = best;
            current_distance = min;
            vertices.push_back(current);
        }
            
        
        // Place-holder for the end of the walk.
                return Face_handle();
        
    }

    /**************************************************************************/
    // ATTENTION: IN PROGRESS.
    // TODO: Refactor.

public:
    Face_handle do_walk( Point p0, Face_handle f = Face_handle(),
                         Strategy s = STANDARD)
    {

        if (f== Face_handle())
            f=dt->infinite_face();

        this->clear();

        this->p0       = p0;
        this->f0       = f;

//        s=RESTRICTED;

        this->strategy = s;

        switch(s)
        {
            case RESTRICTED: return restricted_walk();
            case MONOTONE:   return monotone_walk();
            case STANDARD:   return standard_walk();
        }
    }

    /**************************************************************************/

    void set_dt( T* dt )
    {
        this->dt = dt;
    }

    /**************************************************************************/

    void clear()
    {
        vertices.clear();
        steps.clear();
    }

    /**************************************************************************/
    // Draw a  circle of diameter between two points.

    QGraphicsItemGroup *draw_circle(Point p1, Point p2)
    {

        QPen pen(Qt::red,0);


        QGraphicsItemGroup* g = new QGraphicsItemGroup();

        CGAL::Qt::Converter<Gt> c;

        Point center = Point( (p1.x() + p2.x())/2, (p1.y() + p2.y())/2);

        double rad=0.5*sqrt((p1-p2).squared_length());

        assert(rad>0);


        QGraphicsEllipseItem *e;
        e = new QGraphicsEllipseItem( QRectF( QPointF(-rad,-rad),
                                              QSizeF(2*rad,2*rad) ) );

        e->setPen(pen);
        e->setPos(c(center));

        // This is the point where the cone lines finish
        // relative to the circle.
        Point p_middle =  p1 + (p2-p1)*(0.5+sqrt(2)/4);

        typename T::Geom_traits::Vector_2 v =
                    (p2-p1).perpendicular(CGAL::RIGHT_TURN);

        v = v/sqrt(v.squared_length());

        g->addToGroup(e);

        Point p = p_middle + v*rad*1/sqrt(2);
        Point q = p_middle - v*rad*1/sqrt(2);

        QLineF line1(c(p1), c(p));
        QLineF line2(c(p1), c(q));

        QGraphicsLineItem *l1 = new QGraphicsLineItem(line1);
        QGraphicsLineItem *l2 = new QGraphicsLineItem(line2);

        l1->setPen(pen);
        l2->setPen(pen);

        g->addToGroup(l1);
        g->addToGroup(l2);

        return g;
    }

    /**************************************************************************/

    // Create a graphics item for drawing this triangulation.
    QGraphicsItemGroup*          getGraphics( QPen       pen   = QPen(),
                                              QBrush     brush = QBrush() )
    {
        CGAL::Qt::Converter<Gt> c;
        QPen   e_pen(Qt::blue,0);
        QBrush e_brush(Qt::blue);

        QGraphicsItemGroup* g = new QGraphicsItemGroup();


        Point last_point;

        typename std::vector<Point>::iterator s;

        bool first=true;
        for (s = steps.begin(); s != steps.end(); ++s)
        {
            if (!first)
            {
                if (strategy != MONOTONE)
                {
                    // Calculuate radius of circle/cone.
                    double  d = calculuate_disc_radius(last_point, *s, p0);
                    Point tmp = last_point +
                     (p0-last_point)/sqrt((p0-last_point).squared_length())*2*d;
                    g->addToGroup(draw_circle(last_point, tmp));
                } else {

                    // TODO: FIX THIS BUG PROPERLY.
                    if (f0 == Face_handle())
                        continue;

                    typename T::Geom_traits::Vector_2 direction =
                               (p0-f0->vertex(0)->point());

                    double d = calculuate_disc_radius( last_point,
                                                       *s,
                                                       last_point + direction );
                    Point tmp = last_point +
                              (direction)/sqrt(direction.squared_length())*2*d;

                    g->addToGroup(draw_circle(last_point, tmp));
                }
            }

            first = false;
            last_point = *s;
        }

        typename std::vector<Vertex_handle>::iterator i;

        bool is_first =true;
        for (i = vertices.begin(); i != vertices.end(); ++i)
        {

            if (!is_first)
            {
                QLineF line(c(last_point), c((*i)->point()));
                QGraphicsLineItem* l = new QGraphicsLineItem(line);
                l->setPen(e_pen);
                g->addToGroup(l);

            } else {
                is_first = false;
            }

            QGraphicsEllipseItem *e;
            e = new QGraphicsEllipseItem(QRectF(QPointF(-4,-4), QSizeF(8,8)));

            // Don't resize this point, but do move it!
            e->setFlag(QGraphicsItem::ItemIgnoresTransformations);
            e->setPos(c((*i)->point()));
            e->setBrush(e_brush);
            e->setPen(e_pen);
            g->addToGroup(e);

            last_point = (*i)->point();
        }
        return g;
    }

    /**************************************************************************/

    Point destination()
    {
        return p0;
    }

    /**************************************************************************/

private:
    // Pointer to the triangluation this walk is on.
    T*                           dt;
    // The destination point and start face.
    Point                         p0;
    Face_handle                     f0;

    Strategy strategy;

    // List of faces this walk intersects.
    std::vector<Vertex_handle>   vertices;

    // Temporarily store the vertices
    // that define full steps for restricted
    // vertex walk.
    std::vector<Point>           steps;

};

/******************************************************************************/
#endif
/******************************************************************************/
