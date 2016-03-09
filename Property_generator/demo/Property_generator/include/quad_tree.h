/******************************************************************************/
#ifndef QUAD_TREE_H
#define QUAD_TREE_H
/******************************************************************************/

#include <vector>
#include <QtGui>
#include <CGAL/Qt/Converter.h>

/******************************************************************************/
template <typename T>
class Quad_Tree
{
    typedef typename T::Face                            Face;
    typedef typename T::Face_handle                     Face_handle;
    typedef typename T::Vertex                          Vertex;
    typedef typename T::Vertex_handle                   Vertex_handle;
    typedef typename T::Point                           Point;
    typedef typename T::Geom_traits                     Gt;
    typedef typename Gt::Vector_2						Vector_2;

private:
    CGAL::Qt::Converter<Gt> c;

    float w,h;

    // Each node has its own random seed.
    double rseed;

    // Point at top left of this node.
    Point bottom_left;

    // Maximum number of elements and array to store
    // children if this is a leaf.
    int capacity;

    // Items in this node if leaf.
    // TODO: make this a NULL pointer
    // if not a leaf.
    std::vector<Vertex_handle>* items;

    // CHildren.
    Quad_Tree *children[4];


public:

    /**************************************************************************/
    // Draw the quadtree to a QGraphicsScene object.

    void draw(QGraphicsScene* scene)
    {
        if (children[0]!=NULL)
        {
            QLineF line1( c( bottom_left+Vector_2( w/2, 0   )),
                          c( bottom_left+Vector_2( w/2, h   )));
            QLineF line2( c( bottom_left+Vector_2( 0,   h/2 )),
                          c( bottom_left+Vector_2( w,   h/2 )));

            scene->addItem(new QGraphicsLineItem(line1));
            scene->addItem(new QGraphicsLineItem(line2));

            for (int i=0; i<4; i++)
            {
                children[i]->draw(scene);
            }
        }
    }

    /**************************************************************************/

    // When freeing a quadtree, delete all children first.
    ~Quad_Tree()
    {
        if (children[0] != NULL)
        {
            for (int i=0; i<4; i++)
                delete children[i];
        }
    }

    /**************************************************************************/

    Quad_Tree(Point p, float w, float h, int capacity)
    {
        this->items       = new std::vector<Vertex_handle>();
        this->bottom_left = p;
        this->w           = w ;
        this->h           = h;
        this->capacity    = capacity;
        this->rseed       = (double)rand()/(double)RAND_MAX;

        this->items->reserve(capacity);

        // This is the condition for being a leaf.
        children[0]    = NULL;
    }

    /**************************************************************************/

    Quad_Tree(std::vector<Vertex_handle> *vs, int capacity)
    {
        double min_x = std::numeric_limits<double>::infinity();
        double min_y = std::numeric_limits<double>::infinity();
        double max_x = 0;
        double max_y = 0;

        typename std::vector<Vertex_handle>::const_iterator i;
        for (i=vs->begin(); i!= vs->end(); ++i)
        {
            Point p = (*i)->point();

            if (p.x() < min_x)
                min_x = p.x();
            if (p.y() < min_y)
                min_y = p.y();
            if (p.x() > max_x)
                max_x = p.x();
            if (p.y() > max_y)
                max_y = p.y();
        }

        // There is a bug in Qt... If you set the bounding box for an item to
        // infinity, then it doesn't matter how much you update the bounding box
        // after that, the scene will not draw correctly. So if there are no
        // points we use zeros.
        if (min_x == std::numeric_limits<double>::infinity())
            min_x = 0;

        if (min_y == std::numeric_limits<double>::infinity())
            min_y = 0;

        this->items       = new std::vector<Vertex_handle>();
        this->bottom_left = Point(min_x, min_y);
        this->w           = max_x - min_x ;
        this->h           = max_y - min_y;
        this->capacity    = capacity;
        this->rseed       = (double)rand()/(double)RAND_MAX;

        this->items->reserve(capacity);

        // This is the condition for being a leaf.
        children[0]    = NULL;

        insert(vs);
    }

    /**************************************************************************/

    void insert(std::vector<Vertex_handle> *vs)
    {
        typename std::vector<Vertex_handle>::const_iterator i;

        for (i=vs->begin(); i != vs->end(); ++i)
            this->insert(*i);
    }

    /**************************************************************************/

    double get_w()
    {
        return w;
    }

    /**************************************************************************/

    double get_h()
    {
        return h;
    }

    /**************************************************************************/

    Point get_center()
    {
        return Point(bottom_left->x() + w/2, bottom_left->x() + h/2 );
    }

    /**************************************************************************/

    Point get_bottom_left()
    {
        return Point(bottom_left);

    }

    /**************************************************************************/

    // Add a point to the tree.
    void insert(Vertex_handle v)
    {
        Point p = v->point();

        // If this is a leaf.
        if (children[0] == NULL)
        {
            if (items->size() + 1 > capacity)
            {
                // Create the child nodes.
                children[0] = new Quad_Tree(bottom_left + Vector_2(0  , h/2),
                                            w/2, h/2, capacity );
                children[1] = new Quad_Tree(bottom_left + Vector_2(w/2, h/2),
                                            w/2, h/2, capacity );
                children[2] = new Quad_Tree(bottom_left,
                                            w/2, h/2, capacity );
                children[3] = new Quad_Tree(bottom_left + Vector_2(w/2, 0  ),
                                            w/2, h/2, capacity );

                // Insert the points back into the tree lower down.
                for (int i=0; i< items->size(); i++)
                    this->insert((*items)[i]);

	                // Free the old item list.
                delete items;

                this->insert(v);

            } else  {
                // Just add the point to this leaf.
                items->push_back(v);
            }
        } else {
            // Decide which quadrant to insert this point into.
            int  EW  = value_in_range(p.x(), bottom_left.x(),
                                             bottom_left.x() + w/2)  ?  0:1;

            int  NS  = value_in_range(p.y(), bottom_left.y(),
                                             bottom_left.y() + h/2)  ? 1:0;
            children[ (2*NS) | EW ]->insert(v);
        }
    }

    /**************************************************************************/

    bool value_in_range(float x, float a, float b)
    {
        return (a <= x && x < b);
    }

    /**************************************************************************/

    // Get all the points within the specified rectangle. If the sample density
    // is too great, take a random sample to approximately the right density.
    // Notice use of list for constant time merging.
    std::list<Vertex_handle> *sample_from_rectangle( Point  r_bottom_left,
                                                float  r_w,
                                                float  r_h,
                                                float density           )
    {
        std::list<Vertex_handle> *results = new std::list<Vertex_handle>();

        // Does this node intersect the rectangle? If yes, and we are a leaf
        // node, return points up to a certain density. Else continue down
        // the tree.

        bool C1 = value_in_range( r_bottom_left.x(),   bottom_left.x(),
                                                       bottom_left.x() +  w );

        bool C2 = value_in_range(   bottom_left.x(), r_bottom_left.x(),
                                                     r_bottom_left.x() + r_w );

        bool C3 = value_in_range( r_bottom_left.y(),   bottom_left.y(),
                                                       bottom_left.y() +  h );

        bool C4 = value_in_range(   bottom_left.y(), r_bottom_left.y(),
                                                     r_bottom_left.y() + r_h );

        // If they intersect.
        if ( (C1 || C2) && (C3 || C4) )
        {

            if (children[0] != NULL)
            {
                // Recurse down tree
                for (int i=0; i<4; i++)
                {
                    std::list<Vertex_handle> *l;
                    l = children[i]->sample_from_rectangle( r_bottom_left,
                                                            r_w,
                                                            r_h,
                                                            density        );

                    // Append the list of points from the child to this list.
                    results->splice(results->end(), *l);
                }
            } else {
                // Sample up to density.
                if ((items->size() != 0) &&  (capacity / (w*h) >= density))
                {
                    double val = w*h*density;
                    if (val < 1)
                    {
                        if (rseed < val)
                            results->push_back((*items)[0]);
                    } else {
                        int max = floor(w*h*density);
                        if (max >= items->size())
                            max = items->size();
                        for (int i=0; i<max; i++)
                            results->push_back((*items)[i]);
                    }
                } else {
                    results->insert(results->end(), items->begin(),
                                                    items->end()    );
                }
            }
        }
        return results;
    }
};

/*****************************************************************************/
#endif
/*****************************************************************************/
