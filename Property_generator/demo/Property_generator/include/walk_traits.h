/******************************************************************************/
#ifndef WALK_STATISTICS_H
#define WALK_STATISTICS_H
/******************************************************************************/

#include "walk.h"

namespace CGAL {

/******************************************************************************/

template <typename T>
class Walk_2_face_property_traits
{
public:

    /********************************************************************/

    typedef typename T::Face_handle          Face_handle;
    typedef typename T::Segment				 Segment;
    typedef typename T::Edge				 Edge;
    typedef typename T::Point                Point;

    /********************************************************************/
    // Tedious hacking of enums with string representations.
    /********************************************************************/

    typedef typename T::Face_handle                           Entity_type;
    typedef Walk<T>                                           Container_type;
    typedef typename std::vector<Face_handle>::const_iterator Iterator;

    enum  Property
    {
		CIRCUM_RADIUS,
		CIRCLE_POWER,
		EUCLIDEAN_DISTANCE,
		FACE_INDEX
    };

protected:
    static int get_number_of_properties(){ return 4;}
	
    /********************************************************************/
    
    static std::string property_to_string(Property s)
    {
        switch(s)
        {
            case CIRCUM_RADIUS:      return "Circum Radius";
            case CIRCLE_POWER:  	 return "Circle Power";
            case EUCLIDEAN_DISTANCE: return "Euclidean Distance";
            case FACE_INDEX:         return "Face Index";
        }
    }

    /********************************************************************/
    // Iterators

    static Iterator begin( const Walk<T> &w )
    {
        return w.begin();
    }

    static Iterator end( const Walk<T> &w )
    {
        return w.end();
    }

    /********************************************************************/

    // Return a statistic over one value.
    static 
    double 
    compute_property (       
        const Walk<T>&     w,
        const Face_handle& f,
              Property     s  
    )
    {
        switch(s)
        {
            case CIRCUM_RADIUS:
            {
                // TODO: Deal with infinities.
                if (w.is_infinite(f)) return 0;

                float r2;
                const Point & v = f->vertex(0)->point();

                // The circum radius is the distance from the circum
                // center to any vertex.
                r2 = ( w.circumcenter(f) - v ).squared_length();

                return  sqrt(r2);
            }

            case CIRCLE_POWER:
            {
                if (w.is_infinite(f)) return 0;

                // Find the circumcenter.
                const Point & c = w.circumcenter(f);

                // This is the squared distance between the circumcenter
                // and the destination point.
                double d2 = (c - w.destination()).squared_length();

                // This is the squared radius of the circumball.
                double r2 = (c-(f)->vertex(0)->point()).squared_length();

                // The circle power is d^2 - r^2.
                return d2 - r2;
            }

            case EUCLIDEAN_DISTANCE:
            {
                if (w.is_infinite(f)) return 0;

                // The end point of the walk.
                const Point & p = w.destination();
                return sqrt((w.circumcenter(f)-p).squared_length());
            }

            case FACE_INDEX:
            {
                return w.get_face_index(f);
            }
        }
        return NULL;
    }	

    /**************************************************************************/
	
};

}

/******************************************************************************/
#endif
/******************************************************************************/
