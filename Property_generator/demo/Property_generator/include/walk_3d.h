#ifndef WALK_3D_H
#define WALK_3D_H

// This implements pivot walk and visibility walk on 3 dimensional 
// triangulations. Will probably be refactored into the walk class later.



using namespace std;

/*****************************************************************************/

// Template referencing the triangulation data structure
template <typename T>
class Walk_3d
{
    
    typedef typename T::Cell                            Cell;
    typedef typename T::Cell_handle                     Cell_handle;
    typedef typename T::Geom_traits                     Gt;        
    typedef typename T::Point                           Point;                
    
public:
          Walk_3d(T* dt);
    int do_walk(Point        p0,
                 Cell_handle* f=Cell_handle());
    
    int          getNumOrientationsPerformed();
    int          getNumTrianglesVisited();
    
    // return a pointer to the faces.
    vector<Cell_handle> *getFaces();
    
    
protected:
                
    bool same_side_as_cell(Point p, Cell_handle cell, int index);
    int  get_neighbor(bool cw, int pivot_index, int last_face);
    void addToWalk(Cell_handle f);
    void clear();

    
    CGAL::Orientation orientation(Point p, Point q, Point r, Point s);
    T* dt;    
            
private:
    
    // Reference to the triangulation.
    int o_count;
    vector<Cell_handle>              faces;
    
};

/*****************************************************************************/

template <typename T>
class PivotWalk_3d : public Walk_3d<T>
{
    
    typedef typename T::Cell                            Cell;
    typedef typename T::Cell_handle                     Cell_handle;
    typedef typename T::Geom_traits                     Gt;        
    typedef typename T::Point                           Point;                    
    typedef typename T::Vertex_handle                   Vertex_handle;            
    
    
public:
    PivotWalk_3d(T* dt) : Walk_3d<T>(dt) {}
    
    int do_walk(Point p0, Cell_handle f=Cell_handle())
    {                
        
        this->clear();
        
        CGAL::Random random(time(NULL));
                
        //cout << "Doing pivot walk" << endl;
            
        if (f==Cell_handle())
            f=this->dt->infinite_cell();
    
        Cell_handle c     = f;
        Cell_handle prev;
                
        // Find the starting cell that we are going to walk into.
        // We will then define then 3-pivot to be the point opposite
        // the facet we walked through into this cell.
        // This will define a topological sphere which we can walk
        // around using the 2D pivot walking strategy.
        // (Or even the 2D visibility walking strategy!)
        
        // ** FIND FIRST CELL ** //
        bool found=false;
        for (int x=0; x<4; x++)
        {
            // If this is not the facet we have just walked through
            if (! this->same_side_as_cell(p0, c, x) )
            {
                found = true;
                prev  = c;
                c     = c->neighbor( x );        
                this->addToWalk(c);
                break;
            }
        }        
        // Almost never happens.
        if (!found)
        {
            //cout << "Query point seems to be in first cell" << endl;
            // The first cell contained the point.
            if ( this->dt->locate(p0) == c )
            {
                return 1;
            } else{
                cerr << "TEST FAILED" << endl;
                return -1;
            }                                    
        }
        // ** END OF FIND FIRST CELL ** //

        for (int k=0; k<100000; k++)
        {    
            //cout << "New 3-pivot" << endl;
            
            // The 3-pivot we will use will have the same index
            // as the facet we just walked through.
            int i = c->index(prev);
        
            Vertex_handle pivot_3 = c->vertex(i);
        
            // *** Given the current 3-pivot, we now do a 2D walk *** //    
                            
            // Find the starting face for this lower dimensional walk 
            // (note that we are now effectively starting a 2D walk from scratch).
            
            bool found=false;
                        
            // Start from a random point
            int offset = random.get_int(0,4);

            for (int j=0; j<4; j++ )
            {
                            
                int x = (j+offset)%4;
            
                // We ignore the face which is opposite the pivot.
                if (x == i) continue;            
                if (! this->same_side_as_cell(p0, c, x) )
                {
                    // We've found a direction we can go in.
                    //cout << "Found starting cell for 2D walk." << endl;
                    prev = c;
                    c    = c->neighbor(x);
                    this->addToWalk(c);                    
                    found=true;
                    break;
                }        
            }
            
            if (!found)
            {
                // The point was in the first cell of this walk. 
                if ( this->dt->locate(p0) == c )
                {
                    //cout << "Found query point at start of 2D walk." << endl;
                    return 1;
                } else{
                    cerr << "TEST FAILED" << endl;
                    return -1;
                //    exit(1);
                }                                                    
            }
        
            // Now, walk through the faces until we reach the sink of
            // the star for this 3-pivot.
        
            bool done=false;
            // TODO: unbound the loop.
            for (int j=0; j<100000; j++)
            {
                                
                //cout << "New 2-pivot." << endl;
            
                // This will stop us getting stuck in loops.
                // TODO: remove this and fix the underlying problem..
                if (this->dt->is_infinite(c))
                {                
                    cout << "Walked outside of convex hull" << endl;
                    return -1;
                }    
                                
                // index of the facet we just walked through.
                i = c->index(prev);                            

                // Choose a 2-pivot. It will be the 
                // point opposite the face we walked through
                // to get here.
                Vertex_handle pivot_2 = c->vertex(i);


                // index of the face opposite the 3-pivot.
                int pivot_index = c->index(pivot_3);
                

                                            
                // As usual, we choose a random direction
                // and skip two tests in that direction.                
                bool going_cw = random.get_bool();                

                // The first direction is reversed when we walk into a pivot
                // for the first time.
                
                prev = c;
                //cout << "Skipping first test..." << endl;
                c    = c->neighbor( this->get_neighbor(!going_cw, pivot_index, i) );
                this->addToWalk(c);                
                i           = c->index(prev);
                pivot_index = c->index(pivot_3);                                        

        
                prev = c;
                //cout << "Skipping second test..." << endl;
                c    = c->neighbor( this->get_neighbor(going_cw, pivot_index, i) );                        
                this->addToWalk(c);                
                i           = c->index(prev);                        
                pivot_index = c->index(pivot_3);

/*                        
                prev = c;
                cout << "Skipping third test..." << endl;
                c    = c->neighbor( this->get_neighbor(going_cw, pivot_index, i) );                        
                addToWalk(c);                
                i           = c->index(prev);                        
                pivot_index = c->index(pivot_3);
                
                */
                //cout << "Testing direction." << endl;
                // Now decide which direction we are going to 
                // continue in, we may have to change the direction.            
                if ( this->same_side_as_cell(p0, c, this->get_neighbor(going_cw, pivot_index, i)) )
                {
                    //cout << "Switching direction." << endl;
                    // 'pretend' the previous face was coming from the other direction.                    
                    prev = c->neighbor( this->get_neighbor(going_cw, pivot_index,i));                    
                    going_cw = !going_cw;
                } else {
                    //cout << "Direction seems ok." << endl;
                    prev = c;
                    c    = c->neighbor(this->get_neighbor(going_cw, pivot_index, i));    
                    this->addToWalk(c);                                        
                }

                // Continue in this direction
                for (int j=0; j<1000; j++)
                {
                    i           = c->index(prev);                        
                    pivot_index = c->index(pivot_3);
                
                    if (this->dt->is_infinite(c))
                    {
                        cout << "Walked outside of convex hull whilst going around Pivot" << endl;
                        return -1;
                    }
                            
                    if (! this->same_side_as_cell(p0, c, this->get_neighbor(going_cw, pivot_index, i)) )
                    {
                        //cout << "Walking around 2-pivot" << endl;
                        // Continue in this direction
                        prev = c;
                        c    = c->neighbor(this->get_neighbor(going_cw, pivot_index, i));
                        this->addToWalk(c);                        
                    } else {
                        
                        //cout << "Found 2-sink" << endl;
                        
                        // We found the sink, do one more test
                        // to see if it is the end point of this walk.                        
                        if ( this->same_side_as_cell(p0, c, this->get_neighbor(!going_cw, pivot_index, i)) )
                        {
                            //cout << "Found 3-sink." << endl;
                            // We found the sink of this 2D walk.
                            // Now, either the point is in this cell, or we need to move
                            // to the next 3-pivot:
                            
                            if (this->same_side_as_cell(p0, c, pivot_index))
                            {
                                // We are done //
                                if ( this->dt->locate(p0) == c )
                                {
                                    //cout << "Found query point" << endl;
                                    return 1;
                                } else{
                                    cerr << "TEST FAILED" << endl;
                                    return -1;
                                }                                                                                    
                            }
                            
                            //cout << "End of this 3-pivot" << endl;
                            
                            // Walk to the adjacent cell which does not share the same 3-pivot.
                            prev = c;
                            c    = c->neighbor(pivot_index);
                            this->addToWalk(c);                            
                            done = true;
                            break;
                        }
                        
                        // To get to this point means that we have just finished the 2-pivot, 
                        // and need to walk through to the next 2-pivot.
                        prev = c;                        
                        c    = c->neighbor(c->index(pivot_2));
                        this->addToWalk(c);                        
                        break;                                                      
                    }                                                                                            
                }    
            if (done)
                break;                    
            }            
        }        
        return -1;
    }
};

/*****************************************************************************/
// We pass a cell a facet (defined by a cell and the face index)
// and return true if the point is on the same side as the cell
// relative to the hyperplane deined by the facet.

template <typename T>
bool Walk_3d<T>::same_side_as_cell(Point p, Cell_handle cell, int index)
{
    // We uses cases for efficiency
    switch(index)
    {
        case 0:
        return orientation( cell->vertex(1)->point(),
                            cell->vertex(3)->point(),
                            cell->vertex(2)->point(),        
                            p                            ) == CGAL::POSITIVE;

        case 1:
        return orientation( cell->vertex(0)->point(),
                            cell->vertex(2)->point(),
                            cell->vertex(3)->point(),        
                            p                            ) == CGAL::POSITIVE;
        case 2:
        return orientation( cell->vertex(0)->point(),
                            cell->vertex(3)->point(),
                            cell->vertex(1)->point(),        
                            p                            ) == CGAL::POSITIVE;
        case 3:
        return orientation( cell->vertex(0)->point(),
                            cell->vertex(1)->point(),
                            cell->vertex(2)->point(),        
                            p                            ) == CGAL::POSITIVE;
        
        default:
            std::cerr << "Invalid face index requested." << endl;
            
                exit(1);
    }
}

/*****************************************************************************/
// Given the index of the face opposite the current 3-pivot,
// and the index of the face we just walked through, give the index of
// the ccw face.


template <typename T>
int  Walk_3d<T>::get_neighbor(bool cw, int pivot_index, int last_face)
{
    bool error = false;
    
    if (cw)    
    {            
        // There is probably a much nicer way than this, but 
        // just use conditionals for now.
        if (pivot_index == 3)
        {
            switch(last_face)
            {
                case 0: return 2;
                case 1: return 0;
                case 2: return 1;
                default: error=true;                
            }
        }
            else if (pivot_index == 2)
        {
            switch(last_face)
            {
                case 0: return 1;
                case 1: return 3;
                case 3: return 0;
                default: error=true;                
            }
        }
            else if (pivot_index == 1)
        {
            switch(last_face)
            {
                case 0: return 3;
                case 2: return 0;
                case 3: return 2;
                default: error=true;                
            }
        }
            else if (pivot_index == 0)
        {
            switch(last_face)
            {
                case 1: return 2;
                case 2: return 3;
                case 3: return 1;
                default: error=true;                
            }
        }
        
    } else { // DIRECTION == CW // 
        
        if (pivot_index == 3)
        {
            switch(last_face)
            {
                case 0: return 1;
                case 1: return 2;
                case 2: return 0;
                default: error=true;                
            }
        }
            else if (pivot_index == 2)
        {
            switch(last_face)
            {
                case 0: return 3;
                case 1: return 0;
                case 3: return 1;
                default: error=true;                
            }
        }
            else if (pivot_index == 1)
        {
            switch(last_face)
            {
                case 0: return 2;
                case 2: return 3;
                case 3: return 0;
                default: error=true;                
            }
        }
            else if (pivot_index == 0)
        {
            switch(last_face)
            {
                case 1: return 3;
                case 2: return 1;
                case 3: return 2;
                default: error=true;                
            }
        }
    }

    if (error)
    {
        cerr << "Invalid face request. last face:" << last_face << " pivot: " << pivot_index << endl;
        exit(1);
    }    
    return -1;

    
}

/*****************************************************************************/

template <typename T>
class VisibilityWalk_3d : public Walk_3d<T>
{
    
    typedef typename T::Cell                            Cell;
    typedef typename T::Cell_handle                     Cell_handle;
    typedef typename T::Geom_traits                     Gt;        
    typedef typename T::Point                           Point;                
    
    
public:
    VisibilityWalk_3d(T* dt) : Walk_3d<T>(dt) {}
    
    int do_walk(Point p, Cell_handle f=Cell_handle())
    {

        this->clear();
        // The user did not provide a face handle. So just use the infinite 
        // face.
        // TODO: Fix this, seems to cause loops.
        if (f==Cell_handle())
            f=this->dt->infinite_cell();

        // This is where we store the current face.
        Cell_handle    c = f;  
        Cell_handle prev = c;  

        // Create a random number generator.
        CGAL::Random random(time(NULL));


        // ** FIND FIRST FACE ** //
        for (int x=0; x<4; x++)
        {
            // If this is not the facet we have just walked through
            if (! this->same_side_as_cell(p, c, x) )
            {
                c    = c->neighbor( x );        
                this->addToWalk(c);                
                break;
            }
            // TODO: What if the first cell contains the point?
        }
        // ** END OF FIND FIRST FACE ** //


        // Loop until we find our destination point.
        // Temporarily bound number of iterations.
        while(1) //for (int i=0; i<2000; i++)
        { 
            // This will stop us getting stuck in loops.
            if (this->dt->is_infinite(c))
            {                
                cout << "Walked outside of convex hull" << endl;
                return -1;
            }    

            int i = c->index(prev);

            bool done = true;
            for (int x=0; x<4; x++)
            {
                // If this is not the facet we have just walked through
                if (x!=i)
                {
                    if (! this->same_side_as_cell(p, c, x) )
                    {
                        prev = c;
                        c    = c->neighbor( x );
                        this->addToWalk(c);                        
                        done = false;
                        break;
                    }
                }
            }

            // If we reach this point, we know that we have found
            // the final point.
            if (done)
            {
                if ( this->dt->locate(p) == c )
                {
                    //cout << "Visibility walk found query point." << endl;
                    return 1;
                } else{
                    cerr << "TEST FAILED" << endl;
                    return -1;;
                    //exit(1);
                }
                break;
            }                                    
        }            
    }
};

/*****************************************************************************/

template <typename T>
Walk_3d<T>::Walk_3d(T* dt)
{
    this->dt = dt;
    o_count  = 0;
}

/*****************************************************************************/

template <typename T>
int Walk_3d<T>::getNumOrientationsPerformed()
{
    return o_count;
}

/*****************************************************************************/

template <typename T>
int Walk_3d<T>::getNumTrianglesVisited()
{
    return faces.size();
}

/*****************************************************************************/

template <typename T>  
void Walk_3d<T>::addToWalk(Cell_handle c)
{
    faces.push_back(c);
}

/*****************************************************************************/

template <typename T>
vector<typename T::Cell_handle> * Walk_3d<T>::getFaces()
{
    return &faces;
}

/*****************************************************************************/

// Overload the orientation predicate to enable statistics gathering. 
template <typename T>  
CGAL::Orientation Walk_3d<T>::orientation(Point p, Point q, Point r, Point s)
{
    o_count++;    
    return CGAL::orientation(p,q,r,s);    
}

/*****************************************************************************/
template <typename T>
void Walk_3d<T>::clear()
{
    faces.clear();
    o_count =0;
}
    
/*****************************************************************************/

#endif

