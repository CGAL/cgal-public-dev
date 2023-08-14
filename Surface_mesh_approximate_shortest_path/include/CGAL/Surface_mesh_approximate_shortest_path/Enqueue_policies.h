#ifndef ENQUEUE_POLICIES_H
#define ENQUEUE_POLICIES_H

namespace CGAL {

namespace Surface_mesh_approximate_shortest_path_3 {

// Enqueue strategies for dual queue system
/*
template<class Custom_enqueue_strategy>
class Enqueue_strategy
{
public:
    typedef typename Custom_enqueue_strategy::Input_struct Input_struct;

public:
    Enqueue_strategy() {};

    bool operator() (Input_struct input_struct)
    {
        return Custom_enqueue_strategy::operator() (input_struct);
    }
};
*/

class Always_enqueue_in_A
{
public:
    Always_enqueue_in_A() {};

    bool operator() ()
    {
        return true; // we always enqueue in A
    }
};

class Always_enqueue_in_B
{
public:
    Always_enqueue_in_B() {};

    bool operator() ()
    {
        return false; // we always enqueue in B
    }
};

class Static_speed_limiter// : public Enqueue_strategy
{
public:
    Static_speed_limiter() {};

    bool operator() (double geodesic_dist, double geodesic_radius)
    {
        return (geodesic_dist < geodesic_radius);
    }
};

}

}

#endif // ENQUEUE_POLICIES_H
