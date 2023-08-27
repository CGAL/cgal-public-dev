#ifndef ENQUEUE_POLICIES_H
#define ENQUEUE_POLICIES_H

namespace CGAL {

enum EnqueueResult
{
    ENQUEUE_IN_A = 1,
    ENQUEUE_IN_B = 0,
    DO_NOT_ENQUEUE = -1,
};

namespace Surface_mesh_approximate_shortest_path_3 {

template<class Kernel>
class Always_enqueue_in_A
{
public:
    typedef typename Kernel::FT         FT;
    typedef typename Kernel::Point_3    Point_3;

public:
    Always_enqueue_in_A() {};

    EnqueueResult operator() (FT geodesic_dist = 0., FT geodesic_radius = 0.,
                             Point_3 prev_face_target = Point_3(),
                             Point_3 curr_face_target = Point_3(),
                             Point_3 overall_target = Point_3())
    {
        return ENQUEUE_IN_A;
    }
};

template<class Kernel>
class Always_enqueue_in_B
{
public:
    typedef typename Kernel::FT         FT;
    typedef typename Kernel::Point_3    Point_3;

public:
    Always_enqueue_in_B() {};

    EnqueueResult operator() (FT geodesic_dist, FT geodesic_radius,
                             Point_3 prev_face_target = Point_3(),
                             Point_3 curr_face_target = Point_3(),
                             Point_3 overall_target = Point_3())
    {
        return ENQUEUE_IN_B; // we always enqueue in B
    }
};

template<class Kernel>
class Static_speed_limiter
{
public:
    typedef typename Kernel::FT         FT;
    typedef typename Kernel::Point_3    Point_3;

public:
    Static_speed_limiter() {};

    EnqueueResult operator() (FT geodesic_dist, FT geodesic_radius,
                              Point_3 prev_face_target = Point_3(),
                              Point_3 curr_face_target = Point_3(),
                              Point_3 overall_target = Point_3())
    {
        if (geodesic_dist < geodesic_radius)
        {
            return ENQUEUE_IN_A;
        }
        else
        {
            return ENQUEUE_IN_B;
        }
    }
};

template<class Kernel>
class Embedding_space_distance_limiter
{
public:
    typedef typename Kernel::FT         FT;
    typedef typename Kernel::Point_3    Point_3;

public:
    Embedding_space_distance_limiter() {};

    CGAL::EnqueueResult operator () (FT geodesic_dist, FT geodesic_radius,
                                   Point_3 prev_face_target = Point_3(),
                                   Point_3 curr_face_target = Point_3(),
                                   Point_3 overall_target = Point_3())
    {
        FT prev_sq_embedding_dist = squared_distance(prev_face_target, overall_target);
        FT sq_embedding_dist = squared_distance(curr_face_target, overall_target);
        FT embedding_dist = sqrt(sq_embedding_dist);

        if (sq_embedding_dist <= prev_sq_embedding_dist + 10.) {
            return CGAL::ENQUEUE_IN_A; // keep enqueuing in A => these halfedges will be prioritized
        }
        else if (geodesic_dist + embedding_dist <= geodesic_radius + 10.) {
            return CGAL::ENQUEUE_IN_B; // enqueue in B => explore later, after we have found a first geodesic distance
        }
        else {
            return CGAL::DO_NOT_ENQUEUE; // this can no longer be a sensible geodesic path anyways
        }
    }
};

}

}

#endif // ENQUEUE_POLICIES_H
