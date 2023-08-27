#ifndef SKIP_CONDITIONS_H
#define SKIP_CONDITIONS_H



namespace CGAL {

enum SkipResult
{
    DO_NOT_SKIP = 0,
    SKIP = 1
};

namespace Surface_mesh_approximate_shortest_path_3 {

// skip conditions for dual queue system
/*
class Skip_condition // I'd like to derive all the SkipConditions from one base class, which is the one right here
{
public:
    Skip_condition() {};

    bool operator() ()
    {
        bool to_be_skipped = false;

        return to_be_skipped;
    }
};
*/

class Never_skip_condition
{
public:
    Never_skip_condition() {};

    SkipResult operator() ()
    {
        return CGAL::DO_NOT_SKIP;
    }
};

}

}

#endif // SKIP_CONDITIONS_H
