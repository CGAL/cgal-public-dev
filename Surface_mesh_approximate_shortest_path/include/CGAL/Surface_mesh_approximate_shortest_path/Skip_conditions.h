#ifndef SKIP_CONDITIONS_H
#define SKIP_CONDITIONS_H

namespace CGAL {

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

    bool operator() ()
    {
        return false;
    }
};

}

}

#endif // SKIP_CONDITIONS_H
