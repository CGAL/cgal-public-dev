#ifndef RAY_BBOX_INTERSECTIONS_BBOX_H
#define RAY_BBOX_INTERSECTIONS_BBOX_H

#include "vector3.h"

#include <cassert>

class BBox {
private:

    std::array<Vector3, 2> b;

public:

    BBox(const Vector3 &min, const Vector3 &max) : b{min, max} {
        assert(min.x() <= max.x());
        assert(min.y() <= max.y());
        assert(min.z() <= max.z());
    }

    const decltype(b) &bounds() const { return b; }

};

#endif //RAY_BBOX_INTERSECTIONS_BBOX_H
