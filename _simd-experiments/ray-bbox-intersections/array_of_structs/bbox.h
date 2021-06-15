#ifndef RAY_BBOX_INTERSECTIONS_BBOX_H
#define RAY_BBOX_INTERSECTIONS_BBOX_H

#include "vector3.h"

#include <cassert>

template <typename T>
class BBox {
protected:

    std::array<Vector3<T>, 2> b;

public:

    BBox(const Vector3<T> &min, const Vector3<T> &max) : b{min, max} {
        assert(min.x() <= max.x());
        assert(min.y() <= max.y());
        assert(min.z() <= max.z());
    }

    const decltype(b) &bounds() const { return b; }

};

#endif //RAY_BBOX_INTERSECTIONS_BBOX_H
