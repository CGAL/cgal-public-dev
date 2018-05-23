#include "types.h"

// see, e.g.,
// http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/geometry/simple.html

FT function_ellipsoid(const Point& query)
{
	return query.x() * query.x() + query.y() * query.y()/4.0 + query.z() * query.z()/9.0 - 1.0;
}
