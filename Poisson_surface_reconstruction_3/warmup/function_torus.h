#include "types.h"

// from
// http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/geometry/simple.html
// (x^2 + y^2 + z^2 - (R^2 + r^2))^2 = 4 R^2 (r^2 - z^2)

FT function_torus(const Point& query)
{
	Point origin(CGAL::ORIGIN);
	return pow((4 - std::sqrt(query.x() * query.x() + query.y() * query.y())),2) + query.z()*query.z()  - 1.0; 
	}
