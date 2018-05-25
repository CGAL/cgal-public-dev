#include "types.h"

// a sphere with unit radius
FT function_sphere(const Point& query)
{ 
	Point origin(CGAL::ORIGIN);
	return CGAL::squared_distance(query, origin) - 1.0; 
}