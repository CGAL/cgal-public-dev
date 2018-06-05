#ifndef _RANDOM_
#define _RANDOM_ 1
#define PI 3.14159265
#include <cmath>
inline
double random_double(const double min, const double max)
{
	double range = max - min;
	return min + (double(rand()) / double(RAND_MAX)) * range;
}


template <class Vector>
Vector random_unit_vec(const double scale = 1.0)
{
	double dx = random_double(-scale, scale);
	double dy = random_double(-scale, scale);
	double dz = random_double(-scale, scale);
	Vector vec(dx, dy, dz);
	return vec / std::sqrt(vec * vec);
}


template <typename Point>
Point random_point_on_sphere(double radius){
	double theta = random_double(0, 180);
	double phi = random_double(0, 360);
	double x = radius * std::sin(theta * PI / 180.0) * std::cos(phi * PI / 180.0);
	double y = radius * std::sin(theta * PI / 180.0) * std::sin(phi * PI / 180.0);
	double z = radius * std::cos(theta * PI / 180.0);
	return Point(x, y, z);
}

#endif
