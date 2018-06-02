#ifndef _RANDOM_
#define _RANDOM_ 1

inline
double random_double(const double min, const double max)
{
	double range = max - min;
	return min + (double(rand()) / double(RAND_MAX)) * range;
}


template <class Vector>
Vector random_vec(const double scale)
{
	double dx = random_double(-scale, scale);
	double dy = random_double(-scale, scale);
	double dz = random_double(-scale, scale);
	return Vector(dx, dy, dz);
}


#endif
