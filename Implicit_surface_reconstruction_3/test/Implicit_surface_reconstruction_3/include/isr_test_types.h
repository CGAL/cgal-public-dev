#ifndef ISR_TEST_TYPES_H
#define ISR_TEST_TYPES_H
 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef std::pair<Point, Vector> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Triangle_3 Triangle;
typedef std::list<Point_with_normal> PwnList;
 
#endif // ISR_TEST_TYPES_H