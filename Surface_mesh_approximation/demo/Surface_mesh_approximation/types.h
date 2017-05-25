#ifndef PCA_DEMO_TYPES_H
#define PCA_DEMO_TYPES_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <boost/iterator/transform_iterator.hpp> 

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Segment_3 Segment;
typedef Kernel::Triangle_3 Triangle;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Point_iterator Point_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Edge_iterator Edge_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

#endif // PCA_DEMO_TYPES_H
