#ifndef TYPES_H
#define TYPES_H

// kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/trace.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <deque>
#include <cstdlib>
#include <fstream>
#include <math.h>

#include <boost/foreach.hpp>

#include <CGAL/disable_warnings.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// basic types from kernel
typedef Kernel::FT FT;
typedef CGAL::Bbox_3 Bbox;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Point_3 Point;
typedef Kernel::Sphere_3 Sphere;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Segment_3 Segment;
typedef Kernel::Triangle_3 Triangle;

// 3D mesh generation

// from http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Mesh_3/Chapter_main.html
/*
#include "c3t3.h"
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
*/

typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef std::deque<Point_with_normal> PointList;

// polyhedron
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

// Poisson implicit function
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;


// implicit function

/*typedef FT (Function) (const Point&);

typedef CGAL::Implicit_mesh_domain_3<Function, Kernel> Mesh_domain;

typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
//typedef C3T3<Kernel, Tr> C3t3;

typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
*/

using namespace CGAL::parameters;

#endif // TYPES_H
