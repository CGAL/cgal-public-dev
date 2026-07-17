#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/intersections.h>



// Unified global type mappings for the entire pipeline
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Kernel::Point_3 Point3;
typedef Kernel::Triangle_3 Triangle3;