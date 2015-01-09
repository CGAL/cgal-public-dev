#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2                                     Point;
typedef CGAL::Delaunay_triangulation_2<Kernel>              Delaunay;
typedef Delaunay::Vertex_handle                             Vertex_handle;
typedef Delaunay::Face_handle                               Face_handle;
typedef Delaunay::Edge                                      Edge;