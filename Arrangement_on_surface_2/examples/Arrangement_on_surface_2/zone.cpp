//! \file examples/Arrangement_on_surface_2/zone.cpp
// Computing the zone of a line in an arrangement of five lines.

#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>                   Traits_2;
typedef Traits_2::Point_2                                   Point_2;
typedef Traits_2::Line_2                                    Line_2;
typedef CGAL::Arrangement_2<Traits_2>                       Arrangement_2;

int main()
{
  // Construct the arrangement.
  Arrangement_2 arr;
  insert(arr, Line_2(Point_2(4.2,0), Point_2(5.8,7)));
  insert(arr, Line_2(Point_2(0,0), Point_2(7.985,7)));
  insert(arr, Line_2(Point_2(0,1), Point_2(10,4)));
  insert(arr, Line_2(Point_2(0,5), Point_2(10,2)));
  insert(arr, Line_2(Point_2(5,7), Point_2(10,1)));

  // Compute the zone
  Line_2 zone_line(Point_2(0,4), Point_2(10,3));
  std::list<CGAL::Object> zone_elements;
  CGAL::zone(arr, zone_line, std::back_inserter(zone_elements));

  std::cout << "The zone of the line (" << zone_line << ") cmoprises "
            << zone_elements.size() << " features: " << std::endl;
  std::list<CGAL::Object>::const_iterator it;
  for (it = zone_elements.begin(); it != zone_elements.end(); ++it) {
    CGAL::Object obj = *it;
    typename Arrangement_2::Vertex_handle v;
    typename Arrangement_2::Halfedge_handle  e;
    typename Arrangement_2::Face_handle f;

    if (CGAL::assign(f, obj)) {
      if (f->is_unbounded())
        std::cout << "unbounded face." << std::endl;
      else std::cout << "bounded face." << std::endl;
    }
    else if (CGAL::assign(e, obj))
      std::cout << "edge: " << e->curve() << std::endl;
    else if (CGAL::assign(v, obj)) {
      if (v->is_isolated())
        std::cout << "isolated vertex: " << v->point() << std::endl;
      else std::cout << "vertex: " << v->point() << std::endl;
    }
    else CGAL_error_msg( "Invalid object.");
  }

  return 0;
}
