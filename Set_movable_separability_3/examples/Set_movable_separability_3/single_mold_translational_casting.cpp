#include <list>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/single_mold_translational_casting_3.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <algorithm>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>                          Polyhedron;
typedef Kernel::Direction_3                                 Direction_3;

typedef Polyhedron::Halfedge_handle                         Halfedge_handle;
typedef Polyhedron::Facet_handle                            Facet_handle;
typedef Polyhedron::Vertex_handle                           Vertex_handle;

typedef std::list<Direction_3>                              Direction_patch;
typedef std::pair<boost::graph_traits<Polyhedron>::face_descriptor,
                  Direction_patch>                          Top_facet;

namespace SMS = CGAL::Set_movable_separability_3;

struct Plane_equation {
    template <class Facet>
    typename Facet::Plane_3 operator()( Facet& f) {
        typename Facet::Halfedge_handle h = f.halfedge();
        typedef typename Facet::Plane_3  Plane;
        return Plane( h->vertex()->point(),
                      h->next()->vertex()->point(),
                      h->next()->next()->vertex()->point());
    }
};

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "polyhedron.off";
  std::ifstream input(filename);
  Polyhedron poly;
  if (!input || !(input >> poly) || poly.empty() ) {
    std::cerr << "Not a valid off file." << std::endl;
    return -1;
  }


  CGAL::cgal_bgl_named_params<bool, CGAL::all_default_t, boost::no_property>
    params;
  std::list<Top_facet> top_facets;

  std::transform( poly.facets_begin(), poly.facets_end(), poly.planes_begin(),
                      Plane_equation());
  SMS::single_mold_translational_casting_3(poly, params,
                                           std::back_inserter(top_facets));

  return 0;
}
