#include <list>
#include <fstream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/single_mold_translational_casting_3.h>
#include <CGAL/Set_movable_separability_3/lp_wrapper.h>
//#include "../../include/CGAL/Set_movable_separability_3/lp_wrapper.h"
#include <CGAL/boost/graph/named_function_params.h>
#include <algorithm>


//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef CGAL::Polyhedron_3<Kernel>                          Polyhedron;

typedef Kernel::Direction_3                                 Direction_3;
typedef Kernel::Point_3                                 Point_3;
typedef Kernel::Vector_3                                Vector_3;

typedef Polyhedron::Halfedge_handle                         Halfedge_handle;
typedef Polyhedron::Facet_handle                            Facet_handle;
typedef Polyhedron::Vertex_handle                           Vertex_handle;

typedef std::list<Direction_3>                              Direction_patch;
//typedef std::pair<boost::graph_traits<Polyhedron>::face_descriptor,
//                  Direction_patch>                          Top_facet;
//typedef std::pair< std::vector<boost::graph_traits<Polyhedron>::face_descriptor>,
//                  Direction_3>                          Top_facet;
typedef std::pair< std::vector<Polyhedron::Facet_const_iterator>,
                  Direction_3>                          Top_facet;
//typedef std::pair<unsigned int,
//                  Direction_3>                          Top_facet;
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

#include "../../include/CGAL/Set_movable_separability_3/Utils.h"
#include "../../include/CGAL/Set_movable_separability_3/coveringset_finder.h"

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "polyhedron.off";
  typename CGAL::Exact_predicates_exact_constructions_kernel::Plane_3 pp;
  typename CGAL::Exact_predicates_exact_constructions_kernel::Line_2 ppp;

  std::ifstream input(filename);
  Polyhedron poly;
  if (!input || !(input >> poly) || poly.empty() ) {
    std::cerr <<filename << " is not a valid off file." << std::endl;
    return -1;
  }


  int outLength=0;
  unsigned int outIndexs[6];

  Direction_3 outDirection; //2 hemisphere that the intersection of their boundary is a point that wasn't covered
  bool outDirectionExists;




  CGAL::cgal_bgl_named_params<bool, CGAL::all_default_t, boost::no_property>
    params;
  std::list<Top_facet> top_facets;

  std::transform( poly.facets_begin(), poly.facets_end(), poly.planes_begin(),
                      Plane_equation());

  SMS::single_mold_translational_casting_3(poly, params,
                                           std::back_inserter(top_facets));
  std::cout<<top_facets.size()<<" valid top facets"<<std::endl;

//  for( std::list<Top_facet>::iterator it=top_facets.begin();it!=top_facets.end();++it)
//    {
//      std::cout<<"facet ? direction "<<it->second<<std::endl;
//    }


  return 0;
}
