#include <CGAL/Simple_cartesian.h>

#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Circulator/Circulator_concepts.h>
#include <CGAL/Surface_mesh.h>

#include <boost/foreach.hpp>
#include <boost/concept/assert.hpp>

#include <iostream>
#include <iterator>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>                  PolygonMesh;
typedef boost::graph_traits<PolygonMesh>                     GraphTraits;

typedef GraphTraits::vertex_descriptor                       vertex_descriptor;
typedef GraphTraits::halfedge_descriptor                     halfedge_descriptor;
typedef GraphTraits::edge_descriptor                         edge_descriptor;
typedef GraphTraits::face_descriptor                         face_descriptor;
typedef GraphTraits::in_edge_iterator                        in_edge_iterator;
typedef GraphTraits::out_edge_iterator                       out_edge_iterator;

typedef CGAL::Vertex_around_target_circulator<PolygonMesh>   vertex_around_target_circulator;
typedef CGAL::Halfedge_around_source_circulator<PolygonMesh> halfedge_around_source_circulator;
typedef CGAL::Halfedge_around_target_circulator<PolygonMesh> halfedge_around_target_circulator;
typedef CGAL::Halfedge_around_face_circulator<PolygonMesh>   halfedge_around_face_circulator;
typedef CGAL::Face_around_target_circulator<PolygonMesh>     face_around_target_circulator;
typedef CGAL::Face_around_face_circulator<PolygonMesh, CGAL::VERTEX_INCIDENT>       face_around_face_circulator;

typedef CGAL::Vertex_around_target_iterator<PolygonMesh>     vertex_around_target_iterator;
typedef CGAL::Halfedge_around_source_iterator<PolygonMesh>   halfedge_around_source_iterator;
typedef CGAL::Halfedge_around_target_iterator<PolygonMesh>   halfedge_around_target_iterator;
typedef CGAL::Halfedge_around_face_iterator<PolygonMesh>     halfedge_around_face_iterator;
typedef CGAL::Face_around_face_iterator<PolygonMesh>         face_around_face_iterator;
typedef CGAL::Face_around_face_iterator<PolygonMesh, CGAL::VERTEX_INCIDENT>
                                                             vertex_incident_face_around_face_iterator;

int main(int, char**)
int main(int argc, char** argv)
{
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<halfedge_around_face_circulator>)) CGAL_UNUSED;
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<halfedge_around_target_circulator>)) CGAL_UNUSED;
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<vertex_around_target_circulator>)) CGAL_UNUSED;
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<face_around_target_circulator>)) CGAL_UNUSED;
  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<halfedge_around_source_circulator>)) CGAL_UNUSED;

  BOOST_CONCEPT_ASSERT((CGAL::Concepts::BidirectionalCirculator<halfedge_around_source_circulator>)) CGAL_UNUSED;

   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<face_around_face_iterator>)) CGAL_UNUSED;
   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<halfedge_around_face_iterator>)) CGAL_UNUSED;
   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<halfedge_around_target_iterator>)) CGAL_UNUSED;
   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<vertex_around_target_iterator>)) CGAL_UNUSED;

   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<in_edge_iterator>)) CGAL_UNUSED;
   BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<out_edge_iterator>)) CGAL_UNUSED;

  std::ifstream in("data/generic.off");
  PolygonMesh P;
  in >> P;

  halfedge_descriptor hd = *halfedges(P).first;

  // Circulators
  {
    std::cout << "vertices around target circulator -- "
              << "target: " << get(CGAL::vertex_point, P, target(hd,P)) << std::endl;
    vertex_around_target_circulator vatc(hd,P), done(vatc);

    do {
      std::cout << get(CGAL::vertex_point, P, *vatc) << std::endl;
      ++vatc;
    } while(vatc != done);
  }

  {
    std::cout << "halfedges around source circulator -- "
              << "source: " << get(CGAL::vertex_point, P, source(hd,P)) << std::endl;
    halfedge_around_source_circulator hasc(hd,P), done(hasc);

    vertex_descriptor vd = source(hd,P);
    do {
      halfedge_descriptor hd2 = *hasc;
      assert(source(hd2,P) == vd);
      std::cout << get(CGAL::vertex_point, P, target(*hasc,P)) << std::endl;
      ++hasc;
    } while(hasc != done);
  }

  {
    std::cout << "halfedges around target circulator -- "
              << "target: " << get(CGAL::vertex_point, P, target(hd,P)) << std::endl;
    halfedge_around_target_circulator hatc(hd,P), done(hatc);
    vertex_descriptor vd = target(hd,P);

    do {
      halfedge_descriptor hd2 = *hatc;
      assert(target(hd2,P) == vd);
      std::cout << get(CGAL::vertex_point, P, source(*hatc,P)) << std::endl;
      ++hatc;
    } while(hatc != done);
  }

  {
    std::cout << "halfedges around face circulator --" << std::endl;
    std::cout << "face: " << face(hd,P) << std::endl;
    halfedge_around_face_circulator hafc(hd,P), done(hafc);

    do {
      std::cout << get(CGAL::vertex_point, P, target(*hafc,P)) << std::endl;
      assert(face(*hafc,P) == face(hd,P));
      ++hafc;
    } while(hafc != done);
  }

  {
    std::cout << "face around target circulator -- "
              << "target: " << get(CGAL::vertex_point, P, target(hd,P)) << std::endl;
    face_around_target_circulator fatc(hd,P), done(fatc);
    do {
      face_descriptor fd = *fatc;
      if(fd == boost::graph_traits<PolygonMesh>::null_face())
        std::cout << "null face" << std::endl;
      else
        std::cout << "face: " << fd << std::endl;
      ++fatc;
    } while(fatc != done);
  }

  {
    std::cout << "face around face circulator -- "
              << "face: " << face(hd,P) << std::endl;;
    face_around_face_circulator fafc(hd,P), done(fafc);
    do {
      face_descriptor fd = *fafc;
      if(fd == boost::graph_traits<PolygonMesh>::null_face())
        std::cout << "null face" << std::endl;
      else
        std::cout << "face: " << fd << std::endl;
      ++fafc;
    } while(fafc != done);
  }

  // Iterators
  {
    std::cout << "vertices around target iterator -- "
              << "target: " << get(CGAL::vertex_point, P, target(hd,P)) << std::endl;
    vertex_around_target_iterator vit, end;
    boost::tie(vit, end) = vertices_around_target(hd,P);
    assert(std::distance(vit, end) == 3);
    while(vit != end) {
      vertex_descriptor vd = *vit;
      std::cout << get(CGAL::vertex_point, P, vd) << std::endl;
      ++vit;
    }
  }

  {
    std::cout << "halfedges around source iterator -- "
              << "source: " << get(CGAL::vertex_point, P, source(hd,P)) << std::endl;
    halfedge_around_source_iterator hasit, end;
    vertex_descriptor vd = source(hd,P);
    boost::tie(hasit, end) = halfedges_around_source(hd,P);
    assert(std::distance(hasit, end) == 5);
    while(hasit != end) {
      halfedge_descriptor hd2 = *hasit;
      assert(source(hd2,P) == vd);
      std::cout << get(CGAL::vertex_point, P, target(hd2,P)) << std::endl;
      ++hasit;
    }
  }

  {
    std::cout << "halfedges around target iterator -- "
              << "target: " << get(CGAL::vertex_point, P, target(hd,P)) << std::endl;
    halfedge_around_target_iterator hatit, end;
    vertex_descriptor vd = target(hd,P);
    boost::tie(hatit, end) = halfedges_around_target(hd,P);
    assert(std::distance(hatit, end) == 3);
    while(hatit != end) {
      halfedge_descriptor hd2 = *hatit;
      assert(target(hd2,P) == vd);
      std::cout << get(CGAL::vertex_point, P, source(hd2,P)) << std::endl;
      ++hatit;
    }
  }

  {
    std::cout << "halfedges around face iterator -- "
              << "face: " << face(hd,P) << std::endl;
    halfedge_around_face_iterator hafit, end;
    boost::tie(hafit,end) = halfedges_around_face(hd,P);
    assert(std::distance(hafit, end) == 8);
    while(hafit != end) {
      halfedge_descriptor hd2 = *hafit;
      assert(face(hd2,P) == face(hd,P));
      std::cout << get(CGAL::vertex_point, P, target(hd2,P)) << std::endl;
      ++hafit;
    }
  }

  {
    std::cout << "in_edge iterator -- "
              << "target: " << get(CGAL::vertex_point, P, target(hd,P)) << std::endl;
    in_edge_iterator ohi, end;
    boost::tie(ohi, end) = in_edges(target(hd,P),P);
    assert(std::distance(ohi, end) == 3);
    for(; ohi != end; ++ohi){
      edge_descriptor ed = *ohi;
      halfedge_descriptor hd2 = halfedge(ed,P);
      std::cout << get(CGAL::vertex_point, P, target(hd2,P)) << std::endl;
    }
  }

  {
    std::cout << "out_edge iterator -- "
              << "target: " << get(CGAL::vertex_point, P, target(hd,P)) << std::endl;
    out_edge_iterator ohi, end;
    boost::tie(ohi, end) = out_edges(target(hd,P),P);
    assert(std::distance(ohi, end) == 3);
    for(; ohi != end; ++ohi){
      edge_descriptor ed = *ohi;
      halfedge_descriptor hd2 = halfedge(ed,P);
      std::cout << get(CGAL::vertex_point, P, source(hd2,P)) << std::endl;
    }
  }

  {
    std::cout << "out_edges: " << std::endl;
    BOOST_FOREACH(edge_descriptor ed, out_edges(target(hd,P),P)){
      halfedge_descriptor hd2 = halfedge(ed,P);
      std::cout << get(CGAL::vertex_point, P, target(hd2,P)) << std::endl;
    }
  }

  {
    std::cout << "faces around face iterator -- "
              << "face: " << face(hd,P) << std::endl;
    face_around_face_iterator fafit, end;
    boost::tie(fafit, end) = faces_around_face(hd,P);
    assert(std::distance(fafit, end) == 8);
    while(fafit != end) {
      face_descriptor fd = *fafit;
      if(fd == boost::graph_traits<PolygonMesh>::null_face())
        std::cout << "null face" << std::endl;
      else
        std::cout << "face: " << fd << std::endl;
      ++fafit;
    }
  }

  {
    std::cout << "all faces around face " << face(hd,P) << std::endl;
    vertex_incident_face_around_face_iterator afafit, end;
    boost::tie(afafit, end) = vertex_incident_faces_around_face(hd,P);
    while(afafit != end) {
      face_descriptor fd = *afafit;
      if(fd == boost::graph_traits<PolygonMesh>::null_face())
        std::cout << "~~~~~ null face" << std::endl;
      else
        std::cout << "~~~~~ face: " << fd << std::endl;
      ++afafit;
    }
  }

  return 0;
}
