// based on Apollonius_graph_traits_2 by Menelaos Karavelas

//    (c) 2007-2009 National and Kapodistrian University of Athens
//    (c) 2009-2011 INRIA Nancy
//    (c) 2011-2012 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED

#ifndef CGAL_DELAUNAY_GRAPH_OF_ELLIPSES_2_H
#define CGAL_DELAUNAY_GRAPH_OF_ELLIPSES_2_H

#include <CGAL/Apollonius_graph_2.h>

#ifndef CGAL_APOLLONIUS_GRAPH_PSEUDO_CIRCLE_DESIGN
#ifndef WITH_PSEUDOCIRCLES_PATCH
#warning This code requires the version of the Apollonius graph package \
from the branch: ^/branches/features/Apollonius_graph_2-general_is_hidden-mkaravel, \
otherwise results with hidden sites will be UNPREDICTABLE!
#endif
#endif

namespace CGAL {

template < class Gt,
	   class Agds = Triangulation_data_structure_2 < 
               Apollonius_graph_vertex_base_2<Gt,true>,
               Triangulation_face_base_2<Gt> >,
	   class LTag = Tag_false>
class Delaunay_graph_of_ellipses_2: public Apollonius_graph_2<Gt,Agds,LTag> 
{
    typedef Apollonius_graph_2<Gt,Agds,LTag> Base;
    
    typedef typename Gt::Ellipse_bisector_2 Ellipse_bisector_2;

//    is this possible?
//    template<class T> class VoronoiEllipsesGraphicsItem;
//    friend class VoronoiEllipsesGraphicsItem<Delaunay_graph_of_ellipses_2<Gt> >;

public:
    // TODO: Maybe create access methods a la Apollonius_graph_2

    typedef Gt Geom_traits;

    typedef typename Gt::Point_2                   Point_2;
    typedef typename Gt::Site_2                    Site_2;

    typedef typename Base::Edge                    Edge;

#ifdef WITH_PSEUDOCIRCLES_PATCH
    typedef typename Base::Vertex_handle           Vertex_handle;
    typedef typename Base::Face_handle             Face_handle;
    typedef typename Base::Vertex                  Vertex;
    typedef typename Base::Face                    Face;

    typedef typename Base::Vertex_circulator       Vertex_circulator;
    typedef typename Base::Edge_circulator         Edge_circulator;
    typedef typename Base::Face_circulator         Face_circulator;

    typedef typename Base::Face_iterator           All_faces_iterator;
    typedef typename Base::Vertex_iterator         All_vertices_iterator;
    typedef typename Base::Edge_iterator           All_edges_iterator;

    typedef typename Base::Face_map                Face_map;
    typedef typename Base::Vertex_map              Vertex_map;

    typedef typename Base::Finite_faces_iterator     Finite_faces_iterator;
    typedef typename Base::Finite_vertices_iterator  Finite_vertices_iterator;
    typedef typename Base::Finite_edges_iterator     Finite_edges_iterator;

    typedef typename Base::Conflict_type            Conflict_type;
    typedef typename Base::List                     List;


    int covered(const Site_2 &p, const Site_2 &q, const Site_2 &r) const
    {
      return this->geom_traits().covered_2_object()(p, q, r);
    }
#endif

    typename Gt::Ellipse_bisector_2 dual(const Edge e) const {
    
      CGAL_triangulation_precondition( !is_infinite(e) );

      if ( this->dimension() == 1 ) {
        Site_2 p = (e.first)->vertex(this->cw(e.second))->site();
        Site_2 q = (e.first)->vertex(this->ccw(e.second))->site();

        return Ellipse_bisector_2(p,q);
      }

      // dimension == 2
      // none of the two adjacent faces is infinite
      if( (!is_infinite(e.first)) &&
          (!is_infinite(e.first->neighbor(e.second))) ) {
        Site_2 p = (e.first)->vertex( this->ccw(e.second) )->site();
        Site_2 q = (e.first)->vertex(  this->cw(e.second) )->site();
        Site_2 r = (e.first)->vertex(     e.second  )->site();
        Site_2 s = this->tds().mirror_vertex(e.first, e.second)->site();
        return Ellipse_bisector_2(p,q,r,s);
      }

      // both of the adjacent faces are infinite
      if ( is_infinite(e.first) &&
           is_infinite(e.first->neighbor(e.second)) )  {
        Site_2 p = (e.first)->vertex(this->cw(e.second))->site();
        Site_2 q = (e.first)->vertex(this->ccw(e.second))->site();
        return Ellipse_bisector_2(p,q);
      }

      // only one of the adjacent faces is infinite
      CGAL_triangulation_assertion( is_infinite( e.first ) ||
				    is_infinite( e.first->neighbor(e.second) )
				    );

      CGAL_triangulation_assertion( !(is_infinite( e.first ) &&
				      is_infinite( e.first->neighbor(e.second) )
				      )
				    );

      CGAL_triangulation_assertion
        (  this->is_infinite( e.first->vertex(e.second) ) ||
           is_infinite( this->tds().mirror_vertex(e.first, e.second) )  );

      Edge ee = e;
      if ( is_infinite( e.first->vertex(e.second) )  ) {
        ee = sym_edge(e);
      }
      Site_2 p = ee.first->vertex( this->ccw(ee.second) )->site();
      Site_2 q = ee.first->vertex(  this->cw(ee.second) )->site();
      Site_2 r = ee.first->vertex(     ee.second  )->site();

      return Ellipse_bisector_2(p,q,r);
    }


    // changes propagated to Apollonius_graph_2 package

#ifdef WITH_PSEUDOCIRCLES_PATCH
    Vertex_handle insert_third(const Site_2& p)
    {
      CGAL_triangulation_precondition( this->number_of_vertices() == 2 );

      Vertex_handle v1(this->finite_vertices_begin());
      Vertex_handle v2(++(this->finite_vertices_begin()));

      if ( this->is_hidden(v1->site(), p) ) {
        v1->add_hidden_site(p);
        return Vertex_handle();
      }
      if ( this->is_hidden(v2->site(), p) ) {
        v2->add_hidden_site(p);
        return Vertex_handle();
      }

      bool t1 = this->is_hidden(p, v1->site());
      bool t2 = this->is_hidden(p, v2->site());

      if ( t1 && !t2 ) {
        v1->add_hidden_site(v1->site());
        v1->set_site(p);
        return v1;
      } else if ( !t1 && t2 ) {
        v2->add_hidden_site(v2->site());
        v2->set_site(p);
        return v2;
      } else if ( t1 && t2 ) {
        v1->add_hidden_site(v1->site());
        v1->add_hidden_site(v2->site());
        v1->set_site(p);
        this->remove_second(v2);
        return v1;
      }

      int c = covered(v1->site(), v2->site(), p);
      if (c == 3) {
          v1->add_hidden_site(p);
          return Vertex_handle();
      } else if (c == 2) {
          v2->add_hidden_site(v2->site());
          v2->set_site(p);
          return v2;
      } else if (c == 1) {
          v1->add_hidden_site(v1->site());
          v1->set_site(p);
          return v1;
      }

      Vertex_handle v = this->data_structure().insert_dim_up(this->infinite_vertex());
      v->set_site(p);

      Face_handle f(this->finite_faces_begin());

      Point_2 p1 = f->vertex(0)->site().point();
      Point_2 p2 = f->vertex(1)->site().point();
      Point_2 p3 = f->vertex(2)->site().point();

      Orientation o =
        this->geom_traits().orientation_2_object()(p1, p2, p3);

      if ( o != LEFT_TURN ) {
        f->reorient();
        for (int i = 0; i < 3; i++) {
          f->neighbor(i)->reorient();
        }
      }

      Conflict_type ct =
        this->finite_edge_conflict_type_degenerated(v1->site(), v2->site(), p);

      if ( ct == Base::NO_CONFLICT ) {
        Oriented_side os =
          side_of_bisector(v1->site(), v2->site(), p.point());

        CGAL_assertion( os != ON_ORIENTED_BOUNDARY );
        Vertex_handle vv = ( os == ON_NEGATIVE_SIDE ) ? v1 : v2;

        Face_circulator fc = this->incident_faces(v);
        while ( true ) {
          Face_handle f(fc);
          int k = f->index(v);
          Vertex_handle vh = f->vertex(this->ccw(k));
          if ( vh == vv ) {
        flip(f, this->cw(k));
        break;
          }
          ++fc;
        }
      } else if ( ct == Base::INTERIOR ) {
        Edge_circulator ec = this->incident_edges(v);

        while ( true ) {
          if ( is_infinite(ec) ) {
        flip(*ec);
        break;
          }
          ec++;
        }
      } else if ( ct == Base::ENTIRE_EDGE ) {
        Face_circulator fc = this->incident_faces(v);

        while ( true ) {
          Face_handle f(fc);
          if ( !is_infinite(f) ) {
        flip(f, f->index(v));
        break;
          }
          ++fc;
        }
      } else if ( ct == Base::BOTH_VERTICES ) {


        Conflict_type ct1 =
          this->finite_edge_conflict_type_degenerated(v1->site(), p, v2->site());

        Edge_circulator ec;
        ec = ( ct1 == Base::INTERIOR ) ? this->incident_edges(v2) : this->incident_edges(v1);
        while ( true ) {
          if ( is_infinite(ec) ) {
        flip(*ec);
        break;
          }
          ec++;
        }
      } else {
        CGAL_assertion( ct == Base::RIGHT_VERTEX || ct == Base::LEFT_VERTEX );
        // do nothing here
      }

      //  CGAL_triangulation_assertion( is_valid() );

      return v;
    }

    Vertex_handle insert(const Site_2& p, Vertex_handle vnear)
    {
      if ( this->number_of_vertices() == 0 ) {
        return this->insert_first(p);
      }
      if ( this->number_of_vertices() == 1 ) {
        return this->insert_second(p);
      }
      if ( this->number_of_vertices() == 2 ) {
        return insert_third(p);
      }

      // first find the nearest neighbor
      Vertex_handle vnearest = nearest_neighbor(p.point(), vnear);

      CGAL_assertion( vnearest != Vertex_handle() );


      // check if it is hidden
      Site_2 wp_nearest = vnearest->site();
      if ( is_hidden(wp_nearest, p) ) {
        vnearest->add_hidden_site(p);
        return Vertex_handle();
      }

      // find the first conflict

      // first look for conflict with vertex
      Face_circulator fc_start = incident_faces(vnearest);
      Face_circulator fc = fc_start;
      Face_handle start_f;
      Sign s;
      do {
        Face_handle f(fc);
        s = incircle(f, p);

        if ( s == NEGATIVE ) {
          start_f = f;
          break;
        }
        ++fc;
      } while ( fc != fc_start );

      // we are not in conflict with an Apollonius vertex, so we have to
      // be in conflict with the interior of an Apollonius edge
      if ( s != NEGATIVE ) {
        Edge_circulator ec_start = incident_edges(vnearest);
        Edge_circulator ec = ec_start;

        bool interior_in_conflict(false);
        Edge e;
        do {
          e = *ec;

          Vertex_handle v1(e.first->vertex(this->ccw(e.second)));
          Vertex_handle v2(e.first->vertex(this->cw(e.second)));
          if (!is_infinite(v1) && !is_infinite(v2)) {
              int c = covered(v1->site(), v2->site(), p);
              if (c == 3) {
                  v1->add_hidden_site(p);
                  return Vertex_handle();
              }
          }

          interior_in_conflict = edge_interior(e, p, false);

          if ( interior_in_conflict ) { break; }
          ++ec;
        } while ( ec != ec_start );

        CGAL_assertion( interior_in_conflict );

        return insert_degree_2(e, p);
      }


      // we are in conflict with an Apollonius vertex; start from that and
      // find the entire conflict region and then repair the diagram
      List l;
      Face_map fm;
      Vertex_map vm;

      // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
      // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
      // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT
      // REGION AND EXPANDS THE CONFLICT REGION.
      initialize_conflict_region(start_f, l);
      expand_conflict_region(start_f, p, l, fm, vm, NULL);

      //  retriangulate_conflict_region(v, l, fm, vm);
      Vertex_handle v = retriangulate_conflict_region(p, l, fm, vm);

      fm.clear();
      vm.clear();


      return v;
    }

    Vertex_handle  insert(const Site_2& p) {
      return insert(p, Vertex_handle());
    }
#endif

};

} //namespace CGAL
#endif // CGAL_DELAUNAY_GRAPH_OF_ELLIPSES_2_H
