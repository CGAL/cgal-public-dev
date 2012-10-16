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


#ifndef CGAL_DELAUNAY_GRAPH_OF_ELLIPSES_TRAITS_2_H
#define CGAL_DELAUNAY_GRAPH_OF_ELLIPSES_TRAITS_2_H


#include <CGAL/Ellipse_2.h>
#include <CGAL/Side_of_bisector.h>
#include <CGAL/Bitangent.h>
#include <CGAL/Distance_from_bitangent.h>
#include <CGAL/In_circle.h>
#include <CGAL/Ellipse_triplet.h>
#include <CGAL/Ellipse_bisector_2.h>
#include <CGAL/Voronoi_circle.h>

//#include <CGAL/number_type_basic.h>
//#include <CGAL/Apollonius_graph_2/Kernel_wrapper_2.h>

namespace CGAL {

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// the Traits class
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

#ifdef DGET_VERBOSE
#define RETURN2(pname, s1, s2, r) std::cerr << pname << "(s1 = " << s1 << \
    " s2 = " << s2 << ") = " << r << std::endl; \
    return r
#define RETURN3(pname, s1, s2, s3, r) std::cerr << pname << "(s1 = " << s1 << \
    " s2 = " << s2 << " s3 = " << s3 << ") = " << r << std::endl; \
    return r
#define ENTER4(pname, s1, s2, s3, s4) std::cerr << "ENTER " << pname << \
    "(s1 = " << s1 << " s2 = " << s2 << " s3 = " << s3 << " s4 = " << s4 << \
    ")\n"
#define RETURN4(pname, s1, s2, s3, s4, r) std::cerr << pname << \
    "(s1 = " << s1 << " s2 = " << s2 << " s3 = " << s3 << " s4 = " << s4 << \
    ") = " << r << std::endl; \
    return r
#define RETURN5(pname, s1, s2, s3, s4, s5, r) std::cerr << pname << \
    "(s1 = " << s1 << " s2 = " << s2 << " s3 = " << s3 << " s4 = " << s4 << \
    " s5 = " << s5 << ") = " << r << std::endl; \
    return r
#else
#define RETURN2(pname, s1, s2, r) return r
#define RETURN3(pname, s1, s2, s3, r) return r
#define ENTER4(pname, s1, s2, s3, s4)
#define RETURN4(pname, s1, s2, s3, s4, r) return r
#define RETURN5(pname, s1, s2, s3, s4, s5, r) return r
#endif

template < class ET >
class Delaunay_graph_of_ellipses_traits_2
{
public:
  typedef ET                            Ellipse_traits;
  typedef Ellipse_2<ET>                 Site_2;

//  class EmptyObject {
//  };
  
//  typedef EmptyObject Line_2;
//  typedef EmptyObject Ray_2;
//  typedef EmptyObject Segment_2;
//  typedef EmptyObject Object_2;

  typedef CGAL::Ellipse_bisector_2<ET> Ellipse_bisector_2;
  typedef Ellipse_bisector_2 Object_2;
  typedef Ellipse_bisector_2 Line_2;
  typedef Ellipse_bisector_2 Ray_2;
  typedef Ellipse_bisector_2 Segment_2;

  typedef typename ET::QQ               FT;
  typedef typename ET::ZZ               RT;
  typedef typename ET::Kernel::Point_2  Point_2;
  typedef typename ET::Kernel           Geom_traits;

public:
  // PREDICATES
  //-----------
   class Compare_x_2 {
   public:
       typedef Comparison_result result_type;
       typedef Site_2            argument_type;

       inline result_type operator()(const Site_2& s1, const Site_2& s2) const 
       {
           result_type r = CGAL::compare(s1.x_center(), s2.x_center());
           RETURN2("Compare_x", s1.get_id(), s2.get_id(), r);
       }
   };
   
   class Compare_y_2 {
   public:
       typedef Comparison_result result_type;
       typedef Site_2            argument_type;

       inline result_type operator()(const Site_2& s1, const Site_2& s2) const 
       {
           result_type r = CGAL::compare(s1.y_center(), s2.y_center());
           RETURN2("Compare_y", s1.get_id(), s2.get_id(), r);
       }
   };
   
   class Compare_weight_2 {
   public:
       typedef Comparison_result result_type;
       typedef Site_2            argument_type;

       inline result_type operator()(const Site_2& s1, const Site_2& s2) const 
       {
           result_type r = CGAL::compare(s1.minor_axis() * s1.major_axis(), 
                                s2.minor_axis() * s2.major_axis());
           RETURN2("Compare_weight", s1.get_id(), s2.get_id(), r);
       }
   };

  class Orientation_2 {
  public:
      typedef Orientation result_type;
      typedef Point_2      argument_type;
      
      inline Orientation operator()(const Point_2 &p1, const Point_2& p2, 
                             const Point_2& p3) const {
          result_type r = orientation(p1, p2, p3);
          RETURN3("Orientation", p1, p2, p3, r);
      }
  };

  class Is_hidden_2 {
  public:
      typedef bool     result_type;
      typedef Site_2   argument_type;

      enum { Has_three_argument_operator = true };

      inline bool operator()(const Site_2 &p, const Site_2 &q) const {
          VORELL::Bitangent<ET> bt(p, q);
          result_type r = bt.relative_position() == VORELL::HIDDEN;
          // q is hidden inside p
          RETURN2("Is_hidden", p.get_id(), q.get_id(), r);
      }

//      inline int operator()(const Site_2 &p, const Site_2 &q, const Site_2 &r) const {
//          Ellipse_triplet<ET> triplet(p, q, r);
//          result_type res = (triplet.get_cover_index() == 3);
//          RETURN3("Is_hidden3", p.get_id(), q.get_id(), r.get_id(), res);
//      }
  };

#ifdef WITH_PSEUDOCIRCLES_PATCH
  class Covered_2 {
  public:
      typedef int      result_type;
      typedef Site_2   argument_type;

      inline int operator()(const Site_2 &p, const Site_2 &q, const Site_2 &r) const {
          Ellipse_triplet<ET> triplet(p, q, r);
          result_type res = triplet.get_cover_index();
          RETURN3("Covered", p.get_id(), q.get_id(), r.get_id(), res);
      }
  };
#endif

  class Oriented_side_of_bisector_2 {
  public:
      typedef Oriented_side result_type;
      Oriented_side operator()(const Site_2& p1, const Site_2& p2,
			       const Point_2 &p) const {
          // SMALLER ===> ON_POSITIVE SIDE !!!
          result_type r = - Side_of_bisector<ET>()(p1, p2, p.x(), p.y());
          RETURN3("Oriented_side_of_bisector", p1.get_id(), p2.get_id(), p, r);
      }
  };

  class Vertex_conflict_2 {
  public:
      typedef Sign    result_type;
      typedef Site_2  argument_type;
      inline Sign operator()(const Site_2& p1, const Site_2& p2,
		             const Site_2& q) const {
          result_type r;
          // TODO: check this alt. implem?
//          result_type r = CGAL::POSITIVE;
//          Ellipse_triplet<ET> triplet(p1,p2,q);
//          // TODO: check...
//          // speedup with Distance from bitangent for nonintersecting? or already in triplet?
//          if (triplet.get_cover_index() < 3) {
//              int srt = triplet.get_shadow_region_type();
//              if (srt == Ellipse_triplet<ET>::LEFT_VERTEX ||
//                  srt == Ellipse_triplet<ET>::TWO_VERTICES ||
//                  srt == Ellipse_triplet<ET>::ENTIRE_EDGE) {
//                  r = CGAL::NEGATIVE;
//              }
//          }
          result_type s = Distance_from_bitangent<ET>()(p2, p1, q);
          if (s != CGAL::ZERO) r = s;
          else {
              // TODO: works only for non-intersecting ellipses?
              Ellipse_triplet<ET> triplet(p1, p2, q);
              if (triplet.is_shadow_region_connected())
                  r = CGAL::NEGATIVE;
              else
                  r = CGAL::POSITIVE;
              // TODO: implement new way by wrapping tan. of q...
          }
          RETURN3("Vertex_conflict", p1.get_id(), p2.get_id(), q.get_id(), r);
      }
      
      inline Sign operator()(const Site_2& p1, const Site_2& p2,
		             const Site_2& p3, const Site_2& q) const {
          ENTER4("Vertex_conflict", p1.get_id(), p2.get_id(), p3.get_id(), q.get_id());
          result_type r = In_circle<ET>()(p1, p2, p3, q);
          RETURN4("Vertex_conflict", p1.get_id(), p2.get_id(), p3.get_id(), q.get_id(), r);
      }
  };
  
  // adapted from Finite_edge_test_C2 (Apollonius package)
  class Finite_edge_interior_conflict_2 {
  public:
      typedef bool result_type;
      bool operator()(const Site_2& p1, const Site_2& p2, const Site_2& p3,
	              const Site_2& p4, const Site_2& q, bool b) const {

          result_type res;
          
          Ellipse_triplet<ET> t12q(p1, p2, q);
          Ellipse_triplet<ET> t1q2(p1, q, p2);

          int nvc = t12q.get_num_voronoi_circles();

          // both the ccw and cw circles do not exist
          if ( nvc == 0 ) res = b;
          // the ccw circle exists but not the cw
          else if ( nvc == 1 ) res = b;
          // the cw circle exists but not the ccw
          else if ( nvc == -1 ) res = b;
          else {
              // both circles exist

              // check whether the shadow region is connected, i.e., 
              // whether it is of the form (a, b) or (-oo, a) U (b, +oo)

              if ( t12q.is_shadow_region_connected() ) {
                if ( b ) res = true; // S10
                else {
                    Voronoi_circle<ET> vc142(p1, p4, p2); // a
                    Voronoi_circle<ET> vc1q2(p1, q, p2); // b'

                    Comparison_result r = vc142.order_on_boundary(vc1q2);

                    if ( r != SMALLER ) res = false;  // S14
                    else {
                        Voronoi_circle<ET> vc123(p1, p2, p3); // b
                        Voronoi_circle<ET> vc12q(t12q); // a'

                        r = vc123.order_on_boundary(vc12q);

                        res = ( r == LARGER );
                    }
                }
              } else {
                  // the shadow region is of the form (-oo, a) U (b, +oo)
                  if ( !b ) return res = false; // S15
                  else {
                      Voronoi_circle<ET> vc142(p1, p4, p2); // a
                      Voronoi_circle<ET> vc1q2(t1q2); // b'

                      Comparison_result r =vc142.order_on_boundary(vc1q2);

                      if ( r == LARGER ) res = true; // S20
                      else {
                          Voronoi_circle<ET> vc123(p1, p2, p3); // b
                          Voronoi_circle<ET> vc12q(t12q); // a'

                          r = vc123.order_on_boundary(vc12q);

                          res = ( r == SMALLER );
                      }
                  }
              }
          }
          RETURN5("Finite_edge_interior_conflict", p1.get_id(), p2.get_id(), p3.get_id(), p4.get_id(), q.get_id(), res);
      }
      
      // degenerate
      bool operator()(const Site_2& p1, const Site_2& p2, const Site_2& p3,
	              const Site_2& q, bool b) const {
          
          result_type res;
          
          Ellipse_triplet<ET> t12q(p1, p2, q);

          int nvc = t12q.get_num_voronoi_circles();
          
          // both the ccw and cw circles do not exist
          if ( nvc == 0 ) res = b;
          // the ccw circle exists but not the cw
          else if ( nvc == 1 ) res = b;
          // the cw circle exists but not the ccw
          else if ( nvc == -1 ) res = b;
          else {
              // both circles exist

              // check whether the shadow region is connected, i.e., 
              // whether it is of the form (a, b) or (-oo, a) U (b, +oo)

              if ( t12q.is_shadow_region_connected() ) {
                // the shadow region is of the form (a, b)
                if ( b ) res = true; // S10
                else {
                    Voronoi_circle<ET> vc123(p1, p2, p3); // b
                    Voronoi_circle<ET> vc12q(t12q); // a'

                    Comparison_result r = vc123.order_on_boundary(vc12q);

                    res = ( r == LARGER );
                }
              } else {
                  // the shadow region is of the form (-oo, a) U (b, +oo)
                  if ( !b ) res = false; // S15
                  else {
                      Voronoi_circle<ET> vc123(p1, p2, p3); // b
                      Voronoi_circle<ET> vc12q(t12q); // a'

                      Comparison_result r =vc123.order_on_boundary(vc12q);
                      res = ( r == SMALLER );
                  }
              }
          }
          RETURN5("Finite_edge_interior_conflict+bool", p1.get_id(), p2.get_id(), p3.get_id(), q.get_id(), b, res);
      }

      bool operator()(const Site_2& p1, const Site_2& p2,
           	      const Site_2& q, bool b) const {
          
          result_type res;
          Ellipse_triplet<ET> t12q(p1, p2, q);

          int nvc = t12q.get_num_voronoi_circles();
          
          // both the ccw and cw circles do not exist
          if ( nvc == 0 ) {
              res = b;
          }
          // only the ccw circle exists
          else if ( nvc == 1 ) res = false;
          // only the cw circle exists
          else if ( nvc == -1 ) res = false;
          else {
              // both circles exist

              // check whether the shadow region is connected, i.e., wether it is
              // of the form (a, b) or (-oo, a) U (b, +oo)

              res = !b;
          }
          RETURN3("Finite_edge_interior_conflict", p1.get_id(), p2.get_id(), q.get_id(), res);
      }
  };

  // adapted from Infinite_edge_test_C2 (Apollonius package)
  class Infinite_edge_interior_conflict_2 {
      typedef typename ET::AK::Algebraic_real_1 Root;
    
      Bounded_side bounded_side_of_arc(const Site_2& e1, const Site_2& e2, 
                         const Site_2& e3, const Site_2& q, 
                         const bool first) const {
          
          VORELL::Bitangent<ET> bt21(e2, e1);
          VORELL::Bitangent<ET> bt23(e2, e3);
          VORELL::Bitangent<ET> bt2q(e2, q);
          Root p;

          VORELL::Range<Root> arc(bt21.other_external(), bt23.external());
          if (first) p = bt2q.external(); else p = bt2q.other_external();
          
          if (arc.touches(p)) return ON_BOUNDARY;
          if (arc.strictly_contains(p)) return ON_BOUNDED_SIDE;
          return ON_UNBOUNDED_SIDE;
      }
      
  public:
      typedef bool result_type;
  
      bool operator()(const Site_2& p2, const Site_2& p3,
	              const Site_2& p4, const Site_2& q, bool b) const {
      
          result_type res;
          Bounded_side bs1 = bounded_side_of_arc(p3, p2, p4, q, true);
          if ( b ) {
            if ( bs1 == ON_BOUNDARY ) {
	      Bounded_side bs2 = bounded_side_of_arc(p3, p2, p4, q, false);

	      if ( bs2 != ON_BOUNDARY ) {
	        res = ( bs2 != ON_BOUNDED_SIDE );
	      } else res = !b;
            } else res = ( bs1 != ON_BOUNDED_SIDE );
          } else if ( bs1 == ON_BOUNDARY ) {
            Bounded_side bs2 = bounded_side_of_arc(p3, p2, p4, q, false);
            if ( bs2 != ON_BOUNDARY ) {
	      res = ( bs2 == ON_BOUNDED_SIDE );
            } else res = !b;
          } else res = ( bs1 == ON_BOUNDED_SIDE );
          RETURN4("Infinite_edge_interior_conflict", p2.get_id(), p3.get_id(), p4.get_id(), q.get_id(), res);
      }
      
  };
  
  class Is_degenerate_edge_2 {
  public:
      typedef bool result_type;
      typedef Site_2  argument_type;

      bool operator()(const Site_2& p1, const Site_2& p2,
	    	       const Site_2& p3, const Site_2& p4) const {
      
          result_type res;
          Ellipse_triplet<ET> triplet(p1,p2,p3);
          if (triplet.get_num_voronoi_circles() < 1) res = false;
          else {
              Ellipse_triplet<ET> triplet2(p1, p4, p2);
              if (triplet2.get_num_voronoi_circles() < 1) res = false;
              else {
                  Voronoi_circle<ET> vc123(triplet);
                  Voronoi_circle<ET> vc142(triplet2);
                  res = vc142.order_on_boundary(vc123) == CGAL::ZERO;
              }
          }
          RETURN4("Is_degenerate_edge", p1.get_id(), p2.get_id(), p3.get_id(), p4.get_id(), res);
      }
  };

public:

  // PREDICATES
  //-----------

  Comparison_result
  compare_x_2_object() const {
    return Compare_x_2();
  }

  Comparison_result
  compare_y_2_object() const {
    return Compare_y_2();
  }
  
  Comparison_result
  compare_weight_2_object() const {
    return Compare_weight_2();
  }

  Orientation_2
  orientation_2_object() const {
    return Orientation_2();
  }

  Is_hidden_2
  is_hidden_2_object() const {
    return Is_hidden_2();
  }

#ifdef WITH_PSEUDOCIRCLES_PATCH
  Covered_2
  covered_2_object() const {
    return Covered_2();
  }
#endif

  Oriented_side_of_bisector_2
  oriented_side_of_bisector_2_object() const {
    return Oriented_side_of_bisector_2();
  }

  Vertex_conflict_2
  vertex_conflict_2_object() const {
    return Vertex_conflict_2();
  }

  Finite_edge_interior_conflict_2
  finite_edge_interior_conflict_2_object() const {
    return Finite_edge_interior_conflict_2();
  }

  Infinite_edge_interior_conflict_2
  infinite_edge_interior_conflict_2_object() const {
    return Infinite_edge_interior_conflict_2();
  }

  Is_degenerate_edge_2
  is_degenerate_edge_2_object() const {
    return Is_degenerate_edge_2();
  }

};

} //namespace CGAL

#endif // CGAL_DELAUNAY_GRAPH_OF_ELLIPSES_TRAITS_2_H
