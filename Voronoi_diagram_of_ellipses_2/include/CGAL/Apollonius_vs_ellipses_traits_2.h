// based on Apollonius_graph_traits_2 by Menelaos Karavelas

//    (c) 2007-2009 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>

#ifndef CGAL_APOLLONIUS_VS_ELLIPSES_TRAITS_2_H
#define CGAL_APOLLONIUS_VS_ELLIPSES_TRAITS_2_H


#include <CGAL/Delaunay_graph_of_ellipses_traits_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

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

//   NOTE: This is ugly, because it instantiates an object every time.
//   It could be written with templates/double inheritance,
//   as in cached_predicates...
//   OTOH, the compiler might realize this and optimize... :-)

template < class ET >
class Apollonius_vs_ellipses_traits_2 {
public:
  typedef ET                            Ellipse_traits;
  typedef Ellipse_2<ET>                 Site_2;

  typedef CGAL::Ellipse_bisector_2<ET> Ellipse_bisector_2;
  typedef Ellipse_bisector_2 Object_2;
  typedef Ellipse_bisector_2 Line_2;
  typedef Ellipse_bisector_2 Ray_2;
  typedef Ellipse_bisector_2 Segment_2;

  typedef typename ET::QQ               FT;
  typedef typename ET::ZZ               RT;
  typedef typename ET::Kernel::Point_2  Point_2;


private:   
   typedef Apollonius_graph_traits_2<typename ET::Kernel> AP;
   typedef Delaunay_graph_of_ellipses_traits_2<ET> ELL;
   
   static typename AP::Site_2 convert(const Site_2& e) {
       return typename AP::Site_2(e.point(), e.major_axis());
   }
   
public:

  // PREDICATES
  //-----------
   class Compare_x_2 {
   public:
       typedef Comparison_result result_type;
       typedef Site_2            argument_type;

       inline result_type operator()(const Site_2& s1, const Site_2& s2) const 
       {
           result_type r1 = typename ELL::Compare_x_2()(s1, s2);
           result_type r2 = typename AP::Compare_x_2()(convert(s1), convert(s2));
           CGAL_assertion (r1 == r2);
           return r1;
       }
   };
   
   class Compare_y_2 {
   public:
       typedef Comparison_result result_type;
       typedef Site_2            argument_type;

       inline result_type operator()(const Site_2& s1, const Site_2& s2) const 
       {
           result_type r1 = typename ELL::Compare_y_2()(s1, s2);
           result_type r2 = typename AP::Compare_y_2()(convert(s1), convert(s2));
           CGAL_assertion (r1 == r2);
           return r1;
       }
   };
   
   class Compare_weight_2 {
   public:
       typedef Comparison_result result_type;
       typedef Site_2            argument_type;

       inline result_type operator()(const Site_2& s1, const Site_2& s2) const 
       {
           result_type r1 = typename ELL::Compare_weight_2()(s1, s2);
           result_type r2 = typename AP::Compare_weight_2()(convert(s1), convert(s2));
           CGAL_assertion (r1 == r2);
           return r1;
       }
   };

  class Orientation_2 {
  public:
      typedef Orientation result_type;
      typedef Point_2      argument_type;
      
      inline Orientation operator()(const Point_2 &p1, const Point_2& p2, 
                             const Point_2& p3) const {
          result_type r1 = typename ELL::Orientation_2()(p1, p2, p3);
          result_type r2 = typename AP::Orientation_2()(p1, p2, p3);
          CGAL_assertion (r1 == r2);
          return r1;
      }
  };

  class Is_hidden_2 {
  public:
      inline bool operator()(const Site_2 &p, const Site_2 &q) const {
          bool r1 = typename ELL::Is_hidden_2()(p, q);
          bool r2 = typename AP::Is_hidden_2()(convert(p), convert(q));
          CGAL_assertion (r1 == r2);
          return r1;
      }
  };
  
  class Oriented_side_of_bisector_2 {
  public:
      typedef Oriented_side result_type;
      Oriented_side operator()(const Site_2& p1, const Site_2& p2,
			       const Point_2 &p) const {
          result_type r1 = typename ELL::Oriented_side_of_bisector_2()(
                  p1, p2, p);
          result_type r2 = typename AP::Oriented_side_of_bisector_2()(
                  convert(p1), convert(p2), p);
          CGAL_assertion (r1 == r2);
          return r1;
      }
  };

  class Vertex_conflict_2 {
  public:
      typedef Sign    result_type;
      typedef Site_2  argument_type;
      inline Sign operator()(const Site_2& p1, const Site_2& p2,
		             const Site_2& q) const {
          result_type r1 = typename ELL::Vertex_conflict_2()(p1, p2, q);
          result_type r2 = typename AP::Vertex_conflict_2()(
                  convert(p1), convert(p2), convert(q));
          CGAL_assertion (r1 == r2);
          return r1;
      }
      
      inline Sign operator()(const Site_2& p1, const Site_2& p2,
		             const Site_2& p3, const Site_2& q) const {
          result_type r1 = typename ELL::Vertex_conflict_2()(p1, p2, p3, q);
          result_type r2 = typename AP::Vertex_conflict_2()(
                  convert(p1), convert(p2), convert(p3), convert(q));
          CGAL_assertion (r1 == r2);
          return r1;
      }
  };
  
  // adapted from Finite_edge_test_C2 (Apollonius package)
  class Finite_edge_interior_conflict_2 {
  public:
      typedef bool result_type;
      // TODO: add exact version with predicate Order_on_bisector
      bool operator()(const Site_2& p1, const Site_2& p2, const Site_2& p3,
	              const Site_2& p4, const Site_2& q, bool b) const {

          result_type r1 = typename ELL::Finite_edge_interior_conflict_2()(
                  p1, p2, p3, p4, q, b);
          result_type r2 = typename AP::Finite_edge_interior_conflict_2()(
                  convert(p1), convert(p2), convert(p3), convert(p4), 
                  convert(q), b);
          CGAL_assertion (r1 == r2);
          return r1;
      }
      
      // degenerate
      bool operator()(const Site_2& p1, const Site_2& p2, const Site_2& p3,
	              const Site_2& q, bool b) const {
          
          result_type r1 = typename ELL::Finite_edge_interior_conflict_2()(
                  p1, p2, p3, q, b);
          result_type r2 = typename AP::Finite_edge_interior_conflict_2()(
                  convert(p1), convert(p2), convert(p3), convert(q), b);
          CGAL_assertion (r1 == r2);
          return r1;
      }

      bool operator()(const Site_2& p1, const Site_2& p2,
           	      const Site_2& q, bool b) const {
          
          result_type r1 = typename ELL::Finite_edge_interior_conflict_2()(
                  p1, p2, q, b);
          result_type r2 = typename AP::Finite_edge_interior_conflict_2()(
                  convert(p1), convert(p2), convert(q), b);
          CGAL_assertion (r1 == r2);
          return r1;
      }
  };

  // adapted from Infinite_edge_test_C2 (Apollonius package)
  class Infinite_edge_interior_conflict_2 {
      
  public:
      typedef bool result_type;
  
      bool operator()(const Site_2& p2, const Site_2& p3,
	              const Site_2& p4, const Site_2& q, bool b) const {
          result_type r1 = typename ELL::Infinite_edge_interior_conflict_2()(
                  p2, p3, p4, q, b);
          result_type r2 = typename AP::Infinite_edge_interior_conflict_2()(
                  convert(p2), convert(p3), convert(p4), convert(q), b);
          CGAL_assertion (r1 == r2);
          return r1;
      }
      
  };
  
  class Is_degenerate_edge_2 {
  public:
      typedef bool result_type;
      bool operator()(const Site_2& p1, const Site_2& p2,
	    	       const Site_2& p3, const Site_2& p4) const {
      
          result_type r1 = typename ELL::Is_degenerate_edge_2()(
                  p1, p2, p3, p4);
          result_type r2 = typename AP::Is_degenerate_edge_2()(
                  convert(p1), convert(p2), convert(p3), convert(p4));
          CGAL_assertion (r1 == r2);
          return r1;
      }
  };

public:

  // PREDICATES
  //-----------

  Orientation_2
  orientation_2_object() const {
    return Orientation_2();
  }

  Is_hidden_2
  is_hidden_2_object() const {
    return Is_hidden_2();
  }

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

#endif // CGAL_APOLLONIUS_VS_ELLIPSES_TRAITS_2_H
