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

#ifndef DELAUNAY_GRAPH_OF_ELLIPSES_PROFILED_TRAITS_2_H
#define DELAUNAY_GRAPH_OF_ELLIPSES_PROFILED_TRAITS_2_H


#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <CGAL/Object_cache.h>
#include<CGAL/Timer.h>

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



template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type>
class Profiled_predicate_2_args: public P {
public:
    typedef typename P::result_type result_type;
    static CGAL::Timer timer2;

    result_type operator()(const A1& s1, const A2& s2) {
        timer2.start();
        result_type r = P::operator()(s1, s2);
        timer2.stop();
        return r;
    }
};

template<class P, typename A1, typename A2>
CGAL::Timer Profiled_predicate_2_args<P,A1,A2>::timer2;

template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type>
class Profiled_predicate_3_args: public P {
public:
    typedef typename P::result_type result_type;
    static CGAL::Timer timer3;

    result_type operator()(const A1& s1, const A2& s2, const A3& s3) {
        timer3.start();
        result_type r = P::operator()(s1, s2, s3);
        timer3.stop();
        return r;
    }
};

template<class P, typename A1, typename A2, typename A3>
CGAL::Timer Profiled_predicate_3_args<P,A1,A2,A3>::timer3;

template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type,
                   typename A4 = typename P::argument_type>
class Profiled_predicate_4_args: public P {
public:
    typedef typename P::result_type result_type;
    static CGAL::Timer timer4;

    result_type operator()(const A1& s1, const A2& s2,
                                  const A3& s3, const A4& s4) {
        timer4.start();
        result_type r = P::operator()(s1, s2, s3, s4);
        timer4.stop();
        return r;
    }
};

template<class P, typename A1, typename A2, typename A3, typename A4>
CGAL::Timer Profiled_predicate_4_args<P,A1,A2,A3,A4>::timer4;

template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type,
                   typename A4 = typename P::argument_type,
                   typename A5 = typename P::argument_type>
class Profiled_predicate_5_args: public P {
public:
    typedef typename P::result_type result_type;
    static CGAL::Timer timer5;

    result_type operator()(const A1& s1, const A2& s2, const A3& s3,
                                  const A4& s4, const A5& s5) {
        timer5.start();
        result_type r = P::operator()(s1, s2, s3, s4, s5);
        timer5.stop();
        return r;
    }
};

template<class P, typename A1, typename A2, typename A3, typename A4, typename A5>
CGAL::Timer Profiled_predicate_5_args<P,A1,A2,A3,A4,A5>::timer5;


template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type,
                   typename A4 = typename P::argument_type,
                   typename A5 = typename P::argument_type,
                   typename A6 = typename P::argument_type>
class Profiled_predicate_6_args: public P {
public:
    typedef typename P::result_type result_type;
    static CGAL::Timer timer6;

    result_type operator()(const A1& s1, const A2& s2,
                                  const A3& s3, const A4& s4,
                                  const A5& s5, const A6& s6) {
        timer6.start();
        result_type r = P::operator()(s1, s2, s3, s4, s5, s6);
        timer6.stop();
        return r;
    }
};

template<class P, typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
CGAL::Timer Profiled_predicate_6_args<P,A1,A2,A3,A4,A5,A6>::timer6;


template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type>
class Profiled_predicate_2_3_args: public Profiled_predicate_2_args<P,A1,A2>,
                                 public Profiled_predicate_3_args<P,A1,A2,A3> {
public:
    typedef typename P::result_type result_type;

    inline result_type operator()(const A1& s1, const A2& s2) {
        return Profiled_predicate_2_args<P,A1,A2>::operator ()(s1, s2);
    }

    inline result_type operator()(const A1& s1, const A2& s2, const A3& s3) {
        return Profiled_predicate_3_args<P,A1,A2,A3>::operator ()(s1, s2, s3);
    }
};

template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type,
                   typename A4 = typename P::argument_type>
class Profiled_predicate_3_4_args: public Profiled_predicate_3_args<P,A1,A2,A4>,
                                 public Profiled_predicate_4_args<P,A1,A2,A3,A4> {
public:
    typedef typename P::result_type result_type;

    inline result_type operator()(const A1& s1, const A2& s2, const A4& s3) {
        return Profiled_predicate_3_args<P,A1,A2,A4>::operator ()(s1, s2, s3);
    }

    inline result_type operator()(const A1& s1, const A2& s2,
                                  const A3& s3, const A4& s4) {
        return Profiled_predicate_4_args<P,A1,A2,A3,A4>::operator ()(s1, s2, s3, s4);
    }
};

template< class P, typename A1 = typename P::argument_type, 
                   typename A2 = typename P::argument_type, 
                   typename A3 = typename P::argument_type,
                   typename A4 = typename P::argument_type,
                   typename A5 = typename P::argument_type,
                   typename A6 = typename P::argument_type>
class Profiled_predicate_4_5_6_args: public Profiled_predicate_4_args<P,A1,A2,A3,A6>,
                                   public Profiled_predicate_5_args<P,A1,A2,A3,A4,A6>,
                                   public Profiled_predicate_6_args<P,A1,A2,A3,A4,A5,A6> {
public:
    typedef typename P::result_type result_type;

    inline result_type operator()(const A1& s1, 
                                  const A2& s2,
                                  const A3& s3,
                                  const A6& s4) {
        return Profiled_predicate_4_args<P,A1,A2,A3,A6>::operator ()(s1, s2, s3, s4);
    }

    inline result_type operator()(const A1& s1, const A2& s2,
                                  const A3& s3, const A4& s4,
                                  const A6& s5) {
        return Profiled_predicate_5_args<P,A1,A2,A3,A4,A6>::operator ()(s1, s2, s3, s4, s5);
    }

    inline result_type operator()(const A1& s1, const A2& s2,
                                  const A3& s3, const A4& s4,
                                  const A5& s5, const A6& s6) {
        return Profiled_predicate_6_args<P,A1,A2,A3,A4,A5,A6>::operator ()(s1, s2, s3, s4, s5, s6);
    }
};

template <class DT>
class Delaunay_graph_of_ellipses_profiled_traits_2: public DT {
public:
  typedef DT Base;
  typedef typename Base::Ellipse_traits          Ellipse_traits;
  typedef typename Base::Site_2                  Site_2;

  typedef typename Base::Ellipse_bisector_2 Ellipse_bisector_2;
  typedef Ellipse_bisector_2 Object_2;
  typedef Ellipse_bisector_2 Line_2;
  typedef Ellipse_bisector_2 Ray_2;
  typedef Ellipse_bisector_2 Segment_2;

  typedef typename Base::FT               FT;
  typedef typename Base::RT               RT;
  typedef typename Base::Point_2          Point_2;

//  typedef typename Base::Geom_traits     Geom_traits;

public:

  // PREDICATES
  //-----------

//  Profiled_predicate_3_args<typename Base::Orientation_2>
//  orientation_2_object() const {
//      return Profiled_predicate_3_args<typename Base::Orientation_2>();
//  }

//  Profiled_predicate_2_3_args<typename Base::Is_hidden_2>
//  is_hidden_2_object() const {
//      return Profiled_predicate_2_3_args<typename Base::Is_hidden_2>();
//  }

//  Profiled_predicate_3_args<typename Base::Oriented_side_of_bisector_2, Site_2, Site_2, Point_2>
//  oriented_side_of_bisector_2_object() const {
//    return Profiled_predicate_3_args<typename Base::Oriented_side_of_bisector_2, Site_2, Site_2, Point_2>();
//  }

  Profiled_predicate_3_4_args<typename Base::Vertex_conflict_2>
  vertex_conflict_2_object() const {
    return Profiled_predicate_3_4_args<typename Base::Vertex_conflict_2>();
  }

  Profiled_predicate_4_5_6_args<
            typename Base::Finite_edge_interior_conflict_2,
            Site_2, Site_2, Site_2, Site_2, Site_2, bool>
  finite_edge_interior_conflict_2_object() const {
    return Profiled_predicate_4_5_6_args<
            typename Base::Finite_edge_interior_conflict_2,
            Site_2, Site_2, Site_2, Site_2, Site_2, bool>();
  }

  Profiled_predicate_5_args<
                typename Base::Infinite_edge_interior_conflict_2,
                Site_2, Site_2, Site_2, Site_2, bool>
  infinite_edge_interior_conflict_2_object() const {
    return Profiled_predicate_5_args<
                typename Base::Infinite_edge_interior_conflict_2,
                Site_2, Site_2, Site_2, Site_2, bool>();
  }

  Profiled_predicate_4_args<typename Base::Is_degenerate_edge_2>
  is_degenerate_edge_2_object() const {
    return Profiled_predicate_4_args<typename Base::Is_degenerate_edge_2>();
  }

  static void profile_report() {
      typedef Profiled_predicate_3_4_args<typename Base::Vertex_conflict_2> Vertex_conflict_2;
      typedef Profiled_predicate_4_5_6_args<
                typename Base::Finite_edge_interior_conflict_2,
                Site_2, Site_2, Site_2, Site_2, Site_2, bool> Finite_edge_interior_conflict_2;
      typedef Profiled_predicate_5_args<
              typename Base::Infinite_edge_interior_conflict_2,
              Site_2, Site_2, Site_2, Site_2, bool> Infinite_edge_interior_conflict_2;
      typedef Profiled_predicate_4_args<typename Base::Is_degenerate_edge_2> Is_degenerate_edge_2;

      std::cerr << "Vertex_conflict_2: (3) = " << Vertex_conflict_2::timer3.time() << std::endl;
      std::cerr << "Vertex_conflict_2: (4) = " << Vertex_conflict_2::timer4.time() << std::endl;
      std::cerr << "Finite_edge_interior_conflict_2: (4) = " << Finite_edge_interior_conflict_2::timer4.time() << std::endl;
      std::cerr << "Finite_edge_interior_conflict_2: (5) = " << Finite_edge_interior_conflict_2::timer5.time() << std::endl;
      std::cerr << "Finite_edge_interior_conflict_2: (6) = " << Finite_edge_interior_conflict_2::timer6.time() << std::endl;
      std::cerr << "Infinite_edge_interior_conflict_2: (5) = " << Infinite_edge_interior_conflict_2::timer5.time() << std::endl;
      std::cerr << "Is_degenerate_edge_2: (4) = " << Is_degenerate_edge_2::timer4.time() << std::endl;
  }
};

} //namespace CGAL

#endif // DELAUNAY_GRAPH_OF_ELLIPSES_PROFILED_TRAITS_2_H
