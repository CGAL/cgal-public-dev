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

#ifndef DELAUNAY_GRAPH_OF_ELLIPSES_CACHED_TRAITS_2_H
#define DELAUNAY_GRAPH_OF_ELLIPSES_CACHED_TRAITS_2_H


#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <CGAL/Object_cache.h>

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
class Cached_predicate_2_args: public P {
public:
    typedef typename P::result_type result_type;

private:
    typedef ::std::pair<A1, A2> Key2;
    typedef Value_cache<Key2, result_type> Cache2;
    Cache2 cache2;

public:
    result_type operator()(const A1& s1, const A2& s2) {
        Key2 k = ::std::make_pair(s1, s2);
        typename Cache2::Found_value fv = cache2.read(k);
        if (fv.first) return fv.second;
        result_type r = P::operator()(s1, s2);
        cache2.write(k, r);
        return r;
    }
};


template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type>
class Cached_predicate_3_args: public P {
public:
    typedef typename P::result_type result_type;

private:
    typedef ::boost::tuple<A1, A2, A3> Key3;
    typedef Value_cache<Key3, result_type> Cache3;
    Cache3 cache3;

public:
    result_type operator()(const A1& s1, const A2& s2, const A3& s3) {
        Key3 k = ::boost::make_tuple(s1, s2, s3);
        typename Cache3::Found_value fv = cache3.read(k);
        if (fv.first) return fv.second;
        result_type r = P::operator()(s1, s2, s3);
        cache3.write(k, r);
        return r;
    }
};

template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type,
                   typename A4 = typename P::argument_type>
class Cached_predicate_4_args: public P {
public:
    typedef typename P::result_type result_type;

private:
    typedef ::boost::tuple<A1, A2, A3, A4> Key4;
    typedef Value_cache<Key4, result_type> Cache4;
    Cache4 cache4;

public:
    result_type operator()(const A1& s1, const A2& s2,
                                  const A3& s3, const A4& s4) {
        Key4 k = ::boost::make_tuple(s1, s2, s3, s4);
        typename Cache4::Found_value fv = cache4.read(k);
        if (fv.first) return fv.second;
        result_type r = P::operator()(s1, s2, s3, s4);
        cache4.write(k, r);
        return r;
    }
};

template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type,
                   typename A4 = typename P::argument_type,
                   typename A5 = typename P::argument_type>
class Cached_predicate_5_args: public P {
public:
    typedef typename P::result_type result_type;

private:
    typedef ::boost::tuple<A1, A2, A3, A4, A5> Key5;
    typedef Value_cache<Key5, result_type> Cache5;
    Cache5 cache5;

public:
    result_type operator()(const A1& s1, const A2& s2, const A3& s3,
                                  const A4& s4, const A5& s5) {
        Key5 k = ::boost::make_tuple(s1, s2, s3, s4, s5);
        typename Cache5::Found_value fv = cache5.read(k);
        if (fv.first) return fv.second;
        result_type r = P::operator()(s1, s2, s3, s4, s5);
        cache5.write(k, r);
        return r;
    }
};

template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type,
                   typename A4 = typename P::argument_type,
                   typename A5 = typename P::argument_type,
                   typename A6 = typename P::argument_type>
class Cached_predicate_6_args: public P {
public:
    typedef typename P::result_type result_type;

private:
    typedef ::boost::tuple<A1, A2, A3, A4, A5, A6> Key6;
    typedef Value_cache<Key6, result_type> Cache6;
    Cache6 cache6;

public:

    result_type operator()(const A1& s1, const A2& s2,
                                  const A3& s3, const A4& s4,
                                  const A5& s5, const A6& s6) {
        Key6 k = ::boost::make_tuple(s1, s2, s3, s4, s5, s6);
        typename Cache6::Found_value fv = cache6.read(k);
        if (fv.first) return fv.second;
        result_type r = P::operator()(s1, s2, s3, s4, s5, s6);
        cache6.write(k, r);
        return r;
    }
};

template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type>
class Cached_predicate_2_3_args: public Cached_predicate_2_args<P,A1,A2>,
                                 public Cached_predicate_3_args<P,A1,A2,A3> {
public:
    typedef typename P::result_type result_type;

    inline result_type operator()(const A1& s1, const A2& s2) {
        return Cached_predicate_2_args<P,A1,A2>::operator ()(s1, s2);
    }

    inline result_type operator()(const A1& s1, const A2& s2, const A3& s3) {
        return Cached_predicate_3_args<P,A1,A2,A3>::operator ()(s1, s2, s3);
    }
};

template< class P, typename A1 = typename P::argument_type,
                   typename A2 = typename P::argument_type,
                   typename A3 = typename P::argument_type,
                   typename A4 = typename P::argument_type>
class Cached_predicate_3_4_args: public Cached_predicate_3_args<P,A1,A2,A4>,
                                 public Cached_predicate_4_args<P,A1,A2,A3,A4> {
public:
    typedef typename P::result_type result_type;

    inline result_type operator()(const A1& s1, const A2& s2, const A4& s3) {
        return Cached_predicate_3_args<P,A1,A2,A4>::operator ()(s1, s2, s3);
    }

    inline result_type operator()(const A1& s1, const A2& s2,
                                  const A3& s3, const A4& s4) {
        return Cached_predicate_4_args<P,A1,A2,A3,A4>::operator ()(s1, s2, s3, s4);
    }
};

template< class P, typename A1 = typename P::argument_type, 
                   typename A2 = typename P::argument_type, 
                   typename A3 = typename P::argument_type,
                   typename A4 = typename P::argument_type,
                   typename A5 = typename P::argument_type,
                   typename A6 = typename P::argument_type>
class Cached_predicate_4_5_6_args: public Cached_predicate_4_args<P,A1,A2,A3,A6>,
                                   public Cached_predicate_5_args<P,A1,A2,A3,A4,A6>,
                                   public Cached_predicate_6_args<P,A1,A2,A3,A4,A5,A6> {
public:
    typedef typename P::result_type result_type;

    inline result_type operator()(const A1& s1, 
                                  const A2& s2,
                                  const A3& s3,
                                  const A6& s4) {
        return Cached_predicate_4_args<P,A1,A2,A3,A6>::operator ()(s1, s2, s3, s4);
    }

    inline result_type operator()(const A1& s1, const A2& s2,
                                  const A3& s3, const A4& s4,
                                  const A6& s5) {
        return Cached_predicate_5_args<P,A1,A2,A3,A4,A6>::operator ()(s1, s2, s3, s4, s5);
    }

    inline result_type operator()(const A1& s1, const A2& s2,
                                  const A3& s3, const A4& s4,
                                  const A5& s5, const A6& s6) {
        return Cached_predicate_6_args<P,A1,A2,A3,A4,A5,A6>::operator ()(s1, s2, s3, s4, s5, s6);
    }
};

template <class DT>
class Delaunay_graph_of_ellipses_cached_traits_2: public DT {
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

//  Cached_predicate_3_args<typename Base::Orientation_2>
//  orientation_2_object() const {
//      return Cached_predicate_3_args<typename Base::Orientation_2>();
//  }

//  Cached_predicate_2_3_args<typename Base::Is_hidden_2>
//  is_hidden_2_object() const {
//      return Cached_predicate_2_3_args<typename Base::Is_hidden_2>();
//  }

//  Cached_predicate_3_args<typename Base::Oriented_side_of_bisector_2, Site_2, Site_2, Point_2>
//  oriented_side_of_bisector_2_object() const {
//    return Cached_predicate_3_args<typename Base::Oriented_side_of_bisector_2, Site_2, Site_2, Point_2>();
//  }

  Cached_predicate_3_4_args<typename Base::Vertex_conflict_2>
  vertex_conflict_2_object() const {
    return Cached_predicate_3_4_args<typename Base::Vertex_conflict_2>();
  }

  Cached_predicate_4_5_6_args<
            typename Base::Finite_edge_interior_conflict_2,
            Site_2, Site_2, Site_2, Site_2, Site_2, bool>
  finite_edge_interior_conflict_2_object() const {
    return Cached_predicate_4_5_6_args<
            typename Base::Finite_edge_interior_conflict_2,
            Site_2, Site_2, Site_2, Site_2, Site_2, bool>();
  }

  Cached_predicate_5_args<
                typename Base::Infinite_edge_interior_conflict_2,
                Site_2, Site_2, Site_2, Site_2, bool>
  infinite_edge_interior_conflict_2_object() const {
    return Cached_predicate_5_args<
                typename Base::Infinite_edge_interior_conflict_2,
                Site_2, Site_2, Site_2, Site_2, bool>();
  }

  Cached_predicate_4_args<typename Base::Is_degenerate_edge_2>
  is_degenerate_edge_2_object() const {
    return Cached_predicate_4_args<typename Base::Is_degenerate_edge_2>();
  }
};

} //namespace CGAL

#endif // DELAUNAY_GRAPH_OF_ELLIPSES_CACHED_TRAITS_2_H
