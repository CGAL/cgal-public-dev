#include <CGAL/basic.h>
#include <CGAL/Coercion_traits.h>
#include <cassert>

#if defined(CGAL_NEW_COERCION_TRAITS)
struct Foo { };
struct Bar { Bar(const Foo&) { }  };
struct Baz { explicit Baz(const Foo&) { } };

struct X {}; struct Y {}; 
struct UnrelatedType { 
  UnrelatedType(const X&) { }
  UnrelatedType(const Y&) { }
};

namespace std {
  template<>
  struct common_type<Foo, Baz> {
    typedef Baz type;
  };
  template<>
  struct common_type<Baz, Foo> {
    typedef Baz type;
  };

  template<>
  struct common_type<X, Y> {
    typedef UnrelatedType type;
  };
  template<>
  struct common_type<Y, X> {
    typedef UnrelatedType type;
  };
}
#endif

int main(){
  {
    typedef CGAL::Coercion_traits<int,int> CT;
    BOOST_STATIC_ASSERT(( boost::is_same<CT::Type,int>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_implicit_interoperable,CGAL::Tag_true>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_explicit_interoperable,CGAL::Tag_true>::value));
    assert( 5 == CT::Cast()(5));
  }

#if defined(CGAL_NEW_COERCION_TRAITS)
  //test Coercion_traits in some scenarios
  {
    // implicit conversion
    typedef CGAL::Coercion_traits<Foo, Bar> CT;
    BOOST_STATIC_ASSERT(( boost::is_same<CT::Type, Bar>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_implicit_interoperable,CGAL::Tag_true>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_explicit_interoperable,CGAL::Tag_true>::value));
  }
  {
    // explicit with specialized common_type
    typedef CGAL::Coercion_traits<Foo, Baz> CT;
    BOOST_STATIC_ASSERT(( boost::is_same<CT::Type, Baz>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_implicit_interoperable,CGAL::Tag_false>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_explicit_interoperable,CGAL::Tag_true>::value));
  }
  {
    // explicit to UnrelatedType with specialized common_type
    typedef CGAL::Coercion_traits<X, Y> CT;
    BOOST_STATIC_ASSERT(( boost::is_same<CT::Type, UnrelatedType>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_implicit_interoperable,CGAL::Tag_true>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_explicit_interoperable,CGAL::Tag_true>::value));

    BOOST_STATIC_ASSERT(( boost::is_same< decltype( CT::Cast() ( X() ) ),
                                          UnrelatedType >::value));
    BOOST_STATIC_ASSERT(( boost::is_same< decltype(CT::Cast()(Y() ) ),
                                          UnrelatedType >::value));
  }
#endif

  // only run this test with old coercion traits as the new fail to
  // compile and we cannot test this.
#if !defined(CGAL_NEW_COERCION_TRAITS)
  {
    typedef CGAL::Coercion_traits<CGAL::Tag_true,CGAL::Tag_false> CT;
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_implicit_interoperable,CGAL::Tag_false>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_explicit_interoperable,CGAL::Tag_false>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Cast,CGAL::Null_functor>::value));
  }
#endif
  {
    typedef CGAL::Coercion_traits<int, double> CT;
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Are_implicit_interoperable,CGAL::Tag_true>::value));
    BOOST_STATIC_ASSERT(
      ( boost::is_same<CT::Type,double>::value));
  }
}
