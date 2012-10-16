//    (c) 2011-2012 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED

#ifndef CGAL_POLYNOMIAL_FACTORY_H
#define CGAL_POLYNOMIAL_FACTORY_H

#include<CGAL/Ellipse_2.h>
#include<CGAL/Cache.h>
#include<CGAL/Voronoi_diagram_of_ellipses_2/Generated_code.h>

namespace CGAL {
namespace VORELL {

template<class ET>
class Ellipse_univariate_polynomial_factory {
    typedef typename VORELL::Generated_code<ET> gencode;
    typedef typename ET::upolz_t upolz_t;

    Ellipse_2<ET> e;

    struct Tan_poly_xy_creator {
      typedef Ellipse_2<ET>    argument_type;
      typedef Ellipse_2<ET>    argument1_type;
      typedef upolz_t result_type;
      upolz_t operator()(const Ellipse_2<ET>& a) const {
          return gencode().tan_poly_xy(ELL_PARAM_COEFFS(a),
                                       a.x_center(), a.y_center());
      }
    };

    typedef Cache<Ellipse_2<ET>, upolz_t, Tan_poly_xy_creator> Tan_poly_xy_cache;
    static Tan_poly_xy_cache tan_poly_xy_cache;

public:
    Ellipse_univariate_polynomial_factory(const Ellipse_2<ET>& e1): e(e1) { }
    upolz_t tan_poly_xy() { return tan_poly_xy_cache(e); }
};

template<class ET>
class Ellipse_bivariate_polynomial_factory {
    typedef typename VORELL::Generated_code<ET> gencode;
    typedef typename ET::bpolz_t bpolz_t;

    typedef std::pair<Ellipse_2<ET>, Ellipse_2<ET> > Ellipse_pair;
    Ellipse_pair e;

    struct Polar_creator {
      typedef Ellipse_pair    argument_type;
      typedef Ellipse_pair    argument1_type;
      typedef bpolz_t result_type;
      bpolz_t operator()(const Ellipse_pair& a) const {
          return gencode().polar(ELL_PARAM_COEFFS((a.first)), ELL_PARAM_COEFFS((a.second)));
      }
    };

    struct Tan_poly_cut_creator {
      typedef Ellipse_pair    argument_type;
      typedef Ellipse_pair    argument1_type;
      typedef bpolz_t result_type;
      bpolz_t operator()(const Ellipse_pair& a) const {
          return gencode().tan_poly_cut(ELL_PARAM_COEFFS((a.first)), ELL_PARAM_COEFFS((a.second)));
      }
    };

    typedef Cache<Ellipse_pair, bpolz_t, Polar_creator> Polar_cache;
    static Polar_cache polar_cache;

    typedef Cache<Ellipse_pair, bpolz_t, Tan_poly_cut_creator> Tan_poly_cut_cache;
    static Tan_poly_cut_cache tan_poly_cut_cache;

public:
    Ellipse_bivariate_polynomial_factory(const Ellipse_2<ET>& e1, const Ellipse_2<ET>& e2): e(e1,e2) { }
    bpolz_t polar() { return polar_cache(e); }
    bpolz_t tan_poly_cut() { return tan_poly_cut_cache(e); }
};

template<class ET>
typename Ellipse_univariate_polynomial_factory<ET>::Tan_poly_xy_cache
         Ellipse_univariate_polynomial_factory<ET>::tan_poly_xy_cache;

template<class ET>
typename Ellipse_bivariate_polynomial_factory<ET>::Polar_cache
         Ellipse_bivariate_polynomial_factory<ET>::polar_cache;

template<class ET>
typename Ellipse_bivariate_polynomial_factory<ET>::Tan_poly_cut_cache
         Ellipse_bivariate_polynomial_factory<ET>::tan_poly_cut_cache;


} // namespace
} // namespace CGAL

#endif // CGAL_POLYNOMIAL_FACTORY_H
