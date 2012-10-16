//    (c) 2007-2009 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED

#ifndef CGAL_VORELL_H
#define CGAL_VORELL_H

#include<CGAL/basic.h>
#include<CGAL/Fraction_traits.h>
#include<CGAL/Coercion_traits.h>

namespace CGAL {

namespace VORELL {

template<class POLYQ>
typename CGAL::Fraction_traits<POLYQ>::Numerator_type 
primpart(const POLYQ& polyq) 
{
    typedef typename CGAL::Fraction_traits<POLYQ>::Numerator_type POLYZ;
    typedef CGAL::Polynomial_traits_d<POLYZ> PT;
    POLYZ polyz;
    typename Fraction_traits<POLYQ>::Denominator_type denom;
    typename Fraction_traits<POLYQ>::Decompose()(polyq, polyz, denom);

    Sign s1 = CGAL::sign(typename PT::Innermost_leading_coefficient()(polyz));
    polyz = typename PT::Canonicalize()(polyz);
    Sign s2 = CGAL::sign(typename PT::Innermost_leading_coefficient()(polyz));
    // keep sign!
    if (s1 != s2) return -polyz;
    return polyz;
}

template <class IPOLY, class INT> 
typename CGAL::Interval_traits<INT>::Bound cauchy_bound(const IPOLY &p)
{
    typedef CGAL::Polynomial_traits_d<IPOLY> PT;
    typedef typename CGAL::Interval_traits<INT>::Bound BT;
    typename PT::Innermost_leading_coefficient lcoeff;
    typename PT::Get_coefficient coeff;
    CGAL_assertion ( !CGAL::zero_in(lcoeff(p)) );
    int d = typename PT::Degree()(p);
    if (d == 0) return upper(abs(lcoeff(p)));
    BT m = upper(abs(coeff(p, d-1)));
    for (int i = d-2; i >= 0; i--) if (upper(abs(coeff(p,i))) > m) m = upper(abs(coeff(p,i)));
    m += upper(abs(lcoeff(p)));
    INT a = INT(m) / INT(upper(abs(lcoeff(p))));
    return upper(a);
}

}

} //namespace CGAL
#endif
