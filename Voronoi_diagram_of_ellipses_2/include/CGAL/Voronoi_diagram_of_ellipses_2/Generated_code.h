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

#ifndef CGAL_GENERATED_CODE_H
#define CGAL_GENERATED_CODE_H

#include<CGAL/basic.h>
#include<CGAL/vorell.h>
#include<CGAL/Ellipse_2.h>

namespace CGAL {

namespace VORELL {

template<class ET>
class Generated_code { // TODO: maybe rename to/create Polynomial_factory, or cache...
public:
    typedef typename ET::QQ QQ;
    typedef typename ET::upolz_t upolz_t;
    typedef typename ET::upolq_t upolq_t;
    typedef typename ET::bpolz_t bpolz_t;
    typedef typename ET::bpolq_t bpolq_t;
    typedef typename ET::tpolz_t tpolz_t;
    typedef typename ET::tpolq_t tpolq_t;

private:        
    // side_of_bisector
    QQ eval_R(const Ellipse_2<ET>& e1, const QQ& x, const QQ& y) const;

public:
    // side_of_bisector
    // Delta[v1,v2](s), deg_s = 4
    upolz_t distance_poly(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc, const QQ& v1, const QQ& v2) const;

    // distance_from_bitangent
    // deg = 4
    upolz_t btan_poly(const QQ& ai, const QQ& bi, const QQ& wi, const QQ& xci, const QQ& yci,
             const QQ& aj, const QQ& bj, const QQ& wj, const QQ& xcj, const QQ& ycj) const;
    // deg = 2
    upolq_t tan_poly_xy_q(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc, 
            const QQ& x, const QQ& y) const;
    upolz_t tan_poly_xy(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc, 
            const QQ& x, const QQ& y) const;
    // discriminant of the polar (tangents) from 1 to 2
    // deg = 4 (due to being disc. of quadratic polynomial)
    // when < 0 then there are no tangents from 1 to 2 => pt1 is inside 2
    upolz_t disc_polar_poly(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1,
                            const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const;
    
    // deg = 2
    // intersection of line and ellipse
/*    upolz_t line_ellipse_intersection(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc, 
            const QQ& x1, const QQ& y1, const QQ& x2, const QQ& y2) const;*/

    // deg = 6+6, 1=t, 2=r, B(t,r) -> B(x,y) = y(x)...
    bpolz_t bisector(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                     const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const;
    // deg = 2+2, 1=t, 2=r, P(t,r) -> P(x,y)
    // tangent to Er from (x(t),y(t))
    bpolq_t polar_q(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                  const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const;
    bpolz_t polar(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                  const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const;
    // deg = 2+2
    bpolz_t tan_poly_cut(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                         const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const;
    
    // deg = 4+2, 1=t, 2=r, P(t^4,r^2) -> P(x^2,y^2) = y(x)
    bpolz_t normal_cut(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                       const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const;
    // deg = 12
    upolz_t parnorm(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                    const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const;
    // deg = 12
    upolq_t parnorm_testing(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                            const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const;
    

    // intersection of two ellipses
    // deg = 4
    upolz_t inter_poly(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                       const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const;

    // NORMALS
    void ell_normal(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc, 
            upolq_t& poly_x, upolq_t& poly_y, upolq_t& poly_c);
    // deg = 4+4+4 P(t,r,s) -> P(x,y,z) = z(y(x))...
    tpolz_t ell_normals(const QQ& a1_, const QQ& b1_, const QQ& w1, const QQ& xc1, const QQ& xy1, 
		    const QQ& a2_, const QQ& b2_, const QQ& w2, const QQ& xc2, const QQ& xy2,
		    const QQ& a3_, const QQ& b3_, const QQ& w3, const QQ& xc3, const QQ& xy3);
    

    // deg = 2+2 D(t,r) -> D(x,y) = y(x)
    bpolz_t denom(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                  const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const;
    
    upolq_t numerxt(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1) const;
    upolq_t numeryt(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1) const;
    upolq_t denomt(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1) const;

    // MEDIAL AXIS
    // deg = 2+1 M(t,x) -> M(x,y) = y(x)
    bpolz_t medial_x(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc) const;
};

template<class ET>
typename ET::QQ Generated_code<ET>::eval_R(const Ellipse_2<ET>& e1, const QQ& x, const QQ& y) const 
{
    return e1.invariant_J2()*(x*x+y*y) + 2*(e1.algebraic_coefficient(0)*e1.algebraic_coefficient(4)-
            e1.algebraic_coefficient(1)*e1.algebraic_coefficient(3))*y + 
            2*(e1.algebraic_coefficient(2)*e1.algebraic_coefficient(3)-
            e1.algebraic_coefficient(1)*e1.algebraic_coefficient(4))*x + 
            e1.invariant_J1()*e1.algebraic_coefficient(5) - e1.algebraic_coefficient(3)*e1.algebraic_coefficient(3) - 
            e1.algebraic_coefficient(4)*e1.algebraic_coefficient(4);
}

#if 0
// deg(P) = 4
template<class ET>
typename ET::upolz_t Generated_code<ET>::distance_poly(const Ellipse_2<ET>& e1, const QQ& v1, const QQ& v2) const
{
    std::vector<QQ> poly;

    poly.reserve(5);
    QQ E = e1.evaluate_equation(v1, v2);
    QQ R = eval_R(e1, v1, v2);
    QQ J1 = e1.invariant_J1();
    QQ J2 = e1.invariant_J2();
    QQ t2 = J2 * J2;
    QQ t6 = t2 * E;
    QQ t4 = J2 * E * J1;
    QQ t1 = E * E;
    QQ t3 = t6;
    QQ t5 = R * R;
    poly.push_back(  t1 * (4 * t3 + t5) );
    QQ t14 = t1 * t2;
    poly.push_back( 2 * (-6*t14-t5*E) * J1 + 2*(2*t5 + 9*t3 - t1*J2)*R );
    QQ t27 = J1 * J1;
    QQ t34 = t2 * J1;
    poly.push_back(  t14 - 27 * t2 * t2 + (-18 * t6 - 12 * t5) * J2 + 
           (4 * t4 - 18 * t34) * R + (12 * t3 + t5) * t27 );
    poly.push_back( -2 * J2 * (-9 * t34 + t4 - 6 * J2 * R + 
           (2 * J2 * J1 + R) * t27) );
    poly.push_back( -t2 * (-t27 + 4 * J2) );

    return CGAL::VORELL::primpart(typename ET::upolq_t(poly.begin(), poly.end()));
}
#endif

// deg(P) = 4
template<class ET>
typename ET::upolz_t Generated_code<ET>::distance_poly(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc,
                                                       const QQ& v1, const QQ& v2) const
{
    std::vector<QQ> poly;
    poly.reserve(5);

    QQ t13, t54, t52, t51, t18, t14, t15, t35, c1 , t10, t11, t8 , t36, t43,
       t27, t33, t34, t7 , t9 , t39, t40, t42, c0, t32, t19, F2 , t6 , t48, c3,
       t47, t46, c2 , t4 , t5 , t45, t3 , t44, t41, t2;

    t13 = w*w;
    t54 = 4*(1-t13)*w;
    t52 = (-xc+v1)*t54;
    t51 = (-yc+v2)*t54;
    t18 = (1+t13)*(1+t13);
    t14 = b*b;
    t15 = a*a;
    t35 = t15+t14;
    c1 = t18*t35;
    t10 = v1*v1;
    t11 = t13*t13;
    t8 = xc*xc;
    t36 = -t10-t8;
    t43 = v1*xc;
    t27 = -2*t43-t36;
    t33 = t15+t36;
    t34 = 2*t13;
    t7 = yc*yc;
    t9 = v2*v2;
    t39 = -t9-t7;
    t40 = 2*t11+2;
    t42 = v2*yc;
    c0 = (t39*t11-t52*yc+(-2*t10+4*t43-2*t8-t39)*t34+((-4*t13+t40)*yc
          +t52)*v2+t39)*t15+(t33*t11+t51*xc+(t40*xc-t51)*v1+(-2*t9-2*t7+4*t42+t15+
          t27)*t34+t33)*t14;
    t32 = c1*c0;
    t19 = t18*t18;
    F2 = t15*t14*t19;
    t6 = F2*F2;
    t48 = t6*t32;
    c3 = t15*t15*t14*t14*t18*t19;
    t47 = c3*c0;
    t46 = c1*c3;
    c2 = 2*t42-t27+t35+t39;
    t4 = c2*c2;
    t5 = F2*t6;
    t45 = t5*t4;
    t3 = c1*c1;
    t44 = t6*t3;
    t41 = 9*F2;
    t2 = c0*c0;
    poly.push_back( (t6*t4-4*t47)*t2 );
    poly.push_back( 2*t4*t48-12*t2*t46+2*(-2*t45+t6*t2+t41*t47)*c2 );
    poly.push_back( -12*t45+(t4*t3+t2+4*c2*t32)*t6+(-12*t3*c0-27*c3+18*(c0+c2*c1)*F2)*c3 );
    poly.push_back( 2*t48+2*(-6*t5+t44)*c2+2*(-2*t3+t41)*t46 );
    poly.push_back( t44-4*t5 );

    return CGAL::VORELL::primpart(typename ET::upolq_t(poly.begin(), poly.end()));
}

// deg(P) = 4
template<class ET>
typename ET::upolz_t Generated_code<ET>::btan_poly(const QQ& ai, const QQ& bi, const QQ& wi, const QQ& xci, const QQ& yci,
                 const QQ& aj, const QQ& bj, const QQ& wj, const QQ& xcj, const QQ& ycj) const
{
    std::vector<QQ> poly;
    
    poly.reserve(5);
    QQ t110, t112, t105, t107, t235, t263, t284, t312, t266, t281,
    t106, t232, t254, t288, t125, t306, t290, t315, t109, t238,
    t291, t311, t301, t267, t229, t218, t219, t153, t230, t101,
    t100, t225, t154, t140, t104, t296, t119, t111, t234, t293,
    t123, t314, t255, t313, t108, t236, t308, t113, t220, t283,
    t221, t195, t307, t300, t299, t227, t297, t294, t201, t203,
    t228, t210, t80, t276, t275, t273, t205, t181, t241, t180,
    t237, t212, t216, t268, t259, t253, t252, t250, t249, t244,
    t243, t240, t239, t226, t224, t223, t222, t214, t207, t206,
    t202, t200, t160, t152, t142, t141, t137, t128, t127, t124,
    t121, t118;

    t110 = wj*wj;
    t112 = t110*t110;
    t105 = wi*wi;
    t107 = wi*t105;
    t235 = ycj-yci;
    t263 = 8*t235;
    t284 = t107*t263;
    t312 = -2*xcj+2*xci;
    t266 = 4*t235;
    t281 = t107*t266;
    t106 = t105*t105;
    t232 = xcj-xci;
    t254 = t106*t232;
    t288 = wi*t266;
    t125 = -t288-t281+2*t254+t312;
    t306 = 4*xci-4*xcj;
    t290 = wi*t263;
    t315 = ((t290+t284-4*t254-t306)*t110-t125*t112-t125)*ai;
    t109 = aj*aj;
    t238 = -2*t106-2;
    t291 = t106+1;
    t311 = t232*t235;
    t301 = -t107+wi;
    t267 = 6*t235;
    t229 = xci*xcj;
    t218 = -2*t229;
    t219 = -xci*xci-xcj*xcj;
    t153 = t218-t219;
    t230 = ycj*yci;
    t101 = ycj*ycj;
    t100 = yci*yci;
    t225 = -t101-t100;
    t154 = -2*t230-t225;
    t140 = -t153+t154;
    t104 = bj*bj;
    t296 = t104-t109;
    t119 = 2*t301*(-t140+t296);
    t111 = t110*wj;
    t234 = -t111+wj;
    t293 = t296*t234;
    t123 = -t311-2*t293;
    t314 = (-(-t119-t232*(t105*t267-t291*t235))*t112+(t232*t267+12*t293)*
t105-(t301*(-12*t109+12*t104+4*t140)+(-12*t105-t238)*t311)*t110+t123*
t106+t119+t123)*ai;
    t255 = t112+1;
    t313 = t255*t291;
    t108 = ai*ai;
    t236 = t108*bi;
    t308 = t312*t236;
    t113 = bi*bi;
    t220 = t113*t110;
    t283 = -2*t291*t220;
    t221 = t113*t105;
    t195 = t112*t221;
    t307 = -2*t195+t283;
    t300 = t290-t284;
    t299 = t288-t281;
    t227 = t113*t107;
    t297 = 8*t234*t227;
    t294 = t301*t232;
    t201 = t104*t220;
    t203 = xci*t220;
    t228 = t113*xci;
    t210 = xcj*t228;
    t80 = 2*t210;
    t276 = t112*t80+4*xcj*t203+t80+4*t201;
    t275 = -16*t294;
    t273 = -8*t294;
    t205 = xcj*t227;
    t181 = t221*t230;
    t241 = 4*t104;
    t180 = t105*t210;
    t237 = t113*wi;
    t212 = xcj*t237;
    t216 = yci*t237;
    t268 = -4*xci*t216+8*t181-4*t180+(t306*t227+4*t212)*yci+(-4*
t212+4*t205+(4*wi-4*t107)*t228)*ycj+(t241-2*t219)*t221;
    t259 = 4*t106;
    t253 = 2*t104-t113;
    t252 = t113*t109*t313;
    t250 = ((-2-4*t110)*t221+t307)*t108;
    t249 = 4*bi;
    t244 = -4*t104;
    t243 = 8*t109;
    t240 = 2*t109;
    t239 = t259+4;
    t226 = t109-2*t104;
    t224 = 2*t219;
    t223 = -t104+t240;
    t222 = 2*t225;
    t214 = 4*t105;
    t207 = t111*t237;
    t206 = t110*t236;
    t202 = t105*t220;
    t200 = t109*t221;
    t160 = t222+t253;
    t152 = -4*t230-t222;
    t142 = t152-t153;
    t141 = 4*t229+t154+t224;
    t137 = t255*(-yci-t106*ycj)*t108;
    t128 = 2*(t142+t226)*t105-t219;
    t127 = 4*t255*(t106*yci+ycj)*t236+4*(wi+t107)*(t306*t206+t308*t112+
t308)+8*(-t106+1)*t206*t235;
    t124 = (t141+t223)*t214+t160;
    t121 = (-t108+t219)*t313;
    t118 = -16*t105*t201+t104*t297+t268*t112+t268+t252+8*t104*t207+(20*
t200+16*t181-8*t180)*t110+(t243-8*t104)*wj*t237-2*t200+t276*t106+t276+
t250-t225*(-4*t221-4*t195-8*t202)+(-8*t207-t297+t307)*t109-t219*(4*
t202+t283)-8*t232*(wi*ycj*t220-t110*t216)+8*(-t107*t203+t110*t205)*t235;
    
    poly.push_back( (t315+t121)*t113+t118 );
    poly.push_back( t127+(t137+t314)*t249 );
    poly.push_back( 2*(t218+t153*t106+t128+(-t219*t106+t299*xcj+(t238*xcj-t299)*
xci+t128)*t112+(t240+t244+t300*xcj+2*(-t219+t226)*t106+(-t239*xcj-t300)*xci+(
-5*t109+t241+t142)*t214-t224)*t110+t301*(8*t293+t232*t266))*t113+2*(4*
t230+(-t152+t253)*t106+t124+(t243-4*t100+t244-4*t101+t275*yci+(t223+t225)*
t259+((8+8*t106)*yci-t275)*ycj+8*(5*t104-4*t109+t141)*t105)*t110+(
t273*yci+t160*t106+(t239*yci-t273)*ycj+t124)*t112+t301*(16*t293+t232*t263))*
t108+2*t250-2*t252 );
    poly.push_back( (t137-t314)*t249+t127 );
    poly.push_back( (-t315+t121)*t113+t118 );
    
    return CGAL::VORELL::primpart(typename ET::upolq_t(poly.begin(), poly.end()));
}

// deg(P) = 2
template<class ET>
typename ET::upolq_t Generated_code<ET>::tan_poly_xy_q(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc, 
                    const QQ& x, const QQ& y) const
{
    std::vector<QQ> poly;
    
    poly.reserve(3);
    QQ t6, t5, t4, t3, t8, t9;

    t6 = y-yc;
    t5 = xc-x;
    t3 = w*w;
    t4 = 2*w;
    t9 = (-t5*t3-t6*t4+t5)*b;
    t8 = (t3+1)*b*a;
    poly.push_back( 2*t9+2*t8 );
    poly.push_back( 4*(t6*t3-t5*t4-t6)*a );
    poly.push_back( -2*t9+2*t8 );
    
    return typename ET::upolq_t(poly.begin(), poly.end());
}

// deg(P) = 2
template<class ET>
typename ET::upolz_t Generated_code<ET>::tan_poly_xy(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc, 
                    const QQ& x, const QQ& y) const
{
    return CGAL::VORELL::primpart(tan_poly_xy_q(a, b, w, xc, yc, x, y));
}

// deg(P) = 4
template<class ET>
typename ET::upolz_t Generated_code<ET>::disc_polar_poly(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1,
                        const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const
{
    std::vector<QQ> poly;
    poly.reserve(5);

    QQ t300, t181, t179, t182, t184, t344, t178, t177, t374, t204, t176, t313, t272,
       t302, t194, t342, t196, t180, t304, t311, t315, t359, t399, t297, t398, t397,
       t362, t206, t199, t336, t357, t396, t394, t299, t263, t306, t265, t388, t294,
       t257, t289, t387, t386, t290, t385, t295, t384, t379, t262, t377, t248, t249,
       t291, t370, t230, t231, t366, t253, t254, t296, t280, t365, t232, t363, t185,
       t288, t274, t356, t171, t251, t355, t351, t309, t267, t310, t270, t286, t287,
       t350, t349, t348, t334, t229, t332, t314, t312, t273, t208, t191, t190, t183;

    t300 = yc2-yc1;
    t181 = w1*w1;
    t179 = t181*t181;
    t182 = b2*b2;
    t184 = a2*a2;
    t344 = 4*t300;
    t178 = w2*w2;
    t177 = w2*t178;
    t374 = w2-t177;
    t204 = t374*t344;
    t176 = t178*t178;
    t313 = 2+2*t176;
    t272 = -4*t178+t313;
    t302 = xc1-xc2;
    t194 = t204-t302*t272;
    t342 = 8*t302;
    t196 = t178*t342+t204;
    t180 = w1*t181;
    t304 = w1+t180;
    t311 = -4*t176-4;
    t315 = 8*t178;
    t359 = t300*t178;
    t399 = 16*((t196*t179-t196+t304*(-t300*(t315+t311)+t374*t342))*t184-(
              t194*t179-t204-16*t304*t359+t302*(t272+8*t374*t304))*t182)*a1;
    t297 = t182*t177;
    t398 = t297+t184*w2;
    t397 = t182*w2+t184*t177;
    t362 = 2*t302;
    t206 = t374*t362;
    t199 = t178*t344-t206;
    t336 = t176+1;
    t357 = t336*t300;
    t396 = -((t206+t357-2*t359)*t179-t357+(2*t300-t304*t342)*t178-t374*(
              t362+t304*t344))*t184-(t199*t179-t199+t304*t194)*t182;
    t394 = t180-w1;
    t299 = yc2*t184;
    t263 = yc1*t299;
    t306 = xc1*xc2;
    t265 = t182*t306;
    t388 = -2*t263-2*t265;
    t294 = t182*t181;
    t257 = t177*t294;
    t289 = t184*t181;
    t387 = w2*t289+t257;
    t386 = w2*t294+t177*t289;
    t290 = t184*t179;
    t385 = t184+2*t289+t290;
    t295 = t182*t179;
    t384 = 2*t294+t182+t295;
    t379 = t374*t394;
    t262 = yc2*t289;
    t377 = -8*t387*yc1+8*w2*t262+8*yc2*t257+(8*yc1-8*yc2)*t386;
    t248 = t178*t289;
    t249 = t178*t290;
    t291 = t184*t178;
    t370 = -4*t248-2*t291-2*t249;
    t230 = yc1*t262;
    t231 = t181*t265;
    t366 = -2*t182*t289-4*t231-4*t230+t388*t179+t388;
    t253 = t178*t294;
    t254 = t178*t295;
    t296 = t182*t178;
    t280 = -2*t296;
    t365 = -4*t253+8*t248+4*t291+t280+4*t249-2*t254+t384*t176+t384;
    t232 = yc1*yc2*t296;
    t363 = -8*t291*t306-8*t232+4*(t263+t265)*t178+(-t398*yc1+w2*t299+
              yc2*t297-t300*t397)*(4*xc1-4*xc2);
    t185 = a1*a1;
    t288 = t185*t184;
    t274 = 4*t288;
    t356 = t185*t280+t178*t274;
    t171 = t182*t185;
    t251 = t181*t171;
    t355 = t181*t274-2*t251;
    t351 = (-2*t398+2*t397)*a1;
    t309 = t182*a1;
    t267 = t180*t309;
    t310 = a1*t184;
    t270 = t180*t310;
    t286 = w1*t309;
    t287 = w1*t310;
    t350 = -2*t286+2*t267-2*t270+2*t287;
    t349 = -8*t379;
    t348 = 16*t379;
    t334 = t336*(-t179-1)*t171;
    t229 = -16*t248;
    t332 = -16*t181*t232+t229*t306+t363*t179+(t231+t230)*t315+t370*t182+
              t366*t176+t363+t366+(yc1*yc1+yc2*yc2)*(4*t254+4*t296+8*t253+t370+t385*
              t176+t385)+(t365*xc2-t377)*xc2+(t365*xc1+t377)*xc1;
    t314 = -4-4*t179;
    t312 = 2+2*t179;
    t273 = -4*t181+t312;
    t208 = (-t176*t179-t179-t336)*t182;
    t191 = -12*(t267+t287)*t178+12*(t270+t286)*t178+t350+t351+t351*t179+
              t350*t176+(-12*t386+12*t387)*a1;
    t190 = 320*t178*t251+16*t332+16*t184*t208+16*t185*t229-16*t334+
              16*t356*t179+16*t356+16*t355+16*t355*t176-16*(t288-t171)*t394*(8*
              t177-8*w2);
    t183 = b1*b1;
    poly.push_back( t190-t399 );
    poly.push_back( 64*b1*(t191-t396) );
    poly.push_back( 32*(t208+(t311*t181+(16*t181+t314)*t178+t349)*t185+((40*
              t181+t314)*t178+t273*t176+t273-t348)*t183)*t184+32*((t313*t181+(-20*t181+
              t312)*t178-t349)*t185+((8*t176+8)*t181+(8*t179-32*t181+8)*t178+t348)*
              t183)*t182+32*t332+32*t334 );
    poly.push_back( -64*b1*(t191+t396) );
    poly.push_back( t190+t399 );

    return CGAL::VORELL::primpart(typename ET::upolq_t(poly.begin(), poly.end()));
}

// deg = 2
// intersection of line and ellipse
/*template<class ET>
typename ET::upolz_t Generated_code<ET>::line_ellipse_intersection(
        const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc, 
        const QQ& x1, const QQ& y1, const QQ& x2, const QQ& y2) const
{
    std::vector<QQ> poly;
    
    poly.reserve(3);
    QQ t12, t11, t10, t7, t15, t8;
    
    t12 = -x2+x1;
    t11 = y1-y2;
    t10 = 2*w;
    t7 = w*w;
    t15 = (-t11*t7-t12*t10+t11)*a;
    t8 = (t7+1)*(x1*y2-x2*y1-t12*yc+t11*xc);
    poly.push_back( t15+t8 );
    poly.push_back( -2*b*(-t12*t7+t11*t10+t12) );
    poly.push_back( -t15+t8 );
    
    return CGAL::VORELL::primpart(typename ET::upolq_t(poly.begin(), poly.end()));
}*/


// deg = 6+6
template<class ET>
typename ET::bpolz_t Generated_code<ET>::bisector(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                                                  const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const
{
    std::vector<std::vector<QQ> > c(7);
    std::vector<typename ET::upolq_t> poly;
    
    for (std::size_t i = 0; i < c.size(); i++) c[i].reserve(7);

    QQ t1205, t741, t1404, t1424, t739, t1203, t1117, t1338, t1351, t1401, t737, t740, 
       t1210, t1104, t1372, t1206, t1114, t1400, t1016, t1239, t1211, t1126, t1102, t1402, 
       t1342, t1343, t1366, t1230, t1340, t1346, t1399, t1323, t1227, t1228, t1415, t1406, 
       t1253, t729, t1413, t1311, t1254, t1229, t1261, t1077, t1238, t1141, t1046, t1209, 
       t1119, t1345, t1243, t1157, t732, t1352, t1213, t1138, t1005, t883, t1395, t1068, 
       t1214, t1175, t1049, t1223, t1226, t1414, t1408, t1073, t919, t1252, t1195, t1052, 
       t925, t1356, t927, t1140, t1202, t1108, t1333, t1236, t1139, t1048, t1161, t1003, 
       t882, t1353, t1319, t929, t1216, t1308, t1344, t1258, t1310, t1394, t781, t735, 
       t721, t1381, t727, t829, t1241, t771, t753, t742, t1207, t1099, t979, t447, 
       t1136, t637, t738, t733, t1217, t1208, t1121, t647, t965, t1385, t1215, t718, 
       t1268, t1423, t1293, t755, t1422, t1421, t1420, t1256, t1200, t1109, t936, t1122, 
       t627, t987, t413, t1133, t646, t1246, t715, t1289, t1212, t1120, t1244, t1179, 
       t1127, t1035, t1397, t1204, t1125, t1231, t1105, t1373, t1367, t1324, t734, t1378, 
       t728, t1218, t722, t827, t1220, t1320, t1221, t1172, t1170, t1078, t1059, t935, 
       t1160, t1076, t1142, t1050, t920, t1158, t1341, t1186, t1054, t926, t1354, t1398, 
       t1355, t1061, t923, t1123, t995, t885, t1053, t1391, t1237, t1162, t1017, t1339, 
       t1376, t1164, t1259, t1247, t1137, t1045, t783, t774, t756, t1419, t757, t1418, 
       t736, t1271, t1417, t1416, t778, t1412, t777, t1411, t1410, t1409, t1255, t1201, 
       t634, t635, t458, t457, t1168, t1305, t1196, t1031, t1375, t1060, t1083, t1198, 
       t896, t895, t1392, t1361, t1051, t1313, t785, t775, t1371, t806, t1321, t766, 
       t1403, t805, t1322, t1151, t630, t421, t629, t433, t1284, t1275, t1118, t1396, 
       t1278, t1393, t1389, t1388, t1299, t780, t769, t1298, t1288, t1291, t1260, t1199, 
       t1197, t1369, t1042, t1171, t1043, t1368, t1124, t1032, t1334, t1085, t1087, t1317, 
       t1315, t1314, t1304, t1250, t1302, t1018, t1277, t1296, t1075, t1222, t720, t713, 
       t915, t1294, t1282, t1281, t1280, t1279, t1274, t1272, t897, t1265, t884, t1264, 
       t1225, t716, t1165, t1098, t1074, t914, t814, t811, t810, t809, t799, t798, 
       t797, t796, t795, t794, t793, t792, t770, t767, t765, t764, t762, t761, 
       t760, t751, t750, t749, t748, t747;

    t1205 = -xc2*xc2-xc1*xc1;
    t741  = w1*w1;
    t1404 = 4-4*t741;
    t1424 = t1404*w2;
    t739  = b1*b1;
    t1203 = t741*t739;
    t1117 = b2*t1203;
    t1338 = -2*yc1;
    t1351 = 2*yc2;
    t1401 = 2*t1338+2*t1351;
    t737  = w2*w2;
    t740  = a1*a1;
    t1210 = t737*t740;
    t1104 = t737*t1203;
    t1372 = t1117-b2*t1104;
    t1206 = t741*t740;
    t1114 = b2*t1206;
    t1400 = t739*b2-t1114;
    t1016 = w2*t1114;
    t1239 = w2*t740;
    t1211 = t737*t739;
    t1126 = b2*t1211;
    t1102 = t737*t1206;
    t1402 = t1126-b2*t1102;
    t1342 = 2*xc2;
    t1343 = -2*xc1;
    t1366 = t1342+t1343;
    t1230 = t739*w2;
    t1340 = 2*yc1;
    t1346 = -2*yc2;
    t1399 = 2*t1340+2*t1346;
    t1323 = t1399*(t1016+b2*t1239)+t1366*(-t1402+t1372+t1400+(t1210-t740)*b2)
              +t1401*(w2*t1117+b2*t1230);
    t1227 = yc1*xc2;
    t1228 = xc1*yc2;
    t1415 = t1227+t1228;
    t1406 = t1415*a1;
    t1253 = a1*b2;
    t729  = yc2*yc2;
    t1413 = t1205*b2;
    t1311 = -t729*t1253-t1413*a1;
    t1254 = w2*a1;
    t1229 = xc1*xc2;
    t1261 = w1*b2;
    t1077 = t1229*t1261;
    t1238 = t737*b2;
    t1141 = yc2*t1238;
    t1046 = a1*t1141;
    t1209 = t741*t737;
    t1119 = b2*t1209;
    t1345 = -4*xc2;
    t1243 = yc1*yc2;
    t1157 = a1*t1243;
    t732  = yc1*yc1;
    t1352 = -2*xc2;
    t1213 = t737*a1;
    t1138 = xc2*t1213;
    t1005 = t741*t1138;
    t883  = b2*t1005;
    t1395 = 4*yc2-4*yc1;
    t1068 = b2*t1138;
    t1214 = t741*a1;
    t1175 = b2*t1214;
    t1049 = yc1*t1175;
    t1223 = yc2*xc2;
    t1226 = xc1*yc1;
    t1414 = t1223+t1226;
    t1408 = t1414*b2;
    t1073 = xc2*t1175;
    t919  = w2*t1073;
    t1252 = yc1*b2;
    t1195 = a1*t1252;
    t1052 = t737*t1195;
    t925  = w1*t1052;
    t1356 = 2*b2;
    t927  = w1*t1046;
    t1140 = xc1*t1213;
    t1202 = w1*t1253;
    t1108 = w2*t1202;
    t1333 = -4*t1108;
    t1236 = t741*b2;
    t1139 = xc1*t1236;
    t1048 = a1*t1139;
    t1161 = a1*t1229;
    t1003 = yc1*t1119;
    t882  = a1*t1003;
    t1353 = 4*xc2;
    t1319 = w1+w2;
    t929  = w2*t1048;
    t1216 = t732*b2;
    t1308 = -a1*t1216+t1311;
    t1344 = 2*xc1;
    t1258 = w2*w1;
    t1310 = -t729-t1205;
    t1394 = t732-t1310;
    t781  = t925*t1345+t1068*t1343-t732*t1175+t883*t1344-t1308+t1161*t1356+
            t1311*t741-8*t1077*t1254+t882*t1346+t1046*t1340+t1049*t1351+(-2+8*t1258)*
            b2*t1157+t927*t1353+t1308*t737+t1048*t1352+t1319*(4*t1406*b2-4*t1408*a1)+(
            t1333+a1*t1119)*t1394+(t919-t929-t1140*t1261)*t1395;
    t735  = a1*t740;
    t721  = t735*t1238;
    t1381 = 4*t1258;
    t727  = t735*b2;
    t829  = t727*t1381-t735*t1236+t741*t721+t727-t721;
    t1241 = t739*a1;
    t771  = t781-t829+2*b2*t1241-2*t739*t1333-2*(t1126+t1372)*a1;
    t753  = -t1323-t771;
    t742  = a2*a2;
    t1207 = t737*t742;
    t1099 = t741*t1207;
    t979  = a1*t1099;
    t447  = b2*t979;
    t1136 = a1*t1207;
    t637  = b2*t1136;
    t738  = b2*b2;
    t733  = t738*b2;
    t1217 = t733*w1;
    t1208 = t741*t742;
    t1121 = a1*t1208;
    t647  = b2*t1121;
    t965  = t742*t1108;
    t1385 = 4*t733;
    t1215 = t742*a1;
    t718  = b2*t1215;
    t1268 = -5*t718-20*t965+5*t647+a1*t1385+16*t1217*t1254-5*t447+
            5*t637+(-4*t1214-t1404*t1213)*t733;
    t1423 = t753-t1268;
    t1293 = 4*t965+t447-t637-t647+t718;
    t755  = -t771+t1323;
    t1422 = t1293+t755;
    t1421 = -t1268+t755;
    t1420 = t753+t1293;
    t1256 = a2*w1;
    t1200 = b1*t1256;
    t1109 = w2*t1200;
    t936  = t740*t1109;
    t1122 = a2*t1210;
    t627  = b1*t1122;
    t987  = a2*t1102;
    t413  = b1*t987;
    t1133 = a2*t1206;
    t646  = b1*t1133;
    t1246 = b1*t740;
    t715  = a2*t1246;
    t1289 = 4*t936-t627-t646+t413+t715;
    t1212 = t737*t738;
    t1120 = b1*t1212;
    t1244 = t742*b1;
    t1179 = w1*t1244;
    t1127 = w1*t1207;
    t1035 = b1*t1127;
    t1397 = (-t738+t1207)*b1;
    t1204 = t741*t738;
    t1125 = b1*t1204;
    t1231 = t738*w1;
    t1105 = t737*t1204;
    t1373 = t1120-b1*t1105;
    t1367 = t1352+t1344;
    t1324 = t1367*(t1125+t1397+t1244-t1373+(-t1099-t1208)*b1)+t1399*(t1035+
            t1179)+t1401*(b1*t1231+w1*t1120);
    t734  = a2*t742;
    t1378 = -4*t1258;
    t728  = b1*t734;
    t1218 = t741*b1;
    t722  = t737*t728;
    t827  = t728*t1378+t734*t1218+t722-t741*t722-t728;
    t1220 = t738*a2;
    t1320 = -4*t1109;
    t1221 = t737*a2;
    t1172 = b1*t1221;
    t1170 = b1*t1229;
    t1078 = w2*t1170;
    t1059 = xc1*t1172;
    t935  = w1*t1059;
    t1160 = b1*t1243;
    t1076 = w2*t1160;
    t1142 = yc2*t1218;
    t1050 = a2*t1142;
    t920  = w2*t1050;
    t1158 = a2*t1223;
    t1341 = -4*xc1;
    t1186 = a2*t1218;
    t1054 = yc1*t1186;
    t926  = w2*t1054;
    t1354 = -2*a2;
    t1398 = t1205*a2;
    t1355 = 2*a2;
    t1061 = xc2*t1172;
    t923  = w1*t1061;
    t1123 = a2*t1209;
    t995  = xc1*t1123;
    t885  = b1*t995;
    t1053 = xc2*t1186;
    t1391 = -4*b1;
    t1237 = t729*b1;
    t1162 = b1*t1228;
    t1017 = b1*t1123;
    t1339 = 4*xc1;
    t1376 = 2-2*t737;
    t1164 = b1*t1227;
    t1259 = a2*b1;
    t1247 = t737*b1;
    t1137 = yc1*t1247;
    t1045 = a2*t1137;
    t783  = t1160*t1354+t885*t1342-t1310*t1017+t1053*t1343+t1045*t1351+t1170*
            t1355+t926*t1339+t1059*t1352+t920*t1341+(t1017+t1259-t1186-t1172)*t732+(-8*
            t1078+8*t1076)*t1256+(-1+t737+t741)*(-a2*t1237-t1398*b1)+t1376*yc1*t1050+
            t1394*t1320+t1319*(-4*t1226*t1259+t1158*t1391+4*(t1162+t1164)*a2)+(-t935+w2
            *t1053+t923)*t1395;
    t774  = -t783-t827+2*t738*t1320-2*b1*t1220+2*(t1125+t1373)*a2;
    t756  = t1324-t774;
    t1419 = -t1289+t756;
    t757  = -t1324-t774;
    t1418 = -t757+t1289;
    t736  = b1*t739;
    t1271 = -5*t646+5*t715-5*t627+5*t413+20*t936+(t1404*t1221+(
            -16*t1258-t1404)*a2)*t736;
    t1417 = t757-t1271;
    t1416 = -t1271+t756;
    t778  = t783-t827;
    t1412 = -t1289+t778;
    t777  = t781+t829;
    t1411 = t777-t1293;
    t1410 = -t1271+t778;
    t1409 = t777+t1268;
    t1255 = b1*b2;
    t1201 = w2*t1255;
    t634  = t742*t1201;
    t635  = b2*t1179;
    t458  = t741*t634;
    t457  = b2*t1035;
    t1168 = w1*t1247;
    t1305 = t1168*t1385+5*t635-5*t634+t1217*t1391+5*t458-5*t457+t733*
            b1*t1424;
    t1196 = yc2*t1255;
    t1031 = b1*t1119;
    t1375 = yc1-yc2;
    t1060 = b2*t1142;
    t1083 = b2*t1168;
    t1198 = b1*t1252;
    t896  = yc2*t1031;
    t895  = b1*t1003;
    t1392 = b2*t1237+(t1216+t1413)*b1;
    t1361 = t1076*t1356-t1392*w2+(-2*t1078+t1164)*b2;
    t1051 = b1*t1141;
    t1313 = t1051-b2*t1137;
    t785  = -t1361-t1394*t1083+(t1339+t1345)*t1196*t1258+(t896-t895-t1060-
            t1313)*xc2+(t1376*t1077-t1139*t1375+t1408)*b1-b2*t1162+t1361*t741+(t1031*t1375+
            t1313)*xc1+(t1392+((t1353+t1341)*w2-t1376*yc2)*t1198)*w1;
    t775  = -t785+t1305;
    t1371 = -t896+t1196;
    t806  = t923*t1356+(-t1051-t895+t1198-t1371)*a2+(t1366*t1200-2*t935+
            t1045-t1054+t1050)*b2;
    t1321 = -3*t806;
    t766  = t1321+t775;
    t1403 = -a2*t1104+t1122;
    t805  = (-2*t919+t882+t1052+2*t929-t1046-t1049-t1195)*b1+(t1367*t1201+
            t1060+t1371)*a1;
    t1322 = -3*t805;
    t1151 = b2*t1246;
    t630  = w1*t1151;
    t421  = b1*t1016;
    t629  = w2*t1151;
    t433  = t740*t1083;
    t1284 = -5*t630+5*t433+5*t629-5*t421+((4-4*t737)*w1-t1424)*
            t736*b2;
    t1275 = -t1322+t1284;
    t1118 = a1*t1212;
    t1396 = -t1215+t1118;
    t1278 = -t1322-t1284;
    t1393 = t979+t1121-a1*t1105;
    t1389 = -2*w1;
    t1388 = 2*w2;
    t1299 = -t458-t635+t634+t457;
    t780  = -t785-t1299;
    t769  = t780-t806;
    t1298 = -t630-t421+t629+t433;
    t1288 = -t805+t1298;
    t1291 = t805+t1298;
    t1260 = a2*a1;
    t1199 = xc1*t1260;
    t1197 = xc2*t1260;
    t1369 = -t1199+t1197;
    t1042 = a2*t1140;
    t1171 = a2*t1214;
    t1043 = xc1*t1171;
    t1368 = t1042-t1043;
    t1124 = a2*t1203;
    t1032 = w2*t1124;
    t1334 = 2*t1032;
    t1085 = a2*t1157;
    t1087 = a2*t1161;
    t1317 = 2*w1*t1085+t1087*t1389;
    t1315 = -2*w2*t1085+t1087*t1388;
    t1314 = t1199*t1381+t1197*t1378;
    t1304 = (-t1125+t1397)*a1+(t1393-t1396)*b1;
    t1250 = a2*t740;
    t1302 = (-t1250-t1403)*b2+(t1117+t1400+t1402)*a2;
    t1018 = w1*t1118;
    t1277 = t1375*(-t1136+(t738-t1204)*a1+t1393+t1396)+(w1*t1215+a1*t1127)*
            t1367+t1366*(t1018+a1*t1231);
    t1296 = -t806+t1299;
    t1075 = w2*t1171;
    t1222 = t734*a1;
    t720  = w1*t1222;
    t713  = w2*t1222;
    t915  = -t741*t713-t720+t737*t720+t713;
    t1294 = -2*t738*t1075-t915+t1018*t1355+(t1389+t1388)*a1*t1220;
    t1282 = 6*t1304;
    t1281 = 2*t1302;
    t1280 = 2*t1304;
    t1279 = 6*t1302;
    t1274 = t1321-t1305;
    t1272 = xc2*t1334+t1032*t1343+t1375*(-t1250-t1133+t1124+t987+t1403)+t1367
            *w2*t1133+(t1375*(t739-t1211)+t1367*t1239+t1366*t1230)*a2;
    t897  = a1*t995;
    t1265 = t925*t1355+t927*t1354+(-t897+t1368-t1369)*b2+((t1340+t1346)*t1202
            +t883-t1068+t1073)*a2;
    t884  = a2*t1005;
    t1264 = (t1338+t1351)*t1254*t1259+(-t884+t1368+t1369)*b1+(2*t920-2*
            t926+t1053+t885-t1061)*a1;
    t1225 = t735*a2;
    t716  = w2*t1225;
    t1165 = w1*t1221;
    t1098 = a1*t1158;
    t1074 = a1*t1165;
    t914  = w1*t1225-t716-t735*t1165+t741*t716;
    t814  = 6*t1264;
    t811  = 2*t1264;
    t810  = 6*t1265;
    t809  = 2*t1265;
    t799  = t814-t1282;
    t798  = t814+t1282;
    t797  = t810-t1279;
    t796  = t811-t1280;
    t795  = t810+t1279;
    t794  = t811+t1280;
    t793  = t809-t1281;
    t792  = t809+t1281;
    t770  = t780+t806;
    t767  = t785-t1291;
    t765  = t775-t1321;
    t764  = t785-t1275;
    t762  = -t1414*t1260+(t1098+(t1205*w2-t1414*t737-t1415)*t1260+t1315)*t741+
             t1406*a2+(t897-t1314)*yc2+(t1098+(t1205*w1-t1415)*t1260-t1317)*t737+(t1042+
             t1043+t884+t1314)*yc1-t1315+t1317+(-t1319*t1260+t1074+t1075)*(t729+t732)-t1398*
             t1319*a1;
    t761  = t762+t914;
    t760  = -2*t739*t1074+a1*t1334+t762-t914+2*(-w2+w1)*a2*t1241;
    t751  = t761+t1294;
    t750  = t760-t1272;
    t749  = t760+t1272;
    t748  = t749+t1294;
    t747  = t750+t1294;
    c[0].push_back( t767+t1296 );
    c[1].push_back( t794+t1419 );
    c[2].push_back( t767+t1274 );
    c[3].push_back( 2*t811+2*t1412 );
    c[4].push_back( t766+t1291 );
    c[5].push_back( t796-t1418 );
    c[6].push_back( t769+t1291 );
    c[0].push_back( t793-t1422 );
    c[1].push_back( 4*t748-4*t1277 );
    c[2].push_back( t797-t1421 );
    c[3].push_back( 8*t749+8*t915 );
    c[4].push_back( t797+t1421 );
    c[5].push_back( 4*t748+4*t1277 );
    c[6].push_back( t793+t1422 );
    c[0].push_back( t764+t1296 );
    c[1].push_back( t798+t1416 );
    c[2].push_back( t764+t1274 );
    c[3].push_back( 2*t814+2*t1410 );
    c[4].push_back( t766+t1275 );
    c[5].push_back( t799+t1417 );
    c[6].push_back( t769+t1275 );
    c[0].push_back( 2*t809+2*t1411 );
    c[1].push_back( 8*t751-8*t1277 );
    c[2].push_back( 2*t810+2*t1409 );
    c[3].push_back( 16*t761+16*t915 );
    c[4].push_back( 2*t810-2*t1409 );
    c[5].push_back( 8*t751+8*t1277 );
    c[6].push_back( 2*t809-2*t1411 );
    c[0].push_back( t770-t1278 );
    c[1].push_back( t798-t1416 );
    c[2].push_back( t765-t1278 );
    c[3].push_back( 2*t814-2*t1410 );
    c[4].push_back( -t766+t1278 );
    c[5].push_back( t799-t1417 );
    c[6].push_back( -t769+t1278 );
    c[0].push_back( t792-t1420 );
    c[1].push_back( 4*t747-4*t1277 );
    c[2].push_back( t795-t1423 );
    c[3].push_back( 8*t750+8*t915 );
    c[4].push_back( t795+t1423 );
    c[5].push_back( 4*t747+4*t1277 );
    c[6].push_back( t792+t1420 );
    c[0].push_back( t770+t1288 );
    c[1].push_back( t794-t1419 );
    c[2].push_back( t765+t1288 );
    c[3].push_back( 2*t811-2*t1412 );
    c[4].push_back( -t766-t1288 );
    c[5].push_back( t796+t1418 );
    c[6].push_back( -t769-t1288 );
    
    poly.reserve(7);
    for (std::size_t i = 0; i < c.size(); i++)
        poly.push_back( typename ET::upolq_t(c[i].begin(), c[i].end()) );

    return CGAL::VORELL::primpart(typename ET::bpolq_t(poly.begin(), poly.end()));
}

// deg = 2+2
template<class ET>
typename ET::bpolq_t Generated_code<ET>::polar_q(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                                               const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const
{
    std::vector<std::vector<QQ> > c(3);
    std::vector<typename ET::upolq_t> poly;
    
    for (std::size_t i = 0; i < c.size(); i++) c[i].reserve(3);

    QQ t1, t36, t67, t55, t69, t68, t37, t66, t54, t65, t53, t52, t49, t61, t64, 
       t60, t44, t42, t40, t39, t38;

    t36 = w1*w1;
    t67 = -t36-1;
    t55 = b2*w2;
    t69 = 2*yc1*t55;
    t68 = -2*t36;
    t37 = w2*w2;
    t66 = (-(2*t37-2)*w1-(2+t68)*w2)*a1;
    t54 = t36*b2;
    t65 = b2+t54;
    t53 = t37*b2;
    t52 = t37*t36;
    t49 = b2*t52;
    t61 = (4*w1*t55-t54-t53+t49+b2)*a1;
    t64 = t53+t49;
    t60 = t36*t69+t69+(-t64+t65)*xc1+(-2+t68)*yc2*t55+(t67*b2+t64)*xc2;
    t44 = (t64+t65)*a2;
    t42 = t67*((-2*xc2+2*xc1)*w2+t37*yc1+yc2)*a2;
    t40 = t67*(-yc2*t37-yc1);
    t39 = t44-t60;
    t38 = t44+t60;
    t1  = (t37*w1-w2*t36-w1+w2)*b1*b2;
    c[0].push_back( 2*t39-2*t61 );
    c[0].push_back( -8*t1 );
    c[0].push_back( 2*t39+2*t61 );
    c[1].push_back( 4*(-t66-t40)*a2-4*t42 );
    c[1].push_back( 8*(-1+t37-t52+t36-4*w2*w1)*a2*b1 );
    c[1].push_back( 4*(t66-t40)*a2-4*t42 );
    c[2].push_back( 2*t38+2*t61 );
    c[2].push_back( 8*t1 );
    c[2].push_back( 2*t38-2*t61 );

    poly.reserve(3);
    for (std::size_t i = 0; i < c.size(); i++)
        poly.push_back( typename ET::upolq_t(c[i].begin(), c[i].end()) );

    return typename ET::bpolq_t(poly.begin(), poly.end());
}

// deg = 2+2
template<class ET>
typename ET::bpolz_t Generated_code<ET>::polar(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                                               const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const
{
    return CGAL::VORELL::primpart(polar_q(a1, b1, w1, xc1, yc1, a2, b2, w2, xc2, yc2));
}

// deg = 2+2
template<class ET>
typename ET::bpolz_t Generated_code<ET>::tan_poly_cut(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                                               const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const
{
    std::vector<std::vector<QQ> > c(3);
    std::vector<typename ET::upolq_t> poly;
    
    for (std::size_t i = 0; i < c.size(); i++) c[i].reserve(3);

    QQ t36, t73, t58, t72, t71, t69, t37, t68, t54, t34, t61, t51, 
       t48, t67, t64, t53, t66, t63, t44, t42, t40, t39, t38, t1;

    t36 = w2*w2;
    t73 = 1+t36;
    t58 = b1*w1;
    t72 = 2*yc1*t58;
    t71 = -2*t36;
    t69 = t73*a1;
    t37 = w1*w1;
    t68 = (-(2+t71)*w1-(-2+2*t37)*w2)*a2;
    t54 = b1*t37;
    t34 = b1*a2;
    t61 = w1*w2;
    t51 = t37*t36;
    t48 = b1*t51;
    t67 = -t48+b1*t36;
    t64 = 4*t34*t61+t34+(-t54-t67)*a2;
    t53 = a1*t37;
    t66 = yc2*t53+(2*xc2-2*xc1)*a1*w1;
    t63 = t72+(-t54+b1+t67)*xc1+t36*t72+(t48+(t37-t73)*b1)*xc2+(-2+t71)*yc2*t58;
    t44 = a1*t48+(t53+t69)*b1;
    t42 = t69*yc1+t66*t36+t66;
    t40 = (yc1*t37+yc2)*t73;
    t39 = t44+t63;
    t38 = t44-t63;
    t1 = (w2-w1+w1*t36-t37*w2)*b1*b2;
    c[0].push_back( 2*t39-2*t64 );
    c[0].push_back( 4*(-t68-t40)*a1+4*t42 );
    c[0].push_back( 2*t38+2*t64 );
    c[1].push_back( 8*t1 );
    c[1].push_back( 8*(-t51-1-4*t61+t36+t37)*a1*b2 );
    c[1].push_back( -8*t1 );
    c[2].push_back( 2*t39+2*t64 );
    c[2].push_back( 4*(t68-t40)*a1+4*t42 );
    c[2].push_back( 2*t38-2*t64 );
    
    poly.reserve(3);
    for (std::size_t i = 0; i < c.size(); i++)
        poly.push_back( typename ET::upolq_t(c[i].begin(), c[i].end()) );

    return CGAL::VORELL::primpart(typename ET::bpolq_t(poly.begin(), poly.end()));
}    

// deg = 4+2
template<class ET> typename ET::bpolz_t 
Generated_code<ET>::normal_cut(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                               const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const
{
    std::vector<std::vector<QQ> > c(3);
    std::vector<typename ET::upolq_t> poly;
    
    for (std::size_t i = 0; i < c.size(); i++) c[i].reserve(5);

    QQ t33, t61, t44, t34, t50, t65, t64, t47, t54, t55, t59, t43, t63, t53, 
       t57, t52, t49, t48, t46, t42, t40, t39, t38, t37, t4, t3, t2, t1;

    t33 = w2*w2;
    t61 = t33+1;
    t44 = t61*yc2;
    t34 = w1*w1;
    t50 = t34*t33;
    t65 = -1-t50;
    t64 = t34+t33;
    t47 = a1*t50;
    t54 = a1*t34;
    t55 = w1*w2;
    t59 = (t47-t54+(4*t55+1-t33)*a1)*a2;
    t43 = t61*yc1;
    t63 = (xc1+(2*yc1-2*yc2)*w1)*a1;
    t53 = xc1-xc2;
    t57 = a2*w2;
    t52 = 2*t53;
    t49 = -2*t57;
    t48 = 2*t57;
    t46 = 2*t33;
    t42 = (-t34*t44-t43)*b1;
    t40 = (b1*b1-a1*a1)*(-t64+t65);
    t39 = -(xc2+xc1*t34)*a1*t61+xc2*t54+xc2*t47+t63*t33+t63;
    t38 = t39+t40;
    t37 = t39-t40;
    t4  = (-4*t55+t64+t65)*b2*b1;
    t3  = 8*(w1*t33+w2-w1-t34*w2)*b2*a1;
    t2  = (t48+(t49+t43)*t34+(-2*a2+(a2+t53)*t46+t52)*w1+t44)*b1+t42;
    t1  = (t49+(t48+t43)*t34+(2*a2+(-a2+t53)*t46+t52)*w1+t44)*b1+t42;
    c[0].push_back( t2 );
    c[0].push_back( 2*t38-2*t59 );
    c[0].push_back( 0 );
    c[0].push_back( 2*t37-2*t59 );
    c[0].push_back( -t2 );
    c[1].push_back( -2*t4 );
    c[1].push_back( t3 );
    c[1].push_back( 0 );
    c[1].push_back( t3 );
    c[1].push_back( 2*t4 );
    c[2].push_back( t1 );
    c[2].push_back( 2*t38+2*t59 );
    c[2].push_back( 0 );
    c[2].push_back( 2*t37+2*t59 );
    c[2].push_back( -t1 );

    poly.reserve(3);
    for (std::size_t i = 0; i < c.size(); i++)
        poly.push_back( typename ET::upolq_t(c[i].begin(), c[i].end()) );

    return CGAL::VORELL::primpart(typename ET::bpolq_t(poly.begin(), poly.end()));
}

// deg = 12
template<class ET>
typename ET::upolz_t Generated_code<ET>::parnorm(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                                                 const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const
{
    if (a1 == b1 && a2 == b2) {
        std::vector<QQ> p;
        p.reserve(3);
        
        p.push_back( yc2 - yc1 );
        p.push_back( 2*(xc1-xc2) );
        p.push_back( -p[0] );
        
        return CGAL::VORELL::primpart(typename ET::upolq_t(p.begin(), p.end()));
    } 

    std::vector<QQ> poly;
    poly.reserve(13);
    
    // compute determinant of Sylvester matrix

    std::vector<QQ> poly1, poly2, poly3, poly4, poly5;
    poly1.reserve(3);
    poly2.reserve(3);
    poly3.reserve(5);
    poly4.reserve(5);
    poly5.reserve(5);

    QQ t36, t67, t73, t49, t35, t59, t72, t71, t70, t69, t60, t52, t66, t68,
             t57, t64, t58, t55, t54, t53, t47, t46, t45, t44, t43, t42, t41, t40,
             t39, t6, t5, t4, t3, t2, t1;
    
    t36 = w2*w2;
    t67 = -t36-1;
    t73 = 4*w1*w2;
    t49 = t67*yc1;
    t35 = w1*w1;
    t59 = t35*t36;
    t72 = 1+t59;
    t71 = t67*yc2;
    t70 = -t35-t36;
    t69 = (1+t35)*a1*a1;
    t60 = a1*t35;
    t52 = a1*t59;
    t66 = (-t60+t52+(1-t36+t73)*a1)*a2;
    t68 = (xc1+(-2*yc2+2*yc1)*w1)*a1;
    t57 = xc2-xc1;
    t64 = a2*w2;
    t58 = 2*t57;
    t55 = -2*t64;
    t54 = 2*t64;
    t53 = 2*t36;
    t47 = (t35*t49+t71)*b1;
    t46 = -t35*w2-w1+w1*t36+w2;
    t45 = t73+t70+t72;
    t44 = t69*t36+(t70-t72)*b1*b1+t69;
    t43 = a1*t46;
    t42 = t45*b1;
    t41 = xc2*t52+(xc2+t35*xc1)*t67*a1+xc2*t60+t68*t36+t68;
    t40 = t41+t44;
    t39 = t41-t44;
    t6  = t46*b2*b1;
    t5  = b2*t42;
    t4  = a2*t42;
    t3  = 16*b2*t43;
    t2  = (t55+(t54-t71)*t35+(2*a2+(-a2+t57)*t53+t58)*w1-t49)*b1+t47;
    t1  = (t54+(t55-t71)*t35+(-2*a2+(a2+t57)*t53+t58)*w1-t49)*b1+t47;
    poly1.push_back( -8*t6 );
    poly1.push_back( 8*t45*b2*a1 );
    poly1.push_back( 8*t6 );
    
    poly2.push_back( 8*t4 );
    poly2.push_back( 32*a2*t43 );
    poly2.push_back( -8*t4 );
    
    poly3.push_back( -2*t1 );
    poly3.push_back( 4*t40+4*t66 );
    poly3.push_back( 0 );
    poly3.push_back( 4*t39+4*t66 );
    poly3.push_back( 2*t1 );
    
    poly4.push_back( 4*t5 );
    poly4.push_back( t3 );
    poly4.push_back( 0 );
    poly4.push_back( t3 );
    poly4.push_back( -4*t5 );
    
    poly5.push_back( -2*t2 );
    poly5.push_back( 4*t40-4*t66 );
    poly5.push_back( 0 );
    poly5.push_back( 4*t39-4*t66 );
    poly5.push_back( 2*t2 );
    
    upolq_t p1 = typename ET::upolq_t(poly1.begin(), poly1.end());
    upolq_t p2 = typename ET::upolq_t(poly2.begin(), poly2.end());
    upolq_t p3 = typename ET::upolq_t(poly3.begin(), poly3.end());
    upolq_t p4 = typename ET::upolq_t(poly4.begin(), poly4.end());
    upolq_t p5 = typename ET::upolq_t(poly5.begin(), poly5.end());
                                                                              // = 2*p5 <-- faster? but how to express?
    upolz_t res = CGAL::VORELL::primpart((p3*p5*p2+(p3-p5)*p4*p1)*p2+(p5*p5-p4*p4+(p5+p5+p3)*p3)*p1*p1);
    
    if (a1 == b1) {
        std::vector<typename ET::ZZ> p;
        p.reserve(5);
        // (x^2+1)^2 = x^4+2*x^2+1
        p.push_back( 1 );
        p.push_back( 0 );
        p.push_back( 2 );
        p.push_back( 0 );
        p.push_back( 1 );
        upolz_t t2p1("x^4+2*x^2+1");
        return typename ET::PTZ::Pseudo_division_quotient()(res, typename ET::upolz_t(p.begin(), p.end()));
    }
    
    CGAL_assertion( a2 != b2);
    // TODO: a2 == b2, though not used
    
    return res;

/*    return (p1*p1*p5*p5) - (p1*p4*p2*p5) - (p1*p1*p4*p4) + (2*p3*p1*p1*p5) + 
           (p3*p2*p2*p5) + (p3*p2*p1*p4) + (p3*p3*p1*p1);*/

}

// intersection of two ellipses
// deg = 4
template<class ET>
typename ET::upolz_t Generated_code<ET>::inter_poly(
        const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
        const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const
{
    
    std::vector<QQ> poly;
    poly.reserve(5);

    QQ t218, t423, t219, t213, t349, t215, t376, t238, t214, t220, t355,
       t360, t212, t396, t354, t221, t421, t426, t364, t425, t420, t407, t344,
       t341, t419, t418, t403, t353, t415, t414, t327, t412, t343, t348, t347,
       t345, t310, t309, t409, t408, t322, t346, t321, t406, t342, t300,
       t337, t298, t399, t339, t208, t211, t304, t266, t365, t166, t331, t313,
       t269, t398, t338, t165, t303, t311, t267, t308, t305, t268, t352, t323,
       t397, t216, t394, t391, t390, t297, t296, t389, t388, t387, t318, t386,
       t207, t201, t340, t372, t306, t362, t263, t259, t307, t333, t371,
       t359, t357, t334, t314, t312, t289, t270, t243, t240, t227, t226, t217;

    t218 = w2*w2;
    t423 = xc1-xc2;
    t219 = w1*w1;
    t213 = w1*t219;
    t349 = -t213-w1;
    t215 = w2*t218;
    t376 = t215-w2;
    t238 = 8*t376*t349;
    t214 = b2*b2;
    t220 = t218*t218;
    t355 = -4*t220-4;
    t360 = 8*t218;
    t212 = t219*t219;
    t396 = t376*(4*t212-4);
    t354 = -2-2*t220;
    t221 = a2*a2;
    t421 = yc1-yc2;
    t426 = ((t421*(t396-t349*(t360+t355))+t423*(-(8-8*t212)*t218+t238))*
             t221+(-t421*(-16*t349*t218+t396)+t423*(-t238+(t212-1)*
            (-4*t218-t354)))*t214)*a1;
    t364 = 4*yc1;
    t425 = t364-4*yc2;
    t420 = 4*yc2-4*yc1;
    t407 = xc1*xc1+xc2*xc2;
    t344 = t218*t214;
    t341 = t218*t221;
    t419 = -2*t344+4*t341;
    t418 = t425*w2*t214;
    t403 = t213-w1;
    t353 = 2+2*t212;
    t415 = (t221-t214)*(-t403*(12*t218+t354)-t376*(-12*t219+t353))*a1;
    t414 = t407*t214;
    t327 = 4*t344;
    t412 = -2*t341+t327;
    t343 = t221*t215;
    t348 = yc1*t221;
    t347 = yc2*t221;
    t345 = t215*t214;
    t310 = yc1*t345;
    t309 = yc2*t345;
    t409 = xc1*t327-t418+4*t310-4*t309+t420*t343+(4*t348-4*t347)*w2;
    t408 = t376*t403;
    t322 = xc2*b1*t214;
    t346 = xc1*t214;
    t321 = b1*t346;
    t406 = -2*t321+2*t322;
    t342 = t219*t214;
    t300 = t218*t342;
    t337 = t219*t221;
    t298 = t218*t337;
    t399 = 8*t300-4*t298+t221;
    t339 = t221*t212;
    t208 = yc2*yc2;
    t211 = yc1*yc1;
    t304 = xc1*t342;
    t266 = xc2*t304;
    t365 = -2*xc2;
    t166 = t346*t365;
    t331 = -2*t214*t221;
    t313 = yc1*t347;
    t269 = t219*t313;
    t398 = -4*t269+t219*t331+t208*t339-4*t266+t166+t212*t166+2*t407*
           t342+2*(t208+t211)*t337;
    t338 = t221*t220;
    t165 = -2*t313;
    t303 = yc2*t341;
    t311 = yc1*t344;
    t267 = yc2*t311;
    t308 = xc1*t343;
    t305 = xc1*t341;
    t268 = xc2*t305;
    t352 = t221*w2;
    t323 = xc1*t352;
    t397 = t218*t331+t303*t364-8*t267+t165-8*t268+(t338+t412)*t211+t412*
           t208+t425*t308+(t165+t414)*t220+t418*xc1+t409*xc2+t420*
           (xc1*t345+t323)+t419*t407;
    t216 = a1*a1;
    t394 = t419*t216;
    t391 = t365+2*xc1;
    t390 = 2*xc2-2*xc1;
    t297 = t216*t337;
    t296 = t216*t342;
    t389 = 4*t297-2*t296;
    t388 = 8*t408;
    t387 = -16*t408;
    t318 = b1*t339;
    t386 = yc2*t318+b1*t348-t349*t406;
    t207 = t214*t216;
    t201 = t220*t207;
    t340 = t214*t212;
    t372 = -t216*t340-t212*t201-t207-t201;
    t306 = xc2*t343;
    t362 = 8*xc2;
    t263 = t219*t309;
    t259 = t219*t310;
    t307 = xc1*t337;
    t333 = xc2*t352;
    t371 = (t269+t266)*t360+t398-16*(t267+t268)*t219-8*xc2*t263+(8*t263
            -8*t259)*xc1+t397+t397*t212+(t338+t399)*t208-4*(t221+t407)*t300+
            (t339+t399)*t211+t398*t220+t259*t362+t407*(t340+8*t298)-t421*
            (-8*t215*t307+8*(t307+xc2*t342)*w2-8*w2*t304+
            (-8*t333+8*t306)*t219)+t414;
    t359 = -4-4*t212;
    t357 = -1-t220;
    t334 = -4*t219+t353;
    t314 = b1*t341;
    t312 = b1*t340;
    t289 = -4*b1*t344;
    t270 = t212*t314;
    t243 = (-t212-t220*t212+t357)*t214;
    t240 = (yc2+yc1*t212)*t357*t221;
    t227 = (t289-2*t270)*yc2+t386*t220+t420*t218*t312+(t390*t312+t391*t318)
            *w2+(t390*t318+t391*t312+t406)*t215+t386-t349*(t314*t362+xc2*t289)+
            2*(t270-t314)*yc1+2*(-t322+t321)*w2+(4*t311-2*t323-2*t306-
            t349*(-8*t305+t409)+2*t308+2*t333+2*t303)*b1;
    t226 = t221*t243+t371-t372+t389*t220+t394*t212+t389+(-16*t297+20*t296
           )*t218+t394+8*((-t343+t352)*t216+t376*t207)*t403;
    t217 = b1*b1;
    poly.push_back( t226-t426 );
    poly.push_back( 4*(t240-t415)*b1+4*t227 );
    poly.push_back( 2*(t243+(t355*t219+(16*t219+t359)*t218+t388)*t216+(t334*
              t220+(40*t219+t359)*t218+t334-t387)*t217)*t221+2*
              ((-t354*t219+(-20*t219+t353)*t218-t388)*t216+((8+8*t220)*t219+
              (8+8*t212-32*t219)*t218+t387)*t217)*t214+2*t371+2*t372 );
    poly.push_back( 4*t227+4*(t240+t415)*b1 );
    poly.push_back( t226+t426 );
    
    return CGAL::VORELL::primpart(
            typename ET::upolq_t(poly.begin(), poly.end()));
}

// deg = 12
template<class ET>
typename ET::upolq_t Generated_code<ET>::parnorm_testing(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                                                 const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const
{
    std::vector<QQ> poly;
    poly.reserve(13);
    
    // compute determinant of Sylvester matrix

    std::vector<QQ> poly1, poly2, poly3, poly4, poly5;
    poly1.reserve(3);
    poly2.reserve(3);
    poly3.reserve(5);
    poly4.reserve(5);
    poly5.reserve(5);

    QQ t36, t67, t73, t49, t35, t59, t72, t71, t70, t69, t60, t52, t66, t68,
             t57, t64, t58, t55, t54, t53, t47, t46, t45, t44, t43, t42, t41, t40,
             t39, t6, t5, t4, t3, t2, t1;
    
    t36 = w2*w2;
    t67 = -t36-1;
    t73 = 4*w1*w2;
    t49 = t67*yc1;
    t35 = w1*w1;
    t59 = t35*t36;
    t72 = 1+t59;
    t71 = t67*yc2;
    t70 = -t35-t36;
    t69 = (1+t35)*a1*a1;
    t60 = a1*t35;
    t52 = a1*t59;
    t66 = (-t60+t52+(1-t36+t73)*a1)*a2;
    t68 = (xc1+(-2*yc2+2*yc1)*w1)*a1;
    t57 = xc2-xc1;
    t64 = a2*w2;
    t58 = 2*t57;
    t55 = -2*t64;
    t54 = 2*t64;
    t53 = 2*t36;
    t47 = (t35*t49+t71)*b1;
    t46 = -t35*w2-w1+w1*t36+w2;
    t45 = t73+t70+t72;
    t44 = t69*t36+(t70-t72)*b1*b1+t69;
    t43 = a1*t46;
    t42 = t45*b1;
    t41 = xc2*t52+(xc2+t35*xc1)*t67*a1+xc2*t60+t68*t36+t68;
    t40 = t41+t44;
    t39 = t41-t44;
    t6  = t46*b2*b1;
    t5  = b2*t42;
    t4  = a2*t42;
    t3  = 16*b2*t43;
    t2  = (t55+(t54-t71)*t35+(2*a2+(-a2+t57)*t53+t58)*w1-t49)*b1+t47;
    t1  = (t54+(t55-t71)*t35+(-2*a2+(a2+t57)*t53+t58)*w1-t49)*b1+t47;
    poly1.push_back( -8*t6 );
    poly1.push_back( 8*t45*b2*a1 );
    poly1.push_back( 8*t6 );
    
    poly2.push_back( 8*t4 );
    poly2.push_back( 32*a2*t43 );
    poly2.push_back( -8*t4 );
    
    poly3.push_back( -2*t1 );
    poly3.push_back( 4*t40+4*t66 );
    poly3.push_back( 0 );
    poly3.push_back( 4*t39+4*t66 );
    poly3.push_back( 2*t1 );
    
    poly4.push_back( 4*t5 );
    poly4.push_back( t3 );
    poly4.push_back( 0 );
    poly4.push_back( t3 );
    poly4.push_back( -4*t5 );
    
    poly5.push_back( -2*t2 );
    poly5.push_back( 4*t40-4*t66 );
    poly5.push_back( 0 );
    poly5.push_back( 4*t39-4*t66 );
    poly5.push_back( 2*t2 );
    
    upolq_t p1 = typename ET::upolq_t(poly1.begin(), poly1.end());
    upolq_t p2 = typename ET::upolq_t(poly2.begin(), poly2.end());
    upolq_t p3 = typename ET::upolq_t(poly3.begin(), poly3.end());
    upolq_t p4 = typename ET::upolq_t(poly4.begin(), poly4.end());
    upolq_t p5 = typename ET::upolq_t(poly5.begin(), poly5.end());
                                                                              // = 2*p5 <-- faster? but how to express?
    return (p3*p5*p2+(p3-p5)*p4*p1)*p2+(p5*p5-p4*p4+(p5+p5+p3)*p3)*p1*p1;
}

template<class T> inline T det3(const T& a1, const T& b1, const T& c1, 
                         const T& a2, const T& b2, const T& c2, 
                         const T& a3, const T& b3, const T& c3)
{ 
    return (a1*b2-a2*b1)*c3 + (a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2;
}

template<class ET>
void Generated_code<ET>::ell_normal(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc, 
        upolq_t& poly_x_, upolq_t& poly_y_, upolq_t& poly_c_)
{
    QQ t17, t16, t15, t5, t14, t9, t13, t12, t11, t4, t2;
    
    t17 = a*w;
    t16 = (w+1)*(w-1);
    t15 = b*w;
    t5 = a*xc;
    t14 = t5+2*yc*t17;
    t9 = w*w;
    t13 = b*(yc*t9+2*w*xc-yc);
    t12 = a*a-b*b;
    t11 = b*t16;
    t4 = -8*t17;
    t2 = 4*a*t16;

    std::vector<QQ> poly_x;
    poly_x.reserve(5);
    poly_x.push_back( -4*t15 );
    poly_x.push_back( t2 );
    poly_x.push_back( 0 );
    poly_x.push_back( t2 );
    poly_x.push_back( 4*t15 );
    
    std::vector<QQ> poly_y;
    poly_y.reserve(5);
    poly_y.push_back( -2*t11 );
    poly_y.push_back( t4 );
    poly_y.push_back( 0 );
    poly_y.push_back( t4 );
    poly_y.push_back( 2*t11 );
    
    std::vector<QQ> poly_c;
    poly_c.reserve(5);
    poly_c.push_back( 2*t13 );
    poly_c.push_back( 4*(-t5+t12)*t9+4*t12+4*t14 );
    poly_c.push_back( 0 );
    poly_c.push_back( 4*(-t5-t12)*t9-4*t12+4*t14 );
    poly_c.push_back( -2*t13 );
    
    poly_x_ = typename ET::upolq_t(poly_x.begin(), poly_x.end());
    poly_y_ = typename ET::upolq_t(poly_y.begin(), poly_y.end());
    poly_c_ = typename ET::upolq_t(poly_c.begin(), poly_c.end());
}

template<class ET>
typename ET::tpolz_t Generated_code<ET>::ell_normals(const QQ& a1_, const QQ& b1_, const QQ& w1, const QQ& xc1, const QQ& xy1, 
		    const QQ& a2_, const QQ& b2_, const QQ& w2, const QQ& xc2, const QQ& xy2,
		    const QQ& a3_, const QQ& b3_, const QQ& w3, const QQ& xc3, const QQ& xy3)
{
    typename ET::TPTQ::Move move;
    upolq_t a1, b1, c1, a2, b2, c2, a3, b3, c3;
    
    ell_normal(a1_, b1_, w1, xc1, xy1, a1, b1, c1);
    ell_normal(a2_, b2_, w2, xc2, xy2, a2, b2, c2);
    ell_normal(a3_, b3_, w3, xc3, xy3, a3, b3, c3);
    
    tpolq_t a1t, b1t, c1t;
    tpolq_t a2r, b2r, c2r;
    tpolq_t a3s, b3s, c3s;
    
/*    a1t = move(typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(a1), 0, 2);
    std::cerr << "a1t = " << a1t << std::endl;
    b1t = move(typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(b1), 0, 2);
    c1t = move(typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(c1), 0, 2);*/
    a1t = typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(a1);
//    std::cerr << "a1t = " << a1t << std::endl;
    b1t = typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(b1);
    c1t = typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(c1);
    
    
    a2r = move(typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(a2), 0, 1);
//    std::cerr << "a2r = " << a2r << std::endl;
    b2r = move(typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(b2), 0, 1);
    c2r = move(typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(c2), 0, 1);
    
    /*a3s = typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(a3);
    std::cerr << "a3s = " << a3s << std::endl;
    b3s = typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(b3);
    c3s = typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(c3);*/
    a3s = move(typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(a3), 0, 2);
//    std::cerr << "a3s = " << a3s << std::endl;
    b3s = move(typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(b3), 0, 2);
    c3s = move(typename CGAL::Coercion_traits<tpolq_t, upolq_t>::Cast()(c3), 0, 2);
    
    // determinant ;-)
    tpolq_t Q = CGAL::VORELL::det3( a1t, b1t, c1t, a2r, b2r, c2r, a3s, b3s, c3s);
            // a1t*b2r*c3s - a1t*c2r*b3s + b1t*c2r*a3s - b1t*a2r*c3s + c1t*a2r*b3s - c1t*b2r*a3s;
    return CGAL::VORELL::primpart(Q);
}

// deg = 2+2 D(t,r) -> D(x,y) = y(x)
template<class ET>
typename ET::bpolz_t Generated_code<ET>::denom(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1, 
                                               const QQ& a2, const QQ& b2, const QQ& w2, const QQ& xc2, const QQ& yc2) const
{
    std::vector<std::vector<QQ> > c(3);
    std::vector<typename ET::upolq_t> poly;
    
    for (std::size_t i = 0; i < c.size(); i++) c[i].reserve(3);

    QQ t6, t7, t12, t8, t3, t2, t1;

    t6 = w1*w1;
    t7 = w2*w2;
    t12 = -w2+w1-w1*t7+t6*w2;
    t8 = t7-t6*t7+t6-1-4*w1*w2;
    t3 = t12*b1*b2;
    t2 = t8*b1*a2;
    t1 = t8*a1*b2;
    c[0].push_back( -t3 );
    c[0].push_back( t1 );
    c[0].push_back( t3 );
    c[1].push_back( -t2 );
    c[1].push_back( -4*t12*a1*a2 );
    c[1].push_back( t2 );
    c[2].push_back( t3 );
    c[2].push_back( -t1 );
    c[2].push_back( -t3 );

    poly.reserve(3);
    for (std::size_t i = 0; i < c.size(); i++)
        poly.push_back( typename ET::upolq_t(c[i].begin(), c[i].end()) );

    return CGAL::VORELL::primpart(typename ET::bpolq_t(poly.begin(), poly.end()));
}    

template<class ET>
typename ET::upolq_t Generated_code<ET>::numerxt(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1) const
{
    std::vector<QQ> poly;
    QQ t2, t4, t3;
    
    poly.reserve(3);
    t2 = w1*w1;
    t4 = (1+t2)*xc1;
    t3 = (1-t2)*a1;
    poly.push_back( t3+t4 );
    poly.push_back( -4*b1*w1 );
    poly.push_back( -t3+t4 );

    return typename ET::upolq_t(poly.begin(), poly.end());
}

template<class ET>
typename ET::upolq_t Generated_code<ET>::numeryt(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1) const
{
    std::vector<QQ> poly;
    QQ t2, t4, t3;
    
    poly.reserve(3);
    t4 = a1*w1;
    t2 = w1*w1;
    t3 = (1+t2)*yc1;
    poly.push_back( 2*t4+t3 );
    poly.push_back( 2*(1-t2)*b1 );
    poly.push_back( -2*t4+t3 );

    return typename ET::upolq_t(poly.begin(), poly.end());
}

template<class ET>
typename ET::upolq_t Generated_code<ET>::denomt(const QQ& a1, const QQ& b1, const QQ& w1, const QQ& xc1, const QQ& yc1) const
{
    std::vector<QQ> poly;
    QQ t1;
    
    poly.reserve(3);
    t1 = 1+w1*w1;
    poly.push_back( t1 );
    poly.push_back( 0 );
    poly.push_back( t1 );

    return typename ET::upolq_t(poly.begin(), poly.end());
}


// (-x2r*x2r-y2r*y2r+2*(x*x2r+y2r*y)*dr)*dt*dt+(x1t*x1t+y1t*y1t+2*(-x*x1t-y1t*y)*dt)*dr*dr;

// deg = 2+1 M(t,x) -> M(x,y) = y(x)
template<class ET>
typename ET::bpolz_t Generated_code<ET>::medial_x(const QQ& a, const QQ& b, const QQ& w, const QQ& xc, const QQ& yc) const
{
    std::vector<std::vector<QQ> > c(2);
    std::vector<typename ET::upolq_t> poly;

    for (std::size_t i = 0; i < c.size(); i++) c[i].reserve(3);

    QQ t1, t2, t5, t6, t7, t8;

    t8 = xc*a;
    t7 = a*a-b*b;
    t6 = -t8+t7;
    t5 = -t8-t7;
    t2 = w*w;
    t1 = (t2+1)*a;
    c[0].push_back( t6*t2+t5 );
    c[0].push_back( 0 );
    c[0].push_back( t5*t2+t6 );
    c[1].push_back( t1 );
    c[1].push_back( 0 );
    c[1].push_back( t1 );

    poly.reserve(2);
    for (std::size_t i = 0; i < c.size(); i++)
        poly.push_back( typename ET::upolq_t(c[i].begin(), c[i].end()) );

    return CGAL::VORELL::primpart(typename ET::bpolq_t(poly.begin(), poly.end()));
}


} // namespace VORELL

} //namespace CGAL
#endif
