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

#ifndef CGAL_ELLIPSE_2_H
#define CGAL_ELLIPSE_2_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/determinant.h>

#include <CGAL/OpenMP.h>

namespace CGAL {
    template<class ET> class Ellipse_2; // fwd delaration
}

#include<CGAL/Voronoi_diagram_of_ellipses_2/Generated_code.h>

#define ELL_PARAM_COEFFS(x) x.parametric_coefficient(0), \
                            x.parametric_coefficient(1), \
                            x.parametric_coefficient(2), \
                            x.parametric_coefficient(3), \
                            x.parametric_coefficient(4)
#define ELL_PARAM_COEFF3(x) x.parametric_coefficient(0), \
                            x.parametric_coefficient(1), \
                            x.parametric_coefficient(2)

namespace CGAL {

const int ELLIPSE_DEFAULT_RES = 64;

template<class ET>
class Ellipse_2 {

    typedef typename ET::QQ QQ;
    typedef typename ET::Kernel_scr::Point_2 Point_2;
    typedef typename ET::Kernel_scr::Segment_2 Segment_2;
    typedef typename ET::Root Root;
    typedef typename ET::upolz_t upolz_t;
    typedef typename CGAL::VORELL::Generated_code<ET> gencode;

    QQ pc_[5];   // parametric coeffs: a, b, w, xc, yc
    QQ ac_[6];   // algebraic coeffs: a, b, c, d, e, f (remember 2b, 2d, 2e !)
    QQ J1_, J2_;  // invariants J1 = a+c, J2 = ac - b^2
    
    bool is_circle_;
    
    int res; // draw resolution
    int eid; // ellipse-id
    static int eid_count;
    static OMP_lock mutex;
    
    template<class ET2> 
    friend std::ostream& operator<<(std::ostream& o, 
                                        const Ellipse_2<ET2>& x) {
        return o << x.pc_[0] << ' ' << x.pc_[1] << ' ' << x.pc_[2] << ' ' 
                  << x.pc_[3] << ' ' << x.pc_[4];
    }

    template<class ET2> 
    friend std::istream& operator>>(std::istream& i, 
                                        Ellipse_2<ET2>& x) {
        QQ pc[5];
        i >> pc[0] >> pc[1] >> pc[2] >> pc[3] >> pc[4];
        x = Ellipse_2<ET2>(pc[0], pc[1], pc[2], pc[3], pc[4]);
        return i;
    }


    // maybe TODO: boundary_xy (that is, both)
    QQ boundary_x(const QQ& a, const QQ& b, const QQ& w, 
                    const QQ& xc, const QQ& yc, const QQ& t) const {
        QQ t1, t2, t3, t4;
        // maple gen code
        t4 = -a+xc;
        t3 = a+xc;
        t2 = t*t;
        t1 = w*w;
        return (-4*b*w*t+t4*t2+(t3*t2+t4)*t1+t3)/((1+t1)*(1+t2));
    }

    QQ boundary_y(const QQ& a, const QQ& b, const QQ& w, 
                    const QQ& xc, const QQ& yc, const QQ& t) const {
        QQ t1, t2, t3, t4, t5;
        // maple gen code
        t2 = t*t;
        t5 = 1+t2;
        t4 = b*t;
        t3 = t5*yc;
        t1 = w*w;
        return -(-2*t4+(2*t2-2)*w*a+(2*t4-t3)*t1-t3)/((1+t1)*t5);
    }

#if 0
     // maybe TODO: evolute_xy (that is, both)
    QQ evolute_x(const QQ& a, const QQ& b, const QQ& w, 
                    const QQ& xc, const QQ& yc, const QQ& t) const {
        QQ t20, t3, t19, t7, t6, t18, t11, t17;
        // maple gen code
        t20 = -b*b+a*a;
        t3 = w*w;
        t19 = b*(1+t3);
        t7 = t*t;
        t6 = t7*t7;
        t18 = -1-3*t6;
        t11 = t*t7;
        t17 = -3*t7-t11*t11;
        return ((-t3+1)*t20*(t17-t18)*b+(16*t20*t11*w+xc*(-t17-t18)*t19)*a)/
               (2*t7+t6+1)/a/(1+t7)/t19;
    }
    
    QQ evolute_y(const QQ& a, const QQ& b, const QQ& w, 
                    const QQ& xc, const QQ& yc, const QQ& t) const {
        QQ t21, t5, t19, t9, t12, t8, t6;
        // maple gen code
        t21 = -a*a+b*b;
        t5 = w*w;
        t19 = 1/(t5+1)/a;
        t9 = t*t;
        t12 = t*t9;
        t8 = t9*t9;
        t6 = t12*t12;
        return -(t21*(8*t5-8)*a*t12+(yc*(-3*t9-3*t8-t6-1)/t19+t21*
                (-2*t6-6*t9+6*t8+2)*w)*b)/(2*t9+t8+1)/b/(1+t9)*t19;
    }
#endif
            
    // maybe TODO: major_medial_xy (that is, both)
    QQ major_medial_x(const QQ& a, const QQ& b, const QQ& w, 
                    const QQ& xc, const QQ& yc, const QQ& t) const {
        QQ t8, t7, t6, t5, t4, t3;
        // maple gen code
        t8 = a*a-b*b;
        t7 = a*xc;
        t6 = t7-t8;
        t5 = t7+t8;
        t4 = t*t;
        t3 = w*w;
        return (t6*t3+(t5*t3+t6)*t4+t5)/a/(1+t4+t3+t4*t3);
    }
    
    QQ major_medial_y(const QQ& a, const QQ& b, const QQ& w, 
                    const QQ& xc, const QQ& yc, const QQ& t) const {
        QQ t11, t12, t15;
        // maple gen code
        t11 = w*w;
        t12 = t*t;
        t15 = 1/(-t11-t12*t11-t12-1)/a;
        return ((a*a-b*b)*(2*t12-2)*w+yc/t15)*t15;
    }
    
    // maybe TODO: minor_medial_xy (that is, both)
    QQ minor_medial_x(const QQ& a, const QQ& b, const QQ& w, 
                    const QQ& xc, const QQ& yc, const QQ& t) const {
        QQ t1, t2, t3, t4;
        // maple gen code
        t4 = w*t;
        t1 = w*w;
        t2 = t*t;
        t3 = 1+t2+t1+t1*t2;
        return -(-4*a*a*t4+(4*b*t4-t3*xc)*b)/b/t3;
    }

    QQ minor_medial_y(const QQ& a, const QQ& b, const QQ& w, 
                    const QQ& xc, const QQ& yc, const QQ& t) const {
        QQ t7, t8, t11;
        // maple gen code
        t7 = w*w;
        t8 = t*t;
        t11 = b*(-1-t8-t7-t7*t8);
        return (yc*t11+(-b*b+a*a)*t*(2-2*t7))/t11;
    }
    
    void algebraic_coefficients(const QQ& a, const QQ& b, const QQ& w, 
                        const QQ& xc, const QQ& yc, QQ *coeff) {
        QQ t4, t14, t20, t22, t21, t19, t18, t17, t6, t16,
           t2, t7, t15, t13, t3, t1,
           tmp;

        // maple gen code
        t4 = w*w;
        t14 = (1-t4)*w;
        t20 = 2*t14;
        t22 = yc*t20;
        t21 = xc*t20;
        t19 = 4*t14*yc*xc;
        t18 = -4*t4;
        t17 = 4*t4;
        t6 = t4*t4;
        t16 = 1+t6;
        t2 = yc*yc;
        t7 = b*b;
        t15 = t2-t7;
        t13 = -2*t4+t16;
        t3 = a*a;
        t1 = xc*xc;

        coeff[0] = t16*t7+(4*t3-2*t7)*t4;
        coeff[1] = 2*(a-b)*(a+b)*w*(w-1)*(1+w);
        coeff[2] = t7*t17+t13*t3;
        coeff[3] = (xc*t18+t22)*t3+(-t13*xc-t22)*t7;
        coeff[4] = (yc*t18-t21)*t7+(t21-t13*yc)*t3;
        coeff[5] = (t2*t17+t13*t1+t19)*t7+(t15*t6-t19+
                    2*(2*t1-t7-t2)*t4+t15)*t3;

        // BAD IDEA --> spoils invariants!
    //    for (int i = 0; i < 6; i++) coeff[i] = coeff[i] / ((1+w*w)*(1+w*w));
    }

    void initialize_coefficients(const QQ& a_, const QQ& b_, 
                        const QQ& w_, const QQ& xc_, const QQ& yc_) {
        // set parametric coeffs
        pc_[0] = a_; pc_[1] = b_; pc_[2] = w_; pc_[3] = xc_; pc_[4] = yc_;
#ifdef CGAL_VORELL_CIRCLIFY
        if (a_ > b_) {
            //pc_[0] = b_ + (a_ - b_) / QQ(10);
            pc_[0] = b_;
            pc_[1] = pc_[0];
        } else if (a_ < b_) {
            //pc_[0] = a_ + (b_ - a_) / QQ(10);
            pc_[0] = a_;
            pc_[1] = pc_[0];
        }
#endif
        
        // circle?
        if (pc_[0] == pc_[1]) {
            pc_[2] = 0;
            is_circle_ = true;
        } else {
          is_circle_ = false;
          if (CGAL::abs(pc_[0]) < CGAL::abs(pc_[1])) { // "major" < "minor" !!
              // we rotate parametrization by PI/2 
              // this can be done CW by -(1-w)/(1+w)  (each point moves CW)
              // or CCW by (1+w)/(1-w) 
              pc_[0] = b_;
              pc_[1] = a_;
              if (w_ != 1) pc_[2] = (1+w_) / (1-w_);
              else pc_[2] = 0;
          }
        }

        // set algebraic coeffs
        algebraic_coefficients(pc_[0], pc_[1], pc_[2], 
                                pc_[3], pc_[4], &ac_[0]);

        // (pre)compute invariants
        J1_ = ac_[0] + ac_[2];
        J2_ = ac_[0]*ac_[2] - ac_[1]*ac_[1];
    }

    template< class Stream> void point_separator(Stream &W) const { }
    template< class Stream> void object_separator(Stream &W) const { }
    
    void point_separator(std::ostream &W) const { W << std::endl; }
    void object_separator(std::ostream &W) const { point_separator(W); }

public:

    Ellipse_2() {
        res = 0;
        pc_[0] = pc_[1] = pc_[2] = pc_[3] = pc_[4] = QQ(0);
        ac_[0] = ac_[1] = ac_[2] = ac_[3] = ac_[4] = ac_[5] = QQ(0);
        J1_ = J2_ = QQ(0);
        is_circle_ = true;
        eid = 0;
    }

    Ellipse_2(const QQ& a_, const QQ& b_, const QQ& w_, 
                    const QQ& xc_, const QQ& yc_,
                    const int res_ = ELLIPSE_DEFAULT_RES) { 
        initialize_coefficients(a_, b_, w_, xc_, yc_);
        res = res_;
        OMP_guard<OMP_lock> x_(mutex);
        eid = ++eid_count;
    }

    bool is_circle() const { return is_circle_; }
    
    inline int get_id() const { return eid; }
    
    inline const QQ& parametric_coefficient(const int i) const { return pc_[i]; }
    inline const QQ& major_axis() const { return CGAL::abs(pc_[0]); }
    inline const QQ& minor_axis() const { return CGAL::abs(pc_[1]); }
    inline const QQ& rotation() const { return pc_[2]; }
    inline const QQ& x_center() const { return pc_[3]; }
    inline const QQ& y_center() const { return pc_[4]; }
    inline const QQ& algebraic_coefficient(const int i) const { return ac_[i]; }
    inline const QQ& invariant_J1() const { return J1_; }
    inline const QQ& invariant_J2() const { return J2_; }
    
    const typename ET::Kernel::Point_2 point() const { 
        return typename ET::Kernel::Point_2(x_center(), y_center());
    }
    
    QQ evaluate_equation(const QQ& x, const QQ& y) const {
        return ac_[0]*x*x + ac_[5] + ac_[2]*y*y + 2*ac_[4]*y + 
               2*(ac_[3] + ac_[1]*y)*x;
    }

    inline Bounded_side bounded_side(const QQ& x, const QQ& y) const {
        return static_cast<Bounded_side>(-CGAL::sign(evaluate_equation(x, y)));
    }
    
    void translate(const QQ& x, const QQ& y) {
        pc_[3] = x; pc_[4] = y;
        initialize_coefficients(pc_[0], pc_[1], pc_[2], pc_[3], pc_[4]);
    }
    
    inline QQ boundary_x(const QQ& t) const {
        return boundary_x(pc_[0], pc_[1], pc_[2], pc_[3], pc_[4], t);
    }

    inline QQ boundary_y(const QQ &t) const {
        return boundary_y(pc_[0], pc_[1], pc_[2], pc_[3], pc_[4], t);
    }

    // i-point coordinates
    inline QQ boundary_x_inf() const {
        return boundary_x(-pc_[0], -pc_[1], pc_[2], pc_[3], pc_[4], 0);
    }

    inline QQ boundary_y_inf() const {
        return boundary_y(-pc_[0], -pc_[1], pc_[2], pc_[3], pc_[4], 0);
    }

//    inline QQ evolute_x(const QQ& t) const {
//        return evolute_x(pc_[0], pc_[1], pc_[2], pc_[3], pc_[4], t);
//    }

//    inline QQ evolute_y(const QQ &t) const {
//        return evolute_y(pc_[0], pc_[1], pc_[2], pc_[3], pc_[4], t);
//    }

    QQ medial_x(const QQ& t) const {
        if (is_circle()) return x_center();
        if (major_axis() > minor_axis())
            return major_medial_x(pc_[0], pc_[1], pc_[2], pc_[3], pc_[4], t);
        else
            return minor_medial_x(pc_[0], pc_[1], pc_[2], pc_[3], pc_[4], t);
    }

    QQ medial_y(const QQ& t) const {
        if (is_circle()) return y_center();
        if (major_axis() > minor_axis())
            return major_medial_y(pc_[0], pc_[1], pc_[2], pc_[3], pc_[4], t);
        else
            return minor_medial_y(pc_[0], pc_[1], pc_[2], pc_[3], pc_[4], t);
    }

    inline QQ symmetric(const QQ& t) const {
        if (major_axis() >= minor_axis()) return -t;
        else return 1/t;
    }

    Bounded_side boundary_relpos(const Ellipse_2 &q, const Root& t) const {
        typedef typename ET::AK::Sign_at_1 Sign_at;
        upolz_t poly = gencode().disc_polar_poly(ELL_PARAM_COEFFS(q),
                                                 pc_[0], pc_[1], pc_[2], pc_[3], pc_[4]);
        return static_cast<Bounded_side>(-Sign_at()(poly, t));
    }
    
    QQ squared_distance_from_edges(const QQ& x, const QQ& y) const {
        QQ d0 = CGAL::square(x - boundary_x(-1)) + CGAL::square(y - boundary_y(-1));
        QQ d1 = CGAL::square(x - boundary_x(0)) + CGAL::square(y - boundary_y(0));
        QQ d2 = CGAL::square(x - boundary_x(1)) + CGAL::square(y - boundary_y(1));
        QQ d3 = CGAL::square(x - boundary_x_inf()) + CGAL::square(y - boundary_y_inf());
        QQ m1 = CGAL::min(d1, d3);
        QQ m2 = CGAL::min(d0, d2);
        return CGAL::min(m1, m2);
    }

    bool is_on_axes(const QQ& x, const QQ& y) const {
        return CGAL::determinant(x, y, QQ(1), pc_[3], pc_[4], QQ(1), boundary_x(0), boundary_y(0), QQ(1)) == 0 ||
               CGAL::determinant(x, y, QQ(1), pc_[3], pc_[4], QQ(1), boundary_x(1), boundary_y(1), QQ(1)) == 0;
    }

//    int ccw(QQ x1, QQ y1, QQ x2, QQ y2, Root t) const;
    
//    int hash() const;
    
    void print_equation(std::ostream& o) const {
        o << '(' << ac_[0] << ")*x^2 + " <<
             '(' << 2*ac_[1] << ")*x*y + " <<
             '(' << ac_[2] << ")*y^2 + " <<
             '(' << 2*ac_[3] << ")*x + " <<
             '(' << 2*ac_[4] << ")*y + " <<
             '(' << ac_[5] << ')';
    }
    
    bool operator==(const Ellipse_2 &x) const {
        return (eid == x.eid) || (pc_[0] == x.pc_[0] && pc_[1] == x.pc_[1] &&
               pc_[2] == x.pc_[2] && pc_[3] == x.pc_[3] && 
               pc_[4] == x.pc_[4]);
    }
    
    bool operator<(const Ellipse_2 &x) const {
        if (pc_[0] < x.pc_[0]) return true;
        if (pc_[0] > x.pc_[0]) return false;
        if (pc_[1] < x.pc_[1]) return true;
        if (pc_[1] > x.pc_[1]) return false; 
        if (pc_[2] < x.pc_[2]) return true;
        if (pc_[2] > x.pc_[2]) return false; 
        if (pc_[3] < x.pc_[3]) return true;
        if (pc_[3] > x.pc_[3]) return false; 
        if (pc_[4] < x.pc_[4]) return true;
        return false;
    }
    
    template< class Stream > void draw(Stream &w) const {
        QQ xc = x_center();
        QQ yc = y_center();

        QQ t1 = -1;
        QQ t2 = 1;
        QQ tstep = (t2-t1)/(res/2);
        QQ t;

        for (t = t1; t < t2; t += tstep) {
            Point_2 p(CGAL::to_double(boundary_x(t)), 
                      CGAL::to_double(boundary_y(t)));
            w << p;
            point_separator(w);
        }
        // symmetric part
        for (t = t1; t < t2; t += tstep) {
            Point_2 p(CGAL::to_double(2*x_center() - boundary_x(t)), 
                      CGAL::to_double(2*y_center() - boundary_y(t)));
            w << p;
            point_separator(w);
        }
        // to close the i-point hole
//            t = t1;
//        x = e.get_x_coord(t);
//        y = e.get_y_coord(t);
//        cout << CGAL::to_double(x) << ' ' << CGAL::to_double(y) << endl;
        object_separator(w);
    }
    
    int get_resolution() const { return res; }
    int set_resolution(int res_) { return res = res_; }
};

template<class ET>
int Ellipse_2<ET>::eid_count = 0;

template<class ET>
OMP_lock Ellipse_2<ET>::mutex;

//template< class Stream, class ET >
//inline Stream& operator<<(Stream& s, const Ellipse_2<ET> &e) {
//    e.draw(s);
//    return s;
//}

#undef ELLIPSE_DEFAULT_RES

} //namespace CGAL
#endif // CGAL_ELLIPSE_2_H
