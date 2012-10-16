//    (c) 2007-2009 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED

#ifndef CGAL_RANGE_H
#define CGAL_RANGE_H

#include<CGAL/basic.h>

namespace CGAL {
namespace VORELL {

template<class T>
class Range {
protected:
    T a, b;
    bool empty;

public:
    enum { PRUNE_NONE, PRUNE_LEFT, PRUNE_RIGHT };

    Range() { empty = true; }
    Range(const T &a_, const T &b_) { empty = false; a = a_;  b = b_; CGAL::compare(a, b); }
//    Range(const IntF a_): empty(false), a(a_.lower()), b(a_.upper()) { }
        
    const T& left() const { return a; }
    const T& right() const { return b; }
    bool is_empty() const { return empty; }
    void mkempty() { empty = true; }
    bool is_finite() const { return empty || CGAL::compare(a, b) <= 0; }
    bool is_infinite() const { return !is_finite(); }
    bool is_exact() const { return empty || CGAL::compare(a, b) == CGAL::EQUAL; }

    T width() const { 
        CGAL_assertion (is_finite());
        if (is_empty()) return T(0);
        else return (b - a);
    }
    
    friend std::ostream& operator<<(std::ostream& i, const Range<T>& r) {
        if (r.empty) return (i << "[ ]");
        else return (i << '[' << CGAL::to_double(r.a) << ',' << CGAL::to_double(r.b) << ']');
    }

    Range complement() const { return Range<T>(b, a); }
    Range symmetric() const { return Range<T>(T(-1)/a, T(-1)/b); }
    
    Range hull(const Range<T> &r) const {
        CGAL_precondition (is_finite() && r.is_finite());
        return Range(std::min(a, r.a), std::max(b, r.b));
    }

    template<typename R>
    bool contains(const R q) const {
        if (empty) return false;
        if (is_finite()) return (a <= q &&  b >= q);
        else return (a <= q || b >= q);
    }

    template<typename R>
    bool strictly_contains(const R q) const {
        if (empty) return false;
        if (is_finite())
            return (a < q && b > q);
        else return (a < q || b > q);
    }


    bool subset_of(const Range& r) const {
        if (empty) return true;
        if (is_finite()) {
            if (r.is_finite())
                return a >= r.a && b <= r.b;
            else
                return a >= r.a || b <= r.b;
        } else {
            if (r.is_finite())
                return false;
            else 
                return a >= r.a && b <= r.b;
        }
    }

    bool strict_subset_of(const Range& r) const {
        if (empty) return true;
        if (is_finite()) {
            if (r.is_finite())
                return a > r.a && b < r.b;
            else
                return a > r.a || b < r.b;
        } else {
            if (r.is_finite())
                return false;
            else 
                return a > r.a && b < r.b;
        }
    }

    bool touches(const T q) const {
        if (empty) return false;
        return q == a || q == b;
    }
    
    int order(const T x1, const T x2) const {
        // std::cerr << "this = " << *this << x1 << ' ' << x2 << std::endl;
        // std::cerr << "this = [" << a << ',' << b << ']' << x1 << ' ' << x2 << std::endl;
        assert (contains(x1) && contains(x2));
        if (is_finite()) {
            if (x1 < x2) return -1;
            if (x1 > x2) return 1;
            if (x1 == x2) return 0;
        } else {
            if (x2 <= b) {
                if (x1 < x2) return -1;
                if (x1 > x2 && x1 <= b) return 1;
                if (x1 > x2 && x1 >= a) return -1;
                if (x1 == x2) return 0;
            } else {
                if (x1 > x2) return 1;
                if (x1 < x2 && x1 >= a) return -1;
                if (x1 < x2 && x1 <= b) return 1;
                if (x1 == x2) return 0;
            }
        }
        assert (false);
    }

//  *********************************************
//     compute the intersection of two ranges
//     the result is a range
//     in case of multiple ranges, only one is returned
//  *********************************************
    Range intersection(const Range &r) const {
        if (is_finite()) {
            if (r.is_finite()) {
                T c1 = std::max(a,r.a); 
                T c2 = std::min(b,r.b);
                if (c1 > c2) return Range();
                else return Range(c1,c2);
            } else {
                T c1 = std::max(a,r.a); 
                T c2 = b;
                Range r1;
                if (c1 <= c2) r1 = Range(c1,c2);
                c2 = std::min(b,r.b); c1 = a;
                Range r2;
                if (c1 <= c2) r2 = Range(c1,c2);
                if (r1.is_empty()) {
                        if (r2.is_empty()) return Range();
                        else return r2;
                } else {
                        if (r2.is_empty()) return r1;
                        else {
#if VERBOSE > 0
                            std::cerr << "intersection splits range\n";
                            std::cerr << r1 << " AND " << r2 << std::endl;
#endif
                            // the intersection of two ranges should contain the left endpoint of the second
                            if (r1.touches(r.a)) return r1; // since r is infinite!
                            else return r2;
                        }
                }
            }
        } else {
            if (r.is_finite()) {
                T c1 = std::max(r.a,a); 
                T c2 = r.b;
                Range r1;
                if (c1<=c2) r1 = Range(c1,c2);
                c2 = std::min(r.b,b); c1 = r.a;
                Range r2;
                if (c1<=c2) r2 = Range(c1,c2);
                if (r1.is_empty()) {
                        if (r2.is_empty()) return Range();
                        else return r2;
                } else {
                        if (r2.is_empty()) return r1;
                        else { 
#if VERBOSE > 0
                            std::cerr << "intersection splits range\n";
                            std::cerr << r1 << " AND " << r2 << std::endl;
#endif
                            if (r1.touches(r.a)) return r1;
                            else return r2;
                        }
                }
            } else {
                Range c(complement());
                Range cr(r.complement());
                Range rf, ri;
                if (c.intersection(cr).is_empty()) { // finite part
                    if (cr.b < c.a) rf = Range(cr.b, c.a);
                    else rf = Range(c.b, cr.a);
                }
                // infinite part
                T c1 = std::max(a,r.a); 
                T c2 = std::min(b,r.b);
                ri = Range(c1,c2);
#if VERBOSE > 0
                if (!rf.is_empty()) {
                    std::cerr << "intersection splits range\n";
                    std::cerr << rf << " AND " << ri << std::endl;
                }
#endif
                if (rf.touches(r.a)) return rf;
                else return ri;
            }
        }
    }

    // prune (reject) points in r
    // NOTE: does not prune if r is contained entirely in *this
    Range prune(const Range &r, int alg = PRUNE_NONE) const {
        bool v1 = r.contains(a);
        bool v2 = r.contains(b);
#if VERBOSE > 0
        std::cerr << (*this) << " prune " << r << std::endl;
//        std::cerr << "v1 = " << v1 << " v2 = " << v2 << std::endl;
#endif
        if (v1 && v2) return Range();
        if (v1) return Range(r.b, b);
        if (v2) return Range(a, r.a);
#if VERBOSE > 0
        std::cerr << "pruning splits range\n"; // TODO: unless both endpoints left/right
#endif
        if (alg == PRUNE_LEFT) return Range(r.b, b);
        if (alg == PRUNE_RIGHT) return Range(a, r.a);
        return Range(a, b);
    }
};

} // namespace

} //namespace CGAL
#endif
