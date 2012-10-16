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

#ifndef CGAL_BITANGENT_H
#define CGAL_BITANGENT_H

//#include<signal.h>

#include<iostream>
#include<CGAL/basic.h>

#include<CGAL/Ellipse_2.h>
#include<CGAL/Range.h>
#include<CGAL/vorell.h>
#include<CGAL/Voronoi_diagram_of_ellipses_2/Generated_code.h>

#include<CGAL/Object_cache.h>

namespace CGAL {

namespace VORELL {

enum Relative_position { HIDING = -1, HIDDEN = 0, 
                         PSEUDO_CIRCLES = 2, SEPARATED = 4, NOT_PC = 8 };
// HIDDEN/HIDING are also internally tangent
// SEPARATED is also externally tangent
    
template<class ET>
class Bitangent {
public:
    typedef typename ET::QQ QQ;
    typedef typename ET::AK AK;
    typedef typename ET::upolz_t upolz_t;
    typedef typename ET::upolq_t upolq_t;
    typedef typename ET::Root Root;

    typedef std::pair<int, int> Key_type;
    inline Key_type key() const { return std::make_pair(e1_.get_id(), e2_.get_id()); }

private:
    typedef typename CGAL::VORELL::Generated_code<ET> gencode;
    typedef typename AK::Sign_at_1 Sign_at;
    typedef typename AK::Bound_between_1 Bound_between;
    typedef VORELL::Range<Root> Range;

    Object_cache<Bitangent> cache;
    
    upolz_t poly_;  // tangency points of bitangent
    std::vector<Root> sol;
    enum Bitangent_index {  // sol index
            EXTERNAL, INTERNAL, OTHER_INTERNAL, OTHER_EXTERNAL };
        
    Ellipse_2<ET> e1_, e2_;
    Relative_position relpos;
    bool degenerate; // int/ext tangent
    Root hiding_degeneracy; // common tan. point when e2 hiding e1
    
    void set_sol(const Root& e1, const Root& i1, 
                 const Root& i2, const Root& e2) {
        sol.reserve(4);
        sol.push_back( e1 ); 
        sol.push_back( i1 );
        sol.push_back( i2 );
        sol.push_back( e2 );
    }
    
    void setup_degenerate_pseudocircles(std::vector<Root>& tsol);
    void setup_pseudocircles(std::vector<Root>& tsol, std::vector<Root>& isol);
    void initialize();

public:
    Bitangent(): relpos(HIDDEN) { }

    Bitangent(const Ellipse_2<ET> &e1, const Ellipse_2<ET> &e2): 
            e1_(e1), e2_(e2) {

        const Bitangent *h = cache.read(*this);
        if (h) {
            *this = *h;
            return;
        }

        poly_ = gencode().btan_poly(ELL_PARAM_COEFFS(e1_), 
                                    ELL_PARAM_COEFFS(e2_));
        
        // TODO: check
        if ( typename ET::PTZ::Degree()(poly_) == 0 ) {
            degenerate = 0;
            // coinciding ellipses
            return;
        }
        
        CGAL_assertion( typename ET::PTZ::Degree()(poly_) == 4 );
        
        // TODO: i-points
        
        initialize();

        cache.write(*this);
    }
    
    inline Relative_position relative_position() const { return relpos; }
    inline bool is_degenerate_pair() const { return degenerate; }
    inline Root get_hiding_degeneracy() const { return hiding_degeneracy; }
    
    inline upolz_t poly() const { return poly_; }
    
//       returns the external bitangent line that leaves both ellipses on
//       the left side;
//       if SEPARATED: the other 3 bitangents follow in cyclic order ccw:
//       int, int, ext
    inline Root external() const { 
        CGAL_assertion( relpos >= PSEUDO_CIRCLES );
        return sol[EXTERNAL]; 
    }
    
    // returns the external bitangent line that leaves both ellipses on
    // the right side;
    inline Root other_external() const { 
        CGAL_assertion( relpos >= PSEUDO_CIRCLES );
        return sol[OTHER_EXTERNAL]; 
    }

    inline Root internal() const { 
        CGAL_precondition( relpos >= PSEUDO_CIRCLES );
        return sol[INTERNAL]; 
    }
    
    inline Root other_internal() const { 
        CGAL_precondition( relpos >= PSEUDO_CIRCLES );
        return sol[OTHER_INTERNAL]; 
    }

    // internal bitangent range or intersection range
    inline Range internal_range() const {
        CGAL_precondition( relpos >= PSEUDO_CIRCLES );
        return Range(internal(), other_internal()); 
    }

    inline Range first_range() const {
        CGAL_precondition( relpos >= PSEUDO_CIRCLES );
        return Range(external(), internal()); 
    }

    inline Range last_range() const {
        CGAL_precondition( relpos >= PSEUDO_CIRCLES );
        return Range(other_internal(), other_external());
    }

    inline Range CH_range() const {
        CGAL_precondition( relpos >= PSEUDO_CIRCLES );
        return Range(external(), other_external());
    }

    // relative position =  1 if intersects bitangent (ON_BOUNDED_SIDE), 
    //                      0 if tangent (int/ext) to bitangent (ON_BOUNDARY), 
    //                     -1 otherwise (ON_UNBOUNDED_SIDE)
    Bounded_side relative_position_of_ellipse(
                                        const Bitangent<ET>& bt13) const {
        if (bt13.relative_position() == HIDING) {
            if (compare(external(), bt13.get_hiding_degeneracy()) == EQUAL) {
                return ON_UNBOUNDED_SIDE;
            } else {
                return ON_BOUNDED_SIDE;
            }
        }
        if (bt13.relative_position() == HIDDEN) return ON_UNBOUNDED_SIDE;
        return static_cast<Bounded_side>(Sign_at()(bt13.poly(), external()));
    }
    
    const Ellipse_2<ET>& get_e1() const { return e1_; }
    const Ellipse_2<ET>& get_e2() const { return e2_; }
};

// degenerate pseudo-circles (int. tangent)
template<class ET> 
void Bitangent<ET>::setup_degenerate_pseudocircles(std::vector<Root>& tsol)
{
    CGAL_precondition( tsol.size() == 3 );
    
    std::vector<Root> isol;
    // intersection points of ellipses
    upolz_t ipoly_ = gencode().inter_poly(ELL_PARAM_COEFFS(e1_), 
                                          ELL_PARAM_COEFFS(e2_));
    isol = ET::Real_roots(ipoly_);
    CGAL_assertion( isol.size() == 3 );
    
    int di; // denegerate intersection point
    int dt; // denegerate bitangent point

    for (di = 0; di < isol.size(); di++) {
        for (dt = 0; dt < tsol.size(); dt++) 
            if (isol[di] == tsol[dt]) break;
        if (dt < tsol.size()) break;
    }
    set_sol(tsol[(dt+1)%3], isol[(di+1)%3], isol[(di+2)%3], tsol[(dt+2)%3]);
}

// pseudo-circles
template<class ET> 
void Bitangent<ET>::setup_pseudocircles(std::vector<Root>& tsol, 
                                        std::vector<Root>& isol)
{
    CGAL_precondition( tsol.size() == 2 );
    CGAL_precondition( isol.size() == 2 );

    QQ med = Bound_between()(isol[0], isol[1]);
    
    if (e2_.bounded_side(e1_.boundary_x(med), 
                         e1_.boundary_y(med)) ==  ON_BOUNDED_SIDE) {
        if (Range(tsol[0], tsol[1]).contains(isol[0]))
            set_sol(tsol[0], isol[0], isol[1], tsol[1]);
        else set_sol(tsol[1], isol[0], isol[1], tsol[0]);
    } else {
        if (Range(tsol[0], tsol[1]).contains(isol[0]))
            set_sol(tsol[0], isol[1], isol[0], tsol[1]);
        else set_sol(tsol[1], isol[1], isol[0], tsol[0]);
    }
}

template<class ET> 
void Bitangent<ET>::initialize()
{
    std::vector<Root> tsol;
    tsol = ET::Real_roots(poly_);
    // TODO/CHECK: handle multiplicities? (due to interface change)

    //relpos = static_cast<Relative_position>(tsol.size());
    // TODO: discriminate between twofold internal tangency and P-circles
    // TODO: discriminate between external tangency and P-circles with
    // internal tangency, etc.
    
    //if (relpos < PSEUDO_CIRCLES) return;
    
    upolz_t tan11 = gencode().tan_poly_xy(ELL_PARAM_COEFFS(e1_), 
                                            e1_.x_center(), e1_.y_center());
    upolz_t tan12 = gencode().tan_poly_xy(ELL_PARAM_COEFFS(e1_), 
                                            e2_.x_center(), e2_.y_center());
    
    CGAL_assertion_code( for (int i = 0; i < tsol.size(); i++) { );
    CGAL_assertion( Sign_at()(tan11, tsol[i]) == CGAL::POSITIVE );
    CGAL_assertion_code( } );
    
    Sign btype[4];
    int num_ext = 0, num_int = 0;
    for (int i = 0; i < tsol.size(); i++) {
        btype[i] = Sign_at()(tan12, tsol[i]);
        if (btype[i] > 0) num_ext++;
        if (btype[i] < 0) num_int++;
    }
    
    CGAL_assertion_code( for (int i = 0; i < tsol.size(); i++) { );
    CGAL_assertion( btype[i] != CGAL::ZERO );
    CGAL_assertion_code( } );
    // TODO: doe we need degenerate ellipses? like a segment or point ?
    
    
    switch (tsol.size()) {
    case 4:
        if( num_ext == 2 && num_int == 2) {
            relpos = SEPARATED;
            degenerate = false;
            if (btype[0] == POSITIVE && btype[1] == NEGATIVE) {
            // [1,-1]
                sol = tsol;
            } else if (btype[0] == NEGATIVE && btype[1] == NEGATIVE) {
                // [-1,-1]
                set_sol(tsol[3], tsol[0], tsol[1], tsol[2]);
            } else if (btype[0] == NEGATIVE && btype[1] == POSITIVE) {
                // [-1,1]
                set_sol(tsol[2], tsol[3], tsol[0], tsol[1]);
            } else {
                // [1,1]
                set_sol(tsol[1], tsol[2], tsol[3], tsol[0]);
            }
        } else {
          relpos = NOT_PC;
          degenerate = false;
          sol = tsol;
        }
        break;

    case 3:
        if (num_ext == 2 && num_int == 1 ) {
            relpos = SEPARATED;
            degenerate = true; // externally tangent
            if (btype[0] == POSITIVE && btype[1] == NEGATIVE) {
            // [1,-1]
                set_sol(tsol[0], tsol[1], tsol[1], tsol[2]);
            } else if (btype[0] == CGAL::NEGATIVE && btype[1] == CGAL::POSITIVE) {
                // [-1,1]
                set_sol(tsol[2], tsol[0], tsol[0], tsol[1]);
            } else {
                // [1,1]
                set_sol(tsol[1], tsol[2], tsol[2], tsol[0]);
            }
        } else if (num_ext == 3) {
            relpos = PSEUDO_CIRCLES; 
            degenerate = true; // internally tangent
            // TODO: hidden?
            setup_degenerate_pseudocircles(tsol);
        } else {
            CGAL_assertion( false );
            degenerate = true;
        }
        break;
    case 2: {
            std::vector<Root> isol;
            // intersection points of ellipses
            upolz_t ipoly_ = gencode().inter_poly(ELL_PARAM_COEFFS(e1_), 
                                                  ELL_PARAM_COEFFS(e2_));
            CGAL_assertion( typename ET::PTZ::Degree()(ipoly_) == 4 );
            // TODO: i-points
            
            isol = ET::Real_roots(ipoly_);
            CGAL_assertion( isol.size() == 2 );
            if (isol[0] == tsol[0] && isol[1] == tsol[1]) {
                relpos = HIDDEN; 
                // TODO: hidden vs hiding
                degenerate = true; // internally tangent
            } else {
                relpos = PSEUDO_CIRCLES;
                degenerate = false;
                setup_pseudocircles(tsol, isol);
            }
        }
        break;
    case 1: {
            int root = 0;
            // TODO: trigger this test to see if it really works
            if (compare(tsol[0], Root(0)) == EQUAL) root = 1;
            if (e2_.bounded_side(e1_.boundary_x(root), e1_.boundary_y(root)) == 
                                                          ON_UNBOUNDED_SIDE) {
                // e2 is hidden inside e1
                relpos = HIDDEN;
            } else {
                // e2 is hiding e1 (encloses it)
                relpos = HIDING;
                hiding_degeneracy = tsol[0];
            }
            degenerate = true;
            CGAL_assertion( btype[0] == POSITIVE );
        }
        break;
    case 0:
        if (e2_.bounded_side(e1_.boundary_x(0), e1_.boundary_y(0)) == 
                                                      ON_UNBOUNDED_SIDE) {
            // e2 is hidden inside e1
            relpos = HIDDEN;
        } else {
                // e2 is hiding e1 (encloses it)
            relpos = HIDING;
        }
        degenerate = false;
        break;
    }
}


} //namespace
} //namespace CGAL
#endif

