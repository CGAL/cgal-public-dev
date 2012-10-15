// ============================================================================
//
// Copyright (c) 2001-2006 Algorithmic Solutions Software GmbH (Germany)
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// for the license conditions see the file LICENSE.AS in the EXACUS/ 
// root directory.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : SoX
//
// File          : include/SoX/sweep_curves.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     :  Eric Berberich <eric@mpi-inf.mpg.de>
//                  Susan Hert <hert@mpi-inf.mpg.de>
//                  Kurt Mehlhorn <mehlhorn@mpi-inf.mpg.de>
//                  Stefan Naeher <stefan@mpi-inf.mpg.de>
//                  Michael Seel <seel@mpi-inf.mpg.de>
//
// ============================================================================

/*!\file SoX/sweep_curves.h
 * \brief 
 * Sweep-line function for computing an arrangement of planar curve segments.
 */

/*!\defgroup SoX_sweeping Sweep-line function for arrangement computation
 */

#ifndef SoX_SWEEP_CURVES_H
#define SoX_SWEEP_CURVES_H

#include <SoX/basic.h>

#ifdef SoX_SHOW_TRACE
#define SoX_TRACE(stmt) stmt
#else
#define SoX_TRACE(stmt)
#endif

#ifdef SoX_SHOW_SWEEP_LINE
#define SoX_SWEEP_LINE(stmt) stmt
#else
#define SoX_SWEEP_LINE(stmt)
#endif

#ifdef SoX_CHECK_VALIDITY
#define SoX_VALIDITY_CHECK(stmt) CGAL_assertion(stmt)
#else
#define SoX_VALIDITY_CHECK(stmt)
#endif

#include <LEDA/sortseq.h>
#include <LEDA/graph.h>
#include <LEDA/p_queue.h>
#include <LEDA/map.h>
#include <LEDA/map2.h>
#include <LEDA/std/math.h>
#include <SoX/basic.h>
#include <SoX/LEDA_compare.h>
#include <list>
#include <iostream>
#include <utility>
#include <vector>

#if SoX_LINEAR_REORDER
#include <map>
#include <SoX/Reorder_tree.h>
#endif

#ifdef SoX_SWEEP_AB_TREE
#include <LEDA/impl/ab_tree.h>
#endif // SoX_SWEEP_AB_TREE

// TODO replace in an future implementation the global function ID_Number 
// on points and segments with a memory adress comparison on the objects.
// Currently, we need ID_Number to use LEDA datastructures internally.
// OR NOT - only if one wishes to exchange the datastructures

namespace SoX {

#ifdef SoX_SWEEP_AB_TREE
namespace Intern {

class ab_tree_bugfix : public leda::ab_tree {
public:
    leda::ab_tree_node* min_item() const { return min(); }
    leda::ab_tree_node* max_item() const { return max(); }
};
} // namespace Intern
typedef Intern::ab_tree_bugfix Sortseq_impl;
typedef Sortseq_impl::item Seq_item;

#else

/*!\brief 
 * leda::sortseq implementation chosen for sweep
 *
 * The underlying data structure for the \c leda::sortseq number type
 * (its third template argument) for use within \c sweep_curves()
 * and friends for X- and Y-structure can be  changed by editing
 * \c SoX/sweep_curces.h and modifying this typedef.
 * Don't forget to include the appropriate \c LEDA/impl/ header.
 */
typedef leda::skiplist Sortseq_impl;
//  \brief LEDA item type matching \c Sortseq_impl (set automatically)
typedef leda::seq_item Seq_item;
#endif // SoX_SWEEP_AB_TREE


namespace Intern {

/* \brief
 * stores lexicographically ordered intersection points for a pair of segments
 * and their x-structure items.
*/
template <class CurvePoint_2>
struct Intersect_info {
    typedef std::pair<CurvePoint_2, Seq_item>  Intersect_pair;
    
    int num_to_right;
    std::list<Intersect_pair>  i_points; 
    
    Intersect_info() : num_to_right(-1) {}
};



/* \brief
 * makes \a v the target of e and appends the twin of \a e (its reversal edge)
 * to \a v's adjacency list
 */
template <class CurveSegment_2, class CurvePoint_2>
void link_as_target_and_append(leda::GRAPH<CurvePoint_2, CurveSegment_2>& G, 
                               leda_node v, leda_edge e,
                               CurveSegment_2& curve) {
    G.assign(e, curve);
    leda_edge e_rev = G.reversal(e);
    G.assign(e_rev, curve);
    G.move_edge(e, G.cyclic_adj_pred(e,G.source(e)), v);
    G.move_edge(e_rev, v, G.target(e_rev));
}

/* \brief 
 * returns a newly created edge inserted before the first
 * edge of the adjacency list of \a v. It also creates a reversal edge
 * whose target is \a v.
 */
template <class CurveSegment_2, class CurvePoint_2>
leda_edge 
new_halfedge_pair_at_source(leda::GRAPH<CurvePoint_2, CurveSegment_2>& G,
                            leda_node v) { 
    
    leda_edge e_res, e_rev, e_first = G.first_adj_edge(v);
    if ( e_first == nil ) {
        e_res = G.new_edge(v,v);
        e_rev = G.new_edge(v,v);
    } else {
        e_res = G.new_edge(e_first,v,LEDA::before);
        e_rev = G.new_edge(e_first,v,LEDA::before);
    }
    G.set_reversal(e_res,e_rev);
    return e_res;
}


/* \brief
 * function object used to determine vertical ordering of segments in the
 * Y structure.
 */
template <class CurveSweepTraits_2>
class sweep_cmp_curve : 
        public leda::leda_cmp_base<typename CurveSweepTraits_2::Segment_2> {

public:
    // this instance's template parameter
    typedef CurveSweepTraits_2 Curve_sweep_traits_2;
    
private:
    // internal type
    typedef Curve_sweep_traits_2 CST;

public:
    // type of Point
    typedef typename CST::Point_2 Point_2;
    // type pf Segment
    typedef typename CST::Segment_2 Segment_2;
    
private:
    Point_2   p_sweep;
    Segment_2 c_bottom;
    Segment_2 c_top;
    Curve_sweep_traits_2 traits_;

public:
    /* \brief
     * constructor; \a bottom and \a top are used as sentinels only they have 
     * no geometric meaning.
    */
    sweep_cmp_curve(const Point_2& p, const Segment_2& bottom,
                    const Segment_2& top,
                    const Curve_sweep_traits_2& traits) : 
        p_sweep(p), 
        c_bottom(bottom), 
        c_top(top),
        traits_(traits) {
    }
    
    // places the sweep point a \c p.
    void set_position(const Point_2& p) { 
        p_sweep = p; 
    }
    
    // \pre \a p_sweep is identical to the left endpoint of either c1 or c2. 
    int operator()(const Segment_2& c1, const Segment_2& c2) const { 
        
        typename CST::Source_2 source(
                traits_.source_2_object()
        );
        typename CST::Target_2 target(
                traits_.target_2_object()
        );
        typename CST::Is_degenerate_2 is_degenerate(
                traits_.is_degenerate_2_object()
        );
#if SoX_LINEAR_REORDER
        typename CST::Do_overlap_2 do_overlap(
                traits_.do_overlap_2_object()
        );
#endif
        typename CST::Compare_y_at_x_2 compare_y_at_x(
                traits_.compare_y_at_x_2_object()
        );
        typename CST::Compare_y_right_of_point_2
            compare_y_right_of_point(
                    traits_.compare_y_right_of_point_2_object()
            );
        
        if (ID_Number(c2) == ID_Number(c_top) || 
            ID_Number(c1) == ID_Number(c_bottom)) {
            return -1;
        }
        if (ID_Number(c1) == ID_Number(c_top) || 
            ID_Number(c2) == ID_Number(c_bottom)) {
            return 1;
        }
        if (ID_Number(c1) == ID_Number(c2)) {
            return 0;
        }
        
        SoX_TRACE(
                draw_segment(W, c1, CGAL::GREEN);
                draw_segment(W, c2, CGAL::RED);
        );
        
        int s = 0;
        bool c1_incident = false;
        if (ID_Number(p_sweep) == ID_Number(source(c1))) {
            s = compare_y_at_x(c2, p_sweep);
            c1_incident = true;
        } else {
#if SoX_LINEAR_REORDER
            CGAL_precondition_msg(
                    ID_Number(p_sweep) == ID_Number(source(c2)),
                    "p_sweep not identical to source of either segment."
            );
            s = -compare_y_at_x(c1, p_sweep);
#else
            if (ID_Number(p_sweep) == ID_Number(source(c2))) {
                s = -compare_y_at_x(c1, p_sweep);
            } else {
                if (do_overlap(c1, c2)) {
                    return (ID_Number(c1) > ID_Number(c2) ? +1 : -1);
                }
#else
                // TODO select correct order 
                return -1;
#endif
            }
        }
        
        SoX_TRACE(
                draw_segment(W, c1, CGAL::BLUE);
                draw_segment(W, c2, CGAL::BLUE);
        );
        
        if (s || is_degenerate(c1) || is_degenerate(c2)) return s;
        
        SoX_TRACE(
                draw_segment(W, c1, CGAL::GREEN);
                draw_segment(W, c2, CGAL::RED);
        );
        if (c1_incident) {
            s = compare_y_right_of_point(c1, c2, source(c1));
        } else {
            s = compare_y_right_of_point(c1, c2, source(c2));
        }
        
        SoX_TRACE(
                draw_segment(W, c1, CGAL::BLUE);
                draw_segment(W, c2, CGAL::BLUE);
        );
        
        // overlapping curves will be ordered by their ID_numbers :
        if (s) {
            return s;
        }  else if (ID_Number(c1) > ID_Number(c2)) {
            return +1;
        } else {
            return -1; // identity has been excluded above
        }
    }
};

/* \brief
 * prints the given Y structure's key and information pairs 
 * (one item per line) to stdout.
*/
template <class YStructure>
void  print_y_structure(const YStructure& Y_structure) {
    Seq_item s;
    std::cout << "Y_structure: " << std::endl;
    forall_items(s, Y_structure) {
        std::cout << s << "   " << Y_structure.key(s) << "   " 
                  << Y_structure.inf(s) << std::endl;
    }
    std::cout << std::endl;
}

/* \brief
 * prints the given X structure's key and information pairs 
 * (one item per line) to stdout
*/
template <class XStructure>
void print_x_structure(const XStructure& X_structure) {
    Seq_item s;
    std::cout << "X_structure: " << std::endl;
    forall_items(s, X_structure) {
        std::cout << s << "   " << X_structure.key(s) << "   " 
                  << X_structure.inf(s) << std::endl;
    }
    std::cout << std::endl;
}

/* \brief
 * tests certain invariants of the X and Y structures.
 * 
 * The following tests are performed:
 * \li if two successive segments in the Y structure overlap, the first
 * should point to the second.
 * \li if two successive segments in the Y structure intersect to the right
 * of the point \a p_sweep there should be an item in the X structure
 * corresponding to this point and it should point to one of the segments
 * in the Y structure
 * \li non-overlapping, non-intersecting segments should not have pointers
 * to each other.
 */
template <class CurveSweepTraits_2>
bool 
structures_are_valid(
        const leda::sortseq<typename CurveSweepTraits_2::Segment_2, Seq_item,
        Sortseq_impl>& Y_structure,
        const leda::sortseq<typename CurveSweepTraits_2::Point_2, Seq_item,
        Sortseq_impl>& X_structure,
        const typename CurveSweepTraits_2::Point_2& p_sweep,
        const CurveSweepTraits_2& traits) {

    typedef CurveSweepTraits_2       Curve_sweep_traits_2;
    
    typedef Curve_sweep_traits_2     CST;

    typedef typename CST::Point_2    Point_2;
    typedef typename CST::Segment_2  Segment_2;
    
    typename CST::Intersect_right_of_point_2 intersect_right_of_point(
            traits.intersect_right_of_point_2_object()
    );
    typename CST::Compare_xy_2 compare_xy(
            traits.compare_xy_2_object()
    );
    typename CST::Equal_y_at_x_2 equal_y_at_x(
            traits.equal_y_at_x_2_object()
    );
    typename CST::Do_overlap_2 do_overlap(
            traits.do_overlap_2_object()
    );
    
    Seq_item sit;
    Seq_item sit_succ;
    Point_2 i_point;
    forall_items(sit, Y_structure) {
        if (sit == Y_structure.min_item()) continue;
        if (sit == Y_structure.max_item()) continue;
        
        sit_succ = Y_structure.succ(sit);
        Segment_2 s = Y_structure.key(sit);
        Segment_2 s_succ = Y_structure.key(sit_succ);
        SoX_TRACE(draw_segment(W, s, CGAL::ORANGE););
        SoX_TRACE(draw_segment(W, s_succ, CGAL::VIOLET););
        
#if 0
        if (sit_succ != Y_structure.max_item() &&
            y_order_s,s_succ, p_sweep) != -1) {
            y_order_print_2(s, s_succ, p_sweep);
            std::cerr << "Adjacent segments in Y structure not "
                      << "ordered properly."
                      << std::endl;
            return false;
        } else 
#endif
            if (sit_succ == Y_structure.max_item()) {
                if (Y_structure.inf(sit) != nil) {
                    std::cerr << "Y structure indicates intersection " 
                              << "or overlap with upper sentinel." 
                              << std::endl;
                    return false;
                }
            } else if (do_overlap(s,s_succ)) {
                if (Y_structure.inf(sit) != sit_succ) {  
                    SoX_TRACE(
                            draw_segment(W, s, CGAL::GREEN);
                            draw_segment(W, s_succ, CGAL::RED);
                    );
                    std::cout << "sit = " << sit 
                              << " sit_succ = " << sit_succ 
                              << std::endl;
                    std::cout << "s = " << s << std::endl;
                    std::cout << "s_succ = " << s_succ << std::endl;
                    std::cerr << "Overlapping curves not recorded " 
                              << "overlapping in Y structure." << std::endl;
                    return false;
                }  
            } else if (intersect_right_of_point(s, s_succ, p_sweep, i_point)) {
                if (Y_structure.inf(sit) == nil) {
                    SoX_TRACE(
                            draw_point(W, i_point, CGAL::RED, true);
                            draw_segment(W, s, CGAL::GREEN);
                            draw_segment(W, s_succ, CGAL::RED);
                    );
                    std::cout << "sit = " << sit 
                              << " sit_succ = " << sit_succ 
                              << std::endl;
                    std::cout << "s = " << s << std::endl;
                    std::cout << "s_succ = " << s_succ << std::endl;
                    std::cerr << "Intersection right of p_sweep not "
                              << "recorded in Y structure." << std::endl;
                    return false;
                }
                Seq_item xit = Y_structure.inf(sit);
                if (compare_xy(X_structure.key(xit), i_point)) {
                    SoX_TRACE(
                            draw_point(W, X_structure.key(xit), 
                                       CGAL::GREEN, true);
                            draw_point(W, i_point, CGAL::RED, true);
                            draw_segment(W, s, CGAL::GREEN);
                            draw_segment(W, s_succ, CGAL::RED);
                    );
                    std::cout << "sit = " << sit 
                              << " sit_succ = " << sit_succ 
                              << std::endl;
                    std::cout << "s = " << s << std::endl;
                    std::cout << "s_succ = " << s_succ << std::endl;
                    std::cerr << "Wrong intersection point from "
                              << "X structure recorded in Y structure" 
                              << std::endl;
                    return false;
                }
                if (!equal_y_at_x(Y_structure.key(X_structure.inf(xit)),
                                  X_structure.key(xit))) {
                    
                    SoX_TRACE(
                            draw_point(W, i_point, CGAL::RED, true);
                            draw_segment(W, s, CGAL::GREEN);
                            draw_segment(W, s_succ, CGAL::RED);
                            draw_point(W, X_structure.key(xit), 
                                       CGAL::GREEN, true);
                            draw_segment(W, 
                                         Y_structure.key(
                                                 X_structure.inf(xit)
                                         ), 
                                         CGAL::RED
                            );
                    );
                    std::cout << "sit = " << sit 
                              << " sit_succ = " << sit_succ 
                              << std::endl;
                    std::cout << "s = " << s << std::endl;
                    std::cout << "s_succ = " << s_succ << std::endl;
                    std::cerr << "X structure records wrong Y structure " 
                              << "item" << std::endl;
                    return false;
                }          
            } else if (Y_structure.inf(sit) != nil) {
                // known not to overlap or intersect
                SoX_TRACE(
                        draw_segment(W, s, CGAL::GREEN);
                        draw_segment(W, s_succ, CGAL::RED);
                );
                std::cout << "sit = " << sit << " sit_succ = " << sit_succ 
                          << std::endl;
                std::cout << "s = " << s << std::endl;
                std::cout << "s_succ = " << s_succ << std::endl;
                print_y_structure(Y_structure);
                std::cerr << "Nonoverlapping, nonintersecting segments with "
                          << "non-nil Y structure information " << std::endl;
                return false;
            } else if (Y_structure.inf(sit) == Y_structure.succ(sit)) {
                SoX_TRACE(
                        draw_segment(W, s, CGAL::GREEN);
                        draw_segment(W, s_succ, CGAL::RED);
                );
                std::cout << "sit = " << sit << " sit_succ = " << sit_succ 
                          << std::endl;
                std::cout << "s = " << s << std::endl;
                std::cout << "s_succ = " << s_succ << std::endl;
                std::cerr << "Nonoverlapping curves recorded overlapping in"
                          << " Y structure." <<  std::endl;
                return false;
            }
        SoX_TRACE(draw_segment(W, s, CGAL::BLUE););
        SoX_TRACE(draw_segment(W, s_succ, CGAL::BLUE););
    }
    return true;
}

/* \brief
 * computes the first intersection point to the right of \a p_sweep between
 * Y-structure item \a sit0 and its successor.
 *
 * The intersection point is inserted into the X structure and the information
 * associated with \a sit0 in the Y structure is updated to refer to this
 * intersection point. 
 */   
template <class CurveSweepTraits_2>
static void compute_intersection(
        leda::sortseq<typename CurveSweepTraits_2::Point_2, Seq_item,
        Sortseq_impl>& X_structure,
        leda::sortseq<typename CurveSweepTraits_2::Segment_2, Seq_item,
        Sortseq_impl>& Y_structure,
        const typename CurveSweepTraits_2::Point_2& p_sweep,
        Seq_item sit0, const CurveSweepTraits_2& traits) { 
    
    typedef CurveSweepTraits_2       Curve_sweep_traits_2;
    
    typedef Curve_sweep_traits_2     CST;
    
    typedef typename CST::Point_2    Point_2;
    typedef typename CST::Segment_2  Segment_2;

    typename CST::Intersect_right_of_point_2 intersect_right_of_point(
            traits.intersect_right_of_point_2_object()
    );
    
    Seq_item sit1 = Y_structure.succ(sit0);
    
    if (sit0 == Y_structure.min_item() || sit1 == Y_structure.max_item() ) {
        return;
    }
    
    Segment_2  s0   = Y_structure.key(sit0);
    Segment_2  s1   = Y_structure.key(sit1);
    
    SoX_TRACE(draw_segment(W, s0, CGAL::GREEN););
    SoX_TRACE(draw_segment(W, s1, CGAL::RED););
    
    Point_2 i_point;
    
    if (intersect_right_of_point(s0, s1, p_sweep, i_point)) { 
        Seq_item it = X_structure.insert(i_point, sit0);
        Y_structure.change_inf(sit0, it);
    } 
    SoX_TRACE(draw_segment(W, s0, CGAL::BLUE););
    SoX_TRACE(draw_segment(W, s1, CGAL::BLUE););
}


/* \brief
 * computes the intersection points between Y-structure item \a sit0 and its 
 * successor and inserts the points into the X structure and the pair of 
 * segments into \a inter_dic.
 *  
 * If this pair of segments has not been seen before, all intersection points
 * are computed and those right of the sweep point are inserted into the
 * X structure.  If this pair has been seen before, the intersection 
 * information recorded for this pair in \a inter_dic is examined and update
 * to indicate how many points are now to the right of the sweep point.  Then
 * \a sit0's Y-structure information is updated to point to the (new) first
 * intersection point to the right of \a p_sweep.
 */   
template <class CurveSweepTraits_2>
static void compute_intersection(
        leda::sortseq<typename CurveSweepTraits_2::Point_2, Seq_item, 
                      Sortseq_impl>& X_structure,
        leda::sortseq<typename CurveSweepTraits_2::Segment_2, Seq_item, 
                      Sortseq_impl>& Y_structure,
        const typename CurveSweepTraits_2::Point_2& p_sweep,
        leda_map2<
        typename CurveSweepTraits_2::Segment_2, 
        typename CurveSweepTraits_2::Segment_2, 
        Intern::Intersect_info<typename CurveSweepTraits_2::Point_2> >&
        inter_dic,
        Seq_item sit0, const CurveSweepTraits_2& traits) {
    
    typedef CurveSweepTraits_2       Curve_sweep_traits_2;
    typedef Curve_sweep_traits_2     CST;
    
    typedef typename CST::Point_2    Point_2;
    typedef typename CST::Segment_2  Segment_2;
    
    typename CST::Less_xy_2 less_xy(
            traits.less_xy_2_object()
    );
    typename CST::Intersect_2 intersect(
            traits.intersect_2_object()
    );
    
    typedef typename Intern::Intersect_info<Point_2>::Intersect_pair I_pair;
    
    Seq_item sit1 = Y_structure.succ(sit0);
    
    if (sit0 == Y_structure.min_item() || sit1 == Y_structure.max_item()) {
        return;
    }
    
    Segment_2  s0   = Y_structure.key(sit0);
    Segment_2  s1   = Y_structure.key(sit1);
    SoX_TRACE(draw_segment(W, s0, CGAL::GREEN););
    SoX_TRACE(draw_segment(W, s1, CGAL::RED););
    
    std::list<Point_2>  i_points;
    
    bool update_dict = false;
    Intern::Intersect_info<Point_2> inter_info = inter_dic(s0,s1);
    
    if (inter_info.num_to_right < 0) {
        // this pair has not been considered yet;
        // record the intersection points
        
        inter_info.num_to_right = 0;  
        update_dict = true;
        
        intersect(s0, s1, std::back_inserter(i_points));
        
        typename std::list<Point_2>::iterator ip_it;
        // find first point that is to the right of the sweep point
        for (ip_it = i_points.begin(); ip_it != i_points.end(); ip_it++) {
            SoX_TRACE(draw_point(W, *ip_it, CGAL::ORANGE););
            bool is_smaller = less_xy(p_sweep, *ip_it);
            SoX_TRACE(draw_point(W, *ip_it, CGAL::WHITE););
            if (is_smaller)  {
                break;
            }
        }  
        // insert all remaining points (to right of sweep point) in x-structure
        for (; ip_it != i_points.end(); ip_it++) {
            SoX_TRACE(draw_point(W, *ip_it, CGAL::ORANGE););
            Seq_item it = X_structure.insert(*ip_it, sit0);  
            SoX_TRACE(draw_point(W, *ip_it, CGAL::WHITE););
            inter_info.num_to_right++;
            inter_info.i_points.push_back(I_pair(X_structure.key(it), it));
        }
    }
    // check if an intersection point is to the right 
    // of the current sweep point
    if (inter_info.num_to_right > 0) {
        if (ID_Number(p_sweep) == 
            ID_Number(inter_info.i_points.front().first)) { 
            update_dict = true;
            inter_info.num_to_right--;
            inter_info.i_points.pop_front();
            if (inter_info.num_to_right > 0)
                Y_structure.change_inf(sit0, 
                                       inter_info.i_points.front().second);
        } else {
            while (inter_info.i_points.size() > 0 &&
                   !less_xy(p_sweep, inter_info.i_points.front().first)) {
                update_dict = true;
                inter_info.num_to_right--;
                inter_info.i_points.pop_front();
            }
            if (inter_info.num_to_right > 0)
                Y_structure.change_inf(sit0, 
                                       inter_info.i_points.front().second);
        }
    }
    if (update_dict) {
        inter_dic(s0, s1) = inter_info;
        // need to record symmetric entry because order will switch at 
        // intersection points
        if (inter_info.num_to_right > 0) {
            inter_dic(s1, s0) = inter_info;
        } 
    }
    SoX_TRACE(draw_segment(W, s0, CGAL::BLUE););
    SoX_TRACE(draw_segment(W, s1, CGAL::BLUE););
}

} // namespace Intern

/*!\ingroup SoX_sweeping
 * \brief 
 * computes the arrangement of a set of curve segments,
 * returning the result as an embedded graph.
 *
 * This function computes the arrangement of the curve segments in the
 * iterator range [\a first, \a beyond), storing the result in the directed
 * graph \a G.
 * The vertices of \a G are the endpoints and intersection points of
 * the segments. The edges of \a G are directed and come in pairs of
 * twins (reversal edges of each other) such that each pair represents
 * a piece of curve between two vertices. The edges out of every vertex
 * are ordered cyclically according to the counterclockwise order of segments
 * around the point. Each vertex and edge is labelled with the point
 * or segment it stands for.
 *
 * \pre \a G is empty.
 *
 * If \a use_optimization is \c true, a dictionary will be used to
 * store interesction points of curves such that they need to be computed
 * only once.
 *
 * The point and segment data types to be used and the
 * operations on them are given by the traits class passed as last
 * argument. It must be a model of the (new) \c CurveSweepTraits_2 concept.
 *
 * The implementation is derived from the function \c SWEEP_SEGMENTS()
 * described in Chapter 10 of the LEDA book by Mehlhorn and N&auml;her.
 * It has also drawn on ideas from the generic version of this sweep
 * described in Michael Seel's 2001 research report on
 * the Implementation of Planar Nef Polyhedra.
 * The basic algorithm used is the sweep-line algorithm of Bentley
 * and Ottmann (1979).  The asymptotic running time is
 * O(<I>m</I> + (<I>n</I>+<I>s</I>)log <I>n</I>),
 * where \e n is the number of input segments, \e m is the number
 * of edges in the output graph, and \e s is the number of nodes in
 * the graph.
 * Unlike the original LEDA implementation, this function collapses
 * overlapping parts of segments between two vertices into a single
 * pair of twin edges.
 *
 * Internally, the implementation uses LEDA sorted sequences to maintain the
 * so-called X- and Y-structures. The default implementation of sorted
 * sequences uses a randomized data structure (skip lists), leading to
 * a nondeterministic sequence of predicate evaluations. If you prefer
 * a strictly deterministic data structure (for debugging, say), you can
 * <TT>\#define SoX_SWEEP_AB_TREE</TT> to get one.
 *
 * For debuggung purposes, a partial consistency check of the internal
 * sweep data structures can be enabled with
 * <TT>\#define SoX_CHECK_VALIDITY</TT>.
 */
template <class InputIterator, class CurveSweepTraits_2>
void sweep_curves(InputIterator first, InputIterator beyond, 
                  leda::GRAPH<typename CurveSweepTraits_2::Point_2, 
                  typename std::iterator_traits<InputIterator>::value_type>& G,
                  bool use_optimization, const CurveSweepTraits_2& traits,
                  bool allow_multi_edges = false) {
    
    typedef CurveSweepTraits_2                          Curve_sweep_traits_2;

    typedef Curve_sweep_traits_2                        CST;
    
    typedef typename CST::Point_2                       Point_2;
    typedef typename CST::Segment_2                     Segment_2;
    
    typedef LEDA_compare< typename CST::Compare_xy_2, Point_2 > Compare_xy;
    Compare_xy compare_xy(
            traits.compare_xy_2_object()
    );
    typename CST::Less_xy_2 less_xy(
            traits.less_xy_2_object()
    );
    typename CST::Construct_segment_2 construct_segment(
            traits.construct_segment_2_object()
    );
    typename CST::Source_2 source(
            traits.source_2_object()
    );
    typename CST::Target_2 target(
            traits.target_2_object()
    );
    typename CST::New_endpoints_2 new_endpoints(
            traits.new_endpoints_2_object()
    );
    typename CST::New_endpoints_opposite_2 new_endpoints_opposite(
            traits.new_endpoints_opposite_2_object()
    );
    typename CST::Multiplicity_of_intersection_2 multiplicity_of_intersection(
            traits.multiplicity_of_intersection_2_object()
    );
    typename CST::Do_overlap_2 do_overlap(
            traits.do_overlap_2_object()
    );
    
    CGAL_precondition(G.empty());
    
    if (first == beyond) {
        return;  
    }
    
    leda_list<Segment_2>           internal;
    
    leda_map<Segment_2, Segment_2> original; 
    
    Point_2 p_sweep;
    // these are used solely as sentinels; they have no geometric meaning
    // but must have a value in order to have memory allocated and thus an
    // id.  
    Segment_2 upper_sentinel(construct_segment(p_sweep));
    Segment_2 lower_sentinel(construct_segment(p_sweep));
    
    Intern::sweep_cmp_curve<CST> cmp(p_sweep, lower_sentinel, upper_sentinel,
                                     traits);
    
    leda::sortseq<Point_2, Seq_item, Sortseq_impl>    X_structure(compare_xy);
    
    leda::sortseq<Segment_2, Seq_item, Sortseq_impl>  Y_structure(cmp);
    
    leda_map<Segment_2, leda_edge>                    edge_of(nil);
    
    leda_p_queue<Point_2, Segment_2>                  curve_queue(compare_xy);
    
    leda_map2<Segment_2, Segment_2, Intern::Intersect_info<Point_2> > 
        inter_dic;
    
    SoX_TRACE(
            std::cout << "initializing X structure ... ";
            std::cout.flush();  
    );
    Segment_2 s;
    InputIterator it;
    for (it = first; it != beyond; it++) {
        s = *it;
        Seq_item it1 = X_structure.insert(source(s), Seq_item(nil));
        Seq_item it2 = X_structure.insert(target(s), Seq_item(nil));
        
        if (it1 == it2) {
            continue;  // ignore zero-length segments
        }
        
        Point_2 p = X_structure.key(it1);
        Point_2 q = X_structure.key(it2);
        
        Segment_2 s1; 
        
        if (less_xy(p, q)) {
            s1 = new_endpoints(s,p,q);
        } else {
            s1 = new_endpoints_opposite(s,p,q);
        }
        
        original[s1] = s; 
        internal.append(s1);
        curve_queue.insert(source(s1),s1);
    }
    SoX_TRACE( std::cout << "done." << std::endl; );
    
    
    Y_structure.insert(upper_sentinel, Seq_item(nil));
    Y_structure.insert(lower_sentinel, Seq_item(nil));
    
    Segment_2 next_curve;
    
    // curve_queue is empty if there are only isolated points
    if (curve_queue.empty()) {
        next_curve = upper_sentinel;
    } else {
        next_curve = curve_queue.inf(curve_queue.find_min());
    }
    
    while (!X_structure.empty()) {
        Seq_item event = X_structure.min();
        
        SoX_TRACE(draw_point(W, p_sweep, CGAL::WHITE, true););
        
        p_sweep = X_structure.key(event);
        
        SoX_SWEEP_LINE(
                draw_sweep_line(*W, p_sweep, leda::color("blue"));
                draw_segments(*W, first, beyond);
                draw_graph_nodes(*W);
                W->read_mouse();
                draw_sweep_line(*W, p_sweep, leda::color("white"));
        );
        
        SoX_TRACE(draw_point(W, p_sweep, CGAL::VIOLET, true););
        
        cmp.set_position(p_sweep);
        
        leda_node v = G.new_node(p_sweep);
        
        Seq_item sit = X_structure.inf(event);
        
        if (sit == nil) {
            sit = Y_structure.lookup(construct_segment(p_sweep));
        }
        
        Seq_item sit_succ  =     nil;
        Seq_item sit_pred  =     nil;
        Seq_item sit_pred_succ = nil;
        Seq_item sit_first =     nil;
        
        if (sit != nil) { 
            // walk up
            std::cout << "1" << std::endl;
            Y_structure.inf(sit);
            std::cout << "2" << std::endl;
            std::cout << "equal: " << (sit == Y_structure.max_item()) << std::endl;
            Y_structure.succ(sit);
            std::cout << "3" << std::endl;
            while ( Y_structure.inf(sit) == event || 
                    Y_structure.inf(sit) == Y_structure.succ(sit) ) {
                SoX_TRACE(draw_segment(W, Y_structure.key(sit), 
                                       CGAL::ORANGE););
                SoX_TRACE(draw_segment(W, Y_structure.key(sit), 
                                       CGAL::BLUE););
                sit = Y_structure.succ(sit);
            }
            
            // sit_succ is the curve above those that correspond to this event
            sit_succ = Y_structure.succ(sit);
            SoX_TRACE(draw_segment(W, Y_structure.key(sit_succ), 
                                   CGAL::BLACK););
            
            // walk down 
            bool overlapping;
            do {
                s = Y_structure.key(sit);
                SoX_TRACE(draw_segment(W, s, CGAL::RED););
                Seq_item it = Y_structure.pred(sit);
                overlapping = Y_structure.inf(it) == sit;
                leda_edge e = edge_of[Y_structure.key(sit)];
                CGAL_assertion(e != nil);
                if (allow_multi_edges || !overlapping)
                    Intern::link_as_target_and_append(G, v, e, s);
                
                if (ID_Number(p_sweep) == ID_Number(target(s))) {  
                    // ending segment
                    if (overlapping) {
                        Y_structure.change_inf(it, Y_structure.inf(sit));
                    }
                    Y_structure.del_item(sit);
                } else { // passing segment
                    if (Y_structure.inf(sit) != Y_structure.succ(sit)) {
                        Y_structure.change_inf(sit, Seq_item(nil));
                    }
                }
                SoX_TRACE(draw_segment(W, s, CGAL::BLUE););
                sit = it; 
            } while ( Y_structure.inf(sit) == event || overlapping ||
                      Y_structure.inf(sit) == Y_structure.succ(sit) );
            
            // sit_pred is the curve below those that correspond to this event
            sit_pred = sit;
            SoX_TRACE(draw_segment(W, Y_structure.key(sit_pred), 
                                   CGAL::BLACK););
            sit_first = Y_structure.succ(sit_pred);
            sit_pred_succ = sit_first;  
            
            ///////////////////////////////////////////////////////////////////
            // reverse subsequences of overlapping curves and curves 
            // that have the sweep point as a root of multiplicity > m
            //  overlapping segments should get reversed an odd number of times
            
            sit = Y_structure.succ(sit_pred);
#ifdef SoX_PRINT_VERTICES
            bool none_passing = sit == sit_succ;
#endif // SoX_PRINT_VERTICES
            
            Seq_item sit_next;
            Segment_2 s0;
            Segment_2 s1;
            
            std::vector<int> mults;


#if SoX_LINEAR_REORDER // new reorder
            std::cout << "before reverse" << std::endl;
            print_y_structure(Y_structure);

            std::vector< Segment_2 > items;
            typedef 
                std::map< Segment_2, Seq_item, LiS::Id_less_than< Segment_2 > >
                Ovl_map;
            Ovl_map ovl;
            
            bool overlap_found = false;

            // copy to list that must be reordered
            for (; sit != sit_succ; sit = Y_structure.succ(sit)) {
                sit_next = Y_structure.succ(sit);
                s0 = Y_structure.key(sit);
                s1 = Y_structure.key(sit_next);
                items.push_back(s0);
                if (Y_structure.inf(sit) == sit_next) {
                    // the two segments overlap
                    mults.push_back(1024);
                    ovl[s0] = sit_succ;
                    overlap_found = true;
                } else if (sit_next != sit_succ) {
                    mults.push_back(
                            multiplicity_of_intersection(s0,s1,p_sweep)
                    );
                } 
            }
            
            if (static_cast< int >(items.size() > 1)) {
                
                // TODO overlapping curves
                // and build tree to compute new sequence
                typedef SoX::Reorder_tree< Segment_2 > Reorder_tree;
                Reorder_tree tree(
                        items.begin(), items.end(), mults.begin(), mults.end()
                );
                
                // delete the sequence
                leda::sortseq<Segment_2, Seq_item, Sortseq_impl> seq_pass(cmp);
#if 1
                // copy to list that must be reordered
                for (sit = Y_structure.succ(sit_pred); 
                     sit != sit_succ; sit = Y_structure.succ(sit)) {
                    Y_structure.del_item(sit);
                }
#if 0
                Y_structure.delete_subsequence(Y_structure.succ(sit_pred), 
                                               sit_succ, 
                                               seq_pass);
#endif
#else
                Y_structure.delete_subsequence(Y_structure.succ(sit_pred), 
                                               Y_structure.pred(sit_succ),
                                               seq_pass);
#endif

                // readd ending with tom sit_succ
                Seq_item sit_insert = sit_succ;
                
                // traverse tree
                for (typename Reorder_tree::reverse_iterator trit = 
                         tree.rbegin();
                     trit != tree.rend(); trit++) {
                    //Segment_2 seg = seq_pass.key(*trit);
                    typename Ovl_map::iterator ovlit =
                        ovl.find(*trit);
#if 0
                    Y_structure.insert(*trit, nil);
#else
                    //sit_insert =
                        Y_structure.insert_at(
                                sit_pred, *trit, 
                                nil,
                                leda::after
                        );
#endif
                }
                
                Seq_item tmp, tmp_next;
                if (overlap_found) {
                    for (tmp = Y_structure.succ(sit_pred);
                         tmp != sit_succ; 
                         tmp = Y_structure.succ(tmp)) {
                        tmp_next = Y_structure.succ(tmp);
                        if (tmp_next == Y_structure.max_item()) {
                            break;
                        }
                        if (tmp_next != sit_succ) {
                            s0 = Y_structure.key(tmp);
                            s1 = Y_structure.key(tmp_next);
                            if (do_overlap(s0, s1)) {
                                Y_structure.change_inf(tmp, tmp_next);
                            }
                        }
                    }
                }
            }
            
            std::cout << "after reverse" << std::endl;
            print_y_structure(Y_structure);

            
#else // old reorder
            // determine the multiplicity of intersections between 
            // neighboring pairs of segments in the sequence 
            // [sit_first, sit_succ) using -1 to indicate
            // that two segments overlap
            int max_mult = -1;
            while ( sit != sit_succ ) { 
                sit_next = Y_structure.succ(sit);
                s0 = Y_structure.key(sit);
                s1 = Y_structure.key(sit_next);
                if (Y_structure.inf(sit) == sit_next) {
                    // the two segments overlap
                    mults.push_back(-1);
                } else if (sit_next != sit_succ) {
                    mults.push_back(
                            multiplicity_of_intersection(s0,s1,p_sweep)
                    );
                    if (mults[mults.size()-1] > max_mult)
                        max_mult = mults[mults.size()-1];
                }
                sit = sit_next;
            }
#ifdef SoX_PRINT_VERTICES
            std::cout << "VERTEX maxmult=" << max_mult
                      << ", #passing_segs=" 
                      << (none_passing ? 0 : mults.size()+1)
                      << std::endl;
#endif // SoX_PRINT_VERTICES
            
            // for intersection with sit_succ 
            // (which doesn't go through p_sweep)
            mults.push_back(-2);
            
            Seq_item sub_first;
            Seq_item sub_last;
            Seq_item sub_next;
            std::vector<int>::iterator mult_it, mult_first;
            
            // max_mult should always be an even number to make sure the 
            // overlapping curves are reversed an odd number of times 
            // (and thus remain in the same order after the reversal of 
            // the entire sequence).
            if (max_mult > 0)
                max_mult = max_mult % 2 == 0 ? max_mult: max_mult+1;
            
            do {
                sit = Y_structure.succ(sit_pred);
                mult_it = mults.begin();
                while ( sit != sit_succ ) { 
                    sub_first = sit;
                    sub_last  = sub_first;
                    sub_next  = Y_structure.succ(sub_last);           
                    
                    s0 = Y_structure.key(sub_last);
                    s1 = Y_structure.key(sub_next);
                    
                    mult_first = mult_it;
                    while (*mult_it == -1 || 
                           (sub_next != sit_succ && *mult_it >= max_mult)) {
                        sub_last = sub_next;
                        sub_next = Y_structure.succ(sub_last);           
                        s0 = Y_structure.key(sub_last);
                        s1 = Y_structure.key(sub_next);
                        mult_it++;
                    }
                    if (sub_last != sub_first) {
                        Y_structure.reverse_items(sub_first, sub_last);
                        std::reverse(mult_first, mult_it);
                    }
                    
                    sit = Y_structure.succ(sub_first);
                    mult_it++;
                }
                max_mult--;
            } while (max_mult > 1);
            
            // reverse the entire bundle
            
            if ( sit_first != sit_succ ) {
                Y_structure.reverse_items(Y_structure.succ(sit_pred),
                                          Y_structure.pred(sit_succ));  
            }
#endif
        }
        
        // insert curves starting at p_sweep
        while (!curve_queue.empty() && 
               ID_Number(p_sweep) == ID_Number(source(next_curve))) { 
            SoX_TRACE(draw_segment(W, next_curve, CGAL::RED););
            SoX_TRACE(draw_segment(W, next_curve, CGAL::BLUE););
            Seq_item s_sit = Y_structure.locate(next_curve);
            Seq_item p_sit = Y_structure.pred(s_sit);
            s = Y_structure.key(s_sit);
            if (ID_Number(s) != ID_Number(upper_sentinel) && 
                do_overlap(s,next_curve)) {
                sit = Y_structure.insert_at(s_sit, next_curve, s_sit);
            } else {
                sit = Y_structure.insert_at(s_sit, next_curve, Seq_item(nil));
            }
            
            s = Y_structure.key(p_sit);
            
            if (ID_Number(s) != ID_Number(lower_sentinel) && 
                do_overlap(s,next_curve)) {
                Y_structure.change_inf(p_sit, sit);
            }
            
            if (ID_Number(Y_structure.key(s_sit)) != 
                ID_Number(upper_sentinel)) {
                X_structure.insert(target(next_curve), sit);
            }        
            
            if (sit_succ == nil) { 
                sit_succ = s_sit; 
                sit_pred = p_sit;  
                sit_pred_succ = sit_succ;  
            }
            // delete minimum and assign new minimum to next_curve
            
            curve_queue.del_min();
            if (!curve_queue.empty()) {
                next_curve= curve_queue.inf(curve_queue.find_min());
            }
        }
        
        // insert new edge stubs in the graph not yet linked to their targets
        if (sit_succ != nil) {
            for(Seq_item sitl = Y_structure.pred(sit_succ); 
                sitl != sit_pred;
                sitl = Y_structure.pred(sitl)) {
                if (allow_multi_edges) {
                    edge_of[Y_structure.key(sitl)] = 
                        Intern::new_halfedge_pair_at_source(G, v);
                } else {
                    if (Y_structure.inf(sitl) == Y_structure.succ(sitl)) {
                        // overlapping
                        edge_of[Y_structure.key(sitl)] = 
                            edge_of[ Y_structure.key(Y_structure.succ(sitl))];
                    } else  {
                        edge_of[Y_structure.key(sitl)] = 
                            Intern::new_halfedge_pair_at_source(G, v);
                    }
                }
            }
        }
        
        // must find intersections of all adjacent curves passing 
        // through p_sweep
        if (sit_pred != nil) { 
            if (!use_optimization) { 
                sit = sit_pred;
                while (sit != sit_succ) {
                    if (Y_structure.succ(sit) != Y_structure.inf(sit)) {
                        // not overlapping
                        Y_structure.change_inf(sit,Seq_item(nil));
                        Intern::compute_intersection(X_structure, Y_structure, 
                                                     p_sweep, 
                                                     sit, traits);
                    }
                    sit = Y_structure.succ(sit);
                }
            } else { 
                sit = sit_pred;
                while (sit != sit_succ) {
                    if (Y_structure.succ(sit) != Y_structure.inf(sit)) {
                        // not overlapping
                        Y_structure.change_inf(sit,Seq_item(nil));
                        Intern::compute_intersection(X_structure, Y_structure, 
                                                     p_sweep, 
                                                     inter_dic, sit, traits);
                    }
                    sit = Y_structure.succ(sit);
                }
            }
        }
        
        X_structure.del_item(event);
        SoX_VALIDITY_CHECK(Intern::structures_are_valid(Y_structure, 
                                                        X_structure, 
                                                        p_sweep, traits));
    }
    leda_edge e;
    forall_edges(e, G) G[e] = original[G[e]];
}

} // namespace SoX

#endif // SoX_SWEEP_CURVES_H

// EOF
