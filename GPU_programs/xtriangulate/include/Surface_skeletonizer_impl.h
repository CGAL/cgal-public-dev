// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://asm@scm.gforge.inria.fr/svn/cgal/branches/unsorted-branches/eric/Numerical_algebraic_kernel_d/include/CGAL/Algebraic_kernel_d/Xy_coordinate_2.h $
// $Id: Xy_coordinate_2.h 67816 2012-02-20 13:51:29Z asm $
//
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_SURFACE_SKELETONIZER_IMPL_H
#define CGAL_SURFACE_SKELETONIZER_IMPL_H

#include <CGAL/basic.h>
#include <CGAL/Timer.h>
// our lovely root isolator
// #include <CGAL/Algebraic_curve_kernel_2/IARD_root_isolator_2.h>

#define SKELETONIZER_THRESHOLD      1e-6
#define SKELETONIZER_THRESHOLD_LOG \
    (int)std::floor(-log2(SKELETONIZER_THRESHOLD))

#define ROUND_BFI_PREC 20

#ifndef XTri_OUT
#define XTri_USE_OUT
#endif // XTri_OUT

#ifdef XTri_USE_OUT

#define XTri_OUT(x) std::cerr << x
#define dbl(x) CGAL::to_double(x)
#define TM_START(timer) timer.start()
#define TM_STOP(timer) timer.stop()
#else

#define XTri_OUT(x) static_cast<void>(0)
#define dbl(x) ???
#define TM_START(timer) static_cast<void>(0)
#define TM_STOP(timer) static_cast<void>(0)
#endif // XTri_USE_OUT

#ifndef STILL_ALIVE
#define STILL_ALIVE XTri_OUT(__LINE__ << "\n");
#endif

#ifdef CGAL_XTri_USE_TIMERS
#warning timers!
extern CGAL::Timer tm_external_timer;
#endif

namespace CGAL {

namespace {

// type of vertical stack line
enum { INTERM_LINE = 0, CLIP_LINE, EVENT_LINE };

// characterizes horizontal clipping
enum { NO_CLIP = 0, BOTTOM_CLIP = 1, TOP_CLIP = 2, Z_BOTTOM_CLIP = 4,
    Z_TOP_CLIP = 8 };

// point special flag
enum { INTERIOR_PT = 0, ON_SIL = 1, ON_Z_CURVE = 2 };

struct Sil_point_2 {
    unsigned type; // point special type (silhouette or z-curve)
    // starting index of this point lifts in points array
    unsigned istart;
    unsigned n_zs; // total # of valid z-lifts
    unsigned z_skips_below; // # of lifts skipped below lower z-curve
    unsigned z_skips_above; // # of lifts skipped above upper z-curve
    unsigned idcel; // index in dcel feature map
};

struct Feature_desc { // arrangement feature descriptor
    Feature_desc() { }

    Feature_desc(unsigned _idcel) :
        idcel(_idcel) { }

    unsigned idcel; // dcel index of this feature (face/edge)
};

struct Point_index {

    Point_index() { }
    Point_index(unsigned _x, unsigned _y) : x(_x), y(_y) {  }

    friend bool operator <(const Point_index& p1, const Point_index& p2) {
        return (p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y));
    }

    friend std::ostream& operator <<(std::ostream& os, const Point_index& p) {
        os << "[" << p.x << "; " << p.y << "]";
        return os;
    }

    unsigned x, y;
};

} // anonymous namespace

// NOTE: we can also use SurfaceAnalysisInterface< Kernel_2 >
// to ensure that SFA_interface is consistent with Skeletonizer number types
template < class AlgebraicKernel_2, class SurfaceAnalysisInterface >
class Surface_skeletonizer_impl {

public:

    typedef SurfaceAnalysisInterface Surface_analysis_interface;
    
    //! kernels
    typedef AlgebraicKernel_2 Kernel_2;
    typedef typename Kernel_2::Algebraic_kernel_d_1 Kernel_1;

    //! number types
    typedef typename Kernel_2::Arithmetic_kernel AT;
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;
    typedef typename AT::Bigfloat Bigfloat;
    typedef typename AT::Bigfloat_interval BFI;

    //! polynomials
    typedef CGAL::Polynomial< Integer > Polynomial_1;
    typedef CGAL::Polynomial< Polynomial_1 > Polynomial_2;
    typedef CGAL::Polynomial< Polynomial_2 > Polynomial_3;

    //! analyses
    typedef typename Kernel_2::Curve_analysis_2 CA_2;
    typedef typename Kernel_2::Curve_pair_analysis_2 CPA_2;

    //! status lines 
    typedef typename CA_2::Status_line_1 CA_line_1;
    typedef typename CPA_2::Status_line_1 CPA_line_1;

    //! coordinates
    typedef typename Kernel_2::Algebraic_real_1 X_coordinate_1;
    typedef typename Kernel_2::Algebraic_real_2 Xy_coordinate_2;

    //! vectors
    typedef std::vector< double > Double_vector;
    typedef std::vector< X_coordinate_1 > X_coord_vector;

private:
    //! internal stuff

    typedef CGAL::Adjacencies_3::Adjacency_vector Adjacency_vector;
    typedef CGAL::Adjacencies_3::Adjacency_pair Adj_pair;

    //! container of points along one y-vertical fiber (stack)
    typedef std::vector< Sil_point_2 > Point_stack;
    //! global container of points for all fibers (skeleton data)
    typedef std::vector< Point_stack > Skeleton_data;

    //! stuff for adjacency information
    typedef std::pair< Point_index, Point_index > Point_pair_index;
    typedef std::map< Point_pair_index, unsigned > Adjacency_map;
    typedef std::vector< Adjacency_vector > Adjacency_data;

    //! map from a pair of features (vertices/edges/faces) to an index in
    //! \c Adjacency_data
    typedef std::pair< const void *, const void * > Addr_pair;
    typedef std::map< Addr_pair, unsigned > Feature_adjacency_map;

    //! vector of indices
    typedef std::vector< unsigned > Index_vector;
    //! bivariate polynomial 
    typedef CGAL::Polynomial_type_generator< double, 2 >::Type
        Poly_double_2;

    typename Kernel_1::Algebraic_real_traits::Lower_bound lbound_x;
    typename Kernel_1::Algebraic_real_traits::Upper_bound ubound_x;
    // should be used instead:
    // typename AK_1::Approximate_relative_1 approximate_x;
    
public:
    //! constructors

Surface_skeletonizer_impl(Kernel_2 *_kernel_2 = 0) :
        kernel_2(0), is_valid(false) { }

public:
//! interface:

void reset_parameters(unsigned _n_slices_x = 1, unsigned _n_slices_y = 1,
            double _en_left = 0.1, double _en_right = 0.1,
            double _en_btm = 0.1, double _en_top = 0.1,
            double _z_below = -2.0, double _z_above = 2.0) {

    is_valid = false;
    n_slices_x = _n_slices_x, n_slices_y = _n_slices_y;
    enlarge_left = _en_left, enlarge_right = _en_right;
    enlarge_btm = _en_btm, enlarge_top = _en_top;
    z_below = _z_below, z_above = _z_above;

    g_start_idx = 0;

    int prec = ROUND_BFI_PREC;
    round_bfi(enlarge_left, prec, 0);
    round_bfi(enlarge_right, prec, 1);
    round_bfi(enlarge_btm, prec, 0);
    round_bfi(enlarge_top, prec, 1);
    round_bfi(z_below, prec, 0); // round down
    round_bfi(z_above, prec, 1); // round up

    std::cout << "sx = " << n_slices_x << "\n";
    std::cout << "sy = " << n_slices_y << "\n";
    std::cout << "enlarge_left = " << _en_left << "\n";
    std::cout << "enlarge_right = " << _en_right << "\n";
    std::cout << "enlarge_btm = " << _en_btm << "\n";
    std::cout << "enlarge_top = " << _en_top << "\n";
    std::cout << "z_below = " << z_below << "\n";
    std::cout << "z_above = " << z_above << "\n";
}

bool triangulation_valid() {
    return (is_valid && surf_engine->surface_valid());
}

void set_surface_engine(Surface_analysis_interface *s) {
    surf_engine = s;
    kernel_2 = &surf_engine->kernel_2_inst();
}

void skeletonize(bool use_auto_bounds = true) {

    clear();

#ifdef XTri_USE_OUT
    tm_lift_on_sil.reset();
    tm_lift_in_face.reset();
    tm_arcavoid_in_face.reset();
    tm_create_skeletons.reset();
    tm_tri_silhouette.reset();
    tm_y_range.reset();

//     tm_external_timer.reset();
#endif

    CGAL::set_pretty_mode(std::cerr);
    CGAL::set_pretty_mode(std::cout);

    TM_START(tm_create_skeletons);
    create_skeletons(use_auto_bounds);
    TM_STOP(tm_create_skeletons);

    TM_START(tm_tri_silhouette);
    triangulate_silhouette();
    TM_STOP(tm_tri_silhouette);

#ifdef XTri_USE_OUT
    std::cout << "\ntm_load_surface: " << tm_load_surface.time() <<
        "\ntm_compute_y_range: " << tm_y_range.time() <<
        "\ntm_lift_on_sil: " << tm_lift_on_sil.time() <<
        "\ntm_lift_in_face: " << tm_lift_in_face.time() <<
        "\ntm_arcavoid_in_face: " << tm_arcavoid_in_face.time() <<
        "\ntm_create_skeletons: " << tm_create_skeletons.time() <<
        "\ntm_tri_silhouette: " << tm_tri_silhouette.time() << "\n";
        
//         "\ntm_external_timer: " << tm_external_timer.time() << "\n";
#endif
        
    is_valid = true;
}

bool load_surface(const char *filename) {

    TM_START(tm_load_surface);
    // note that load surface does not modify is_valid flag
    bool res = surf_engine->load_surface(filename);
    TM_STOP(tm_load_surface);
    return res;
}

bool overlay_surface(const char *filename) {
//     return surf_engine->overlay_surface(filename);
    return false; // NYI
}

const std::vector< unsigned >& get_zlift_ranges() const {
    return zlift_ranges;
}

const Triangles_vector& get_triangles() const {
    return triangles;
}

const Triangles_vector& get_sil_triangles() const {
    return sil_triangles;
}

const Point_vec_3d& get_vertices() const {
    return points;
}

double get_z_silhouette() const {
    return CGAL::to_double(z_silhouette);
}

void get_surface_bounds(double& left, double& right,
                        double& btm, double& top) const {

    if(!is_valid)
        return;
    left = CGAL::to_double(g_left), right = CGAL::to_double(g_right);
    btm = CGAL::to_double(g_btm), top = CGAL::to_double(g_top);
}

void compute_normals(const Point_vec_3d& in,
        Point_vec_3d& normals, Point_3d& centroid, bool is_inversed) const {

    surf_engine->compute_normals(in, normals, is_inversed);
    fix_normals(in, triangles, normals, centroid);
}

void fix_normals(const Point_vec_3d& in,
        const Triangles_vector& tris, Point_vec_3d& normals,
                    Point_3d& centroid) const {

    Point_vec_3d tri_norms(in.size());
    std::map< int, int > vmap;

    typedef std::map< int, int >::iterator iterator;
    std::pair< iterator, bool > res;

    double ux, uy, uz, vx, vy, vz;
    Point_3d n;

    unsigned n_centroid = 0;
    centroid[0] = 0; centroid[1] = 0; centroid[2] = 0;
    for(unsigned i = 0; i < tris.size(); i++) {

        const Tri_index& t = tris[i];
        const Point_3d& p1 = in[t[0]],
            p2 = in[t[1]], p3 = in[t[2]];
        
        ux = p2[0] - p1[0], uy = p2[1] - p1[1], uz = p2[2] - p1[2];
        vx = p3[0] - p1[0], vy = p3[1] - p1[1], vz = p3[2] - p1[2];

        n[0] = uy*vz - uz*vy, n[1] = uz*vx - ux*vz, n[2] = ux*vy - uy*vx;
        double det = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];

        det = (det >= 1e-13 ? 1.0 / std::sqrt(det) : 0.0);
        n[0] *= det, n[1] *= det, n[2] *= det;

        for(int j = 0; j < 3; j++) {
            const Point_3d& pv = in[t[j]];
            centroid[0] += pv[0];
            centroid[1] += pv[1];
            centroid[2] += pv[2];
            n_centroid++;

            Point_3d& p = tri_norms[t[j]];
            res = vmap.insert(std::make_pair(t[j], 1));
            if(res.second) { // inserted
                p[0] = n[0], p[1] = n[1], p[2] = n[2];
            } else {
                int nv = res.first->second++; // # of adjacent faces
                det = 1.0 / (double)(nv + 1);
                p[0] = (p[0] * nv + n[0]) * det;
                p[1] = (p[1] * nv + n[1]) * det;
                p[2] = (p[2] * nv + n[2]) * det;
            }
        }
    }

    if(n_centroid != 0) {
        centroid[0] /= n_centroid;
        centroid[1] /= n_centroid;
        centroid[2] /= n_centroid;
    }

    for(unsigned i = 0; i < in.size(); i++) {
        Point_3d& n = normals[i], tri_n = tri_norms[i];

        double det = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        if(det < 1e-10) {
            n = tri_n; // replace degenerate normal with triangle's normal
            continue;
        }

        det = n[0] * tri_n[0] + n[1] * tri_n[1] + n[2] * tri_n[2];

        // NOTE NOTE NOTE: investigate wrong_triangle_colors bug..
        // NOTE: if angle is close to 90 degrees => makes sense
        // to replace by triangle's normal
        if(std::abs(det) < 0.05) {
            n = tri_n; // replace degenerate normal with triangle's normal
            continue;
        }
        if(det < 0) {
            n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
        }
    }

    // once you've computed normals, do the following:
    // if original normal is 0 or points wrong direction, reverse it!!
}

~Surface_skeletonizer_impl() {
    clear();
}

private:
//! private stuff    
        
void _refine_x(const X_coordinate_1& x, const Rational& bound) {
    while(ubound_x(x) - lbound_x(x) >= bound)
        x.refine();
}

void _refine_xy(const Xy_coordinate_2& xy, const Rational& bound) {

    while(xy.upper_bound_y() - xy.lower_bound_y() >= bound)
        xy.refine_y();
}

//! rounds a rational \c x to \c prec bits up or down
void round_bfi(Rational& x, int prec, int dir = 0) {

    prec = std::max(prec, 2);
    int old_prec = CGAL::set_precision(BFI(), prec);
    BFI xf = CGAL::convert_to_bfi(x);
    CGAL::set_precision(BFI(), old_prec);
    x = (dir == 0 ? CGAL::lower(xf) : CGAL::upper(xf));
}

//! returns \c 1 if features \c i1 and \c i2 share the same face,
//! \c -1 if one of the features does not belong to the face
unsigned share_same_face(unsigned i1, unsigned i2) {

    std::cerr << "i2: " << i2 << "; dcel_vec: " << dcel_vec.size() << "\n";
    CGAL_assertion(i1 < dcel_vec.size());
    CGAL_assertion(i2 < dcel_vec.size());

    if(dcel_vec[i1].first != CGAL::FACE ||
            dcel_vec[i2].first != CGAL::FACE)
        return -1; // wrong feature

    if(dcel_vec[i1].second == dcel_vec[i2].second)
        return 1;
    return 0;
}

//! tests whether points [i1; 1] and [i2; 0] lie in the same face
//! or on one side of z-curve w.r.t. the point [i1; 0] which must be either
//! silhouette or z-point
bool curve_on_one_side(unsigned i1, unsigned i2) {
//NOTE: if a point belongs simultaneously to z-curve and to silhouette,
// check *z-curve* first

    Poly_double_2 ppoly = g_z_poly;
    if(sk_data[i1][0].type == ON_SIL) { // on sil but *not* on z-curve

        unsigned res = share_same_face(sk_data[i1][1].idcel,
                    sk_data[i2][0].idcel);
        if(res != -1u) // this means that one of the points does not lie
            return res;  // in a face

        XTri_OUT("curve_on_one_side: checking (" << i1 << "; 1) and (" <<
                i2 << "; 0) degenerate case\n");
        ppoly = g_sil_poly;
    } else
        CGAL_assertion(sk_data[i1][0].type & ON_Z_CURVE);

    const Point_3d& p1 = points[sk_data[i1][1].istart],
                p2 = points[sk_data[i2][0].istart],
                p1_low = points[sk_data[i1][0].istart];

    CGAL::Polynomial_traits_d< Poly_double_2 >::Sign_at sign_at;

    // we choose a middle point to avoid the situation where
    // i2 point lies on the curve
    double s1[] = {p1[0], p1[1]}, s2[] = {(p1_low[0] + p2[0]) / 2,
                (p1_low[1] + p2[1]) / 2};
    CGAL::Sign x1 = sign_at(ppoly, s1, s1 + 2),
          x2 = sign_at(ppoly, s2, s2 + 2);

//     XTri_OUT("checking curve side: " << ppoly << "; (" <<
//             p1[0] << ", " << p1[1] << ") and (" << p2[0] << ", " << p2[1]
//             << "); x1 = " << x1 << "; x2 = " << x2 << "\n");

    return (x1 * x2 > 0);
}

//! computes range of y-coordinates for a surface using event points of
//! silhouette and z-curves
void compute_y_range(const CPA_2& cpa,
    double& y_min, double& y_max, bool& is_set) {

    is_set = false;
    surf_engine->compute_y_range(y_min, y_max, is_set);

    CA_2 ca1 = cpa.curve_analysis(0),
        ca2 = cpa.curve_analysis(1);

//     std::cerr << ca1.polynomial_2() << "; " << ca2.polynomial_2() << "\n";

    for(int i = 0; i < cpa.number_of_status_lines_with_event(); i++) {

        XTri_OUT("event: " << i << "\n");

        // NOTE: there is no need to take event points exactly:
        // can use status_line_of_interval instead..
        CPA_line_1 cpaline = cpa.status_line_at_event(i);

        if(!cpaline.is_intersection())
            continue;

        X_coordinate_1 x = cpaline.x();
        int m = cpaline.number_of_events();
        for(int j = 0; j < m; j++) {

            typename CPA_line_1::Arc_pair arcs = cpaline.curves_at_event(j);
            if(arcs.first == -1 || arcs.second == -1)
                continue;

            CA_line_1 line = ca1.status_line_at_exact_x(x);
            const std::pair< double, double>& xy =
                line.xy_coordinate_2(arcs.first).to_double();

            XTri_OUT("i1: " << i << "; j: " << j << "; xy: " <<
                xy.first << "; " << xy.second << "\n");
            if(!is_set || xy.second < y_min)
                y_min = xy.second;
            if(!is_set || xy.second > y_max) {
                y_max = xy.second;
                is_set = true;
            }
//             if(m == 1)  // only one event: done
//                 break;
//             arcs = cpaline.curves_at_event(m-1);
        }
    }

    // if range collapses, use intermediate lines additionally
    if(is_set && fabs(y_min - y_max) > 1e-5)
        return;

    // NOTE: look at how the topology graph is computed..
    for(int i = 0; i < cpa.number_of_status_lines_with_event(); i++) {

        // NOTE: there is no need to take event points exactly:
        // can use status_line_of_interval instead..
        CPA_line_1 cpaline = cpa.status_line_of_interval(i);

//         if(!cpaline.is_intersection())
//             continue;

        X_coordinate_1 x = cpaline.x();
        int m = cpaline.number_of_events();
        for(int j = 0; j < m; j++) {

            typename CPA_line_1::Arc_pair arcs = cpaline.curves_at_event(j);
//             if(arcs.first == -1 || arcs.second == -1)
//                 continue;
            CA_line_1 line;
            int arc = arcs.first;
            if(arcs.first != -1) {
                line = ca1.status_line_at_exact_x(x);
            } else {
                line = ca2.status_line_at_exact_x(x);
                arc = arcs.second;
            }
            const std::pair< double, double>& xy =
                line.xy_coordinate_2(arc).to_double();

            XTri_OUT("i2: " << i << "; j: " << j << "; xy: " <<
                xy.first << "; " << xy.second << "\n");
            if(!is_set || xy.second < y_min)
                y_min = xy.second;
            if(!is_set || xy.second > y_max) {
                y_max = xy.second;
                is_set = true;
            }
//             if(m == 1)  // only one event: done
//                 break;
//             arcs = cpaline.curves_at_event(m-1);
        }
    }
}

//! returns starting index
void add_z_lift(const double& x, const double& y,
        const Double_vector& zs, Sil_point_2& sil_pt) {

    sil_pt.istart = g_start_idx, sil_pt.n_zs = zs.size();
    g_start_idx += sil_pt.n_zs + 1;

    Point_3d pt;
    // first add dummy points and special z-point for silhouette plane

//     XTri_OUT("add_z_lift at (" << x << "; " << y << "); " <<
//             sil_pt.z_skips_below << ": [");

    pt[0] = x, pt[1] = y, pt[2] = 0;
//     for(int i = 0; i < sil_pt.z_skips_below; i++) {
//         points.push_back(pt);
//     }
    for(int i = 0; i < static_cast<int>(zs.size()); i++) {
        pt[2] = zs[i];
        points.push_back(pt);
//         XTri_OUT(zs[i] << " ");
    }
//     XTri_OUT("]\n");
    // add 1 special z-point to show silhouette triangulation
    pt[2] = CGAL::to_double(z_silhouette);
    points.push_back(pt);
}

//! \c break_on_hit - indicates whether to exit the loop if an event point
//! is very close to the boundary, otherwise this point will be skipped
//! return \c true if \c break_on_hit fired
//! \c last_skipped: last event index being skipped because as was
//! too close to the curve
bool first_event_above_y(const CA_line_1& sline, int& idx,
        int highest_idx, const Rational& y_bnd, bool break_on_hit) {

    bool get_next = true;
    Xy_coordinate_2 xy;
    int last_skipped = -1;

    XTri_OUT(idx << ": first_event_above_y y = " << dbl(y_bnd) << "\n");

    while(idx <= highest_idx) {
        if(get_next) {
            xy = sline.algebraic_real_2(idx);
            get_next = false;
        }

        if(xy.lower_bound_y() > y_bnd) {
            XTri_OUT("lower_bound_y: " << dbl(xy.lower_bound_y()) << "\n");
            break;
        }
            
        if(xy.upper_bound_y() < y_bnd) {
           XTri_OUT("skipping " << idx << ": [" <<
                dbl(xy.lower_bound_y()) << "; " <<
                dbl(xy.upper_bound_y()) << "]\n");
            idx++, get_next = true;
            continue;
        }

        // NOTE NOTE: be carefull here
        if(xy.upper_bound_y() - xy.lower_bound_y()
                    < SKELETONIZER_THRESHOLD / 4) {

           XTri_OUT(dbl(y_bnd) << " hit at: " << idx << ": [" <<
                dbl(xy.lower_bound_y()) << "; " <<
                dbl(xy.upper_bound_y()) << "]\n");
            
            if(break_on_hit)
                return true;
            last_skipped = idx++;
            get_next = true;
        } else
            xy.refine_y();
    }
/*
    XTri_OUT("first_event_above_y with y_bnd = " << y_bnd << " found: " <<
             idx  << "\n");*/
    return false;
}

//! attention: \c idx is decreasing, \c lowest_idx - is a lowest valid event
//! index
bool last_event_below_y(const CA_line_1& sline, int& idx,
        int lowest_idx, const Rational& y_bnd, bool break_on_hit) {

    bool get_next = true;
    Xy_coordinate_2 xy;
    int last_skipped = -1;

    XTri_OUT(idx << ": last_event_below_y y = " << dbl(y_bnd) << "\n");

    while(idx >= lowest_idx) {
        if(get_next) {
            xy = sline.algebraic_real_2(idx);
            get_next = false;
        }

        if(xy.upper_bound_y() < y_bnd) {
            XTri_OUT("upper_bound_y: " << dbl(xy.upper_bound_y()) << "\n");
            break;
        }
        if(xy.lower_bound_y() > y_bnd) {
            XTri_OUT("skipping " << idx << ": [" <<
                dbl(xy.lower_bound_y()) << "; " <<
                dbl(xy.upper_bound_y()) << "]\n");
            idx--, get_next = true;
            continue;
        }

        if(xy.upper_bound_y() - xy.lower_bound_y()
                < SKELETONIZER_THRESHOLD / 4) {
            
           XTri_OUT(dbl(y_bnd) << " hit at: " << idx << ": [" <<
                dbl(xy.lower_bound_y()) << "; " <<
                dbl(xy.upper_bound_y()) << "]\n");
            
            if(break_on_hit)
                return true;
            last_skipped = idx--;
            get_next = true;
        } else
            xy.refine_y();
    }
    return false;
}

//! \c dry_run: if true, point-location query and exact root isolation
//! is run for each face/edge encountered at the vertical line, the contents
//! of \c features get updated
//! \c z_line is event line of z-curve at x
//! \c intersection - indicates 2-curve intersection at this x (account
//! for different number of branches)
void run_vertical_line(const CA_line_1& sline,
        const CA_line_1& zline, unsigned line_idx,
     const Rational& max_step_y, Point_stack& stack, bool dry_run = false,
        int what_line = INTERM_LINE, bool intersection = false) {

    int m = sline.number_of_events(), mz = zline.number_of_events();
    stack.clear();

    X_coordinate_1 x = sline.x();
    Rational rat_x = lbound_x(x); // TODO x must be sufficiently refined !!!
    double xd = CGAL::to_double(x), yd;
    Xy_coordinate_2 xy;

    if(what_line == INTERM_LINE)
        XTri_OUT("\n************ intermediate line at: ");
    else if(what_line & EVENT_LINE)
        XTri_OUT("\n************ EVENT line at: ");
    else
        XTri_OUT("\n************ cline at: ");

    XTri_OUT(xd << "; line_idx: " << line_idx << "; # events = " << m <<
        "; # z events = " << mz << "\n\n");


    int idx_btm = 0, idx_top = m - 1, z_idx_btm = 0, z_idx_top = mz - 1;
    unsigned pt_idx = 0, clip_flag = NO_CLIP;

    bool break_on_hit = ((what_line & CLIP_LINE) != 0);
    if(first_event_above_y(sline, idx_btm, m - 1, g_btm, break_on_hit)) {
        clip_flag = BOTTOM_CLIP;
        XTri_OUT("clip_line hit at bottom boundary\n");
    }

    if(last_event_below_y(sline, idx_top, idx_btm, g_top, break_on_hit)) {
        clip_flag |= TOP_CLIP;
        XTri_OUT("clip_line hit at top boundary\n");
    }

    XTri_OUT("@@@@@@@@@@@@@@@@@@@@ computing for z-curve\n");
    //break_on_hit = ((what_line & CLIP_LINE) != 0);
    //! NOTE NOTE NOTE
    break_on_hit = true;
    //! NOTE NOTE NOTE
    if(first_event_above_y(zline, z_idx_btm, mz - 1, g_btm, break_on_hit)) {
        clip_flag = Z_BOTTOM_CLIP;
        XTri_OUT("clip_line hit z-line at bottom boundary\n");
    }
    if(last_event_below_y(zline, z_idx_top, z_idx_btm, g_top, break_on_hit)) {
        clip_flag |= Z_TOP_CLIP;
        XTri_OUT("clip_line hit z-line at top boundary\n");
    }

    int idx = idx_btm, z_idx = z_idx_btm, n_pts = 0;
    Rational cur_y = g_btm, y = cur_y, next_y;

    XTri_OUT("idx_btm = " << idx_btm << "; idx_top = " << idx_top <<
        "\n");
    XTri_OUT("z_idx_btm = " << z_idx_btm << "; z_idx_top = " << z_idx_top <<
        "\n");

    Polynomial_2 poly_at_x = surf_engine->poly_at_rational_x(rat_x);

/** *     CGAL::IARD_root_isolator_2< Integer, CGAL::Interval_nt< true > >
             puppy_isolator(poly_at_x);
*/

    std::vector< Feature_desc >::const_iterator fdi = features.begin();
    if(dry_run)
        features.clear();

    Rational beg(z_below), end(z_above);
    Rational tolerance = (end - beg) * SKELETONIZER_THRESHOLD;
    beg -= tolerance, end += tolerance;

    Index_vector sil_map_l, sil_map_r;
    // FIXME some faces can have no z-lifts: need to indentify and avoid
    // processing them
    while(y < g_top) { //!??? check loop condition

        XTri_OUT("========== idx = " << idx << " y = " << dbl(y) << "\n");

        bool do_z_point = false, z_sil_intersect = false;
        // point idx is the next point above y
        next_y = g_top; // set to the topmost point
        if(idx <= idx_top) {
            if(idx == idx_btm && (clip_flag & BOTTOM_CLIP)) {
                next_y = g_btm;
                XTri_OUT("bottom clip event");
            } else if(idx == idx_top && (clip_flag & TOP_CLIP))
                next_y = g_top;
            else {
                xy = sline.algebraic_real_2(idx);
                _refine_xy(xy, SKELETONIZER_THRESHOLD);
                next_y = xy.lower_bound_y();

                XTri_OUT("xy event\n");

                if(intersection && z_idx <= z_idx_top &&
                        xy == zline.algebraic_real_2(z_idx)) {
                    z_sil_intersect = true, z_idx++;
                    XTri_OUT("%%%%%%%%%%%%%%% z-curve hit!!!\n");
                }
            }
        }

        if(z_idx == z_idx_btm && (clip_flag & Z_BOTTOM_CLIP)) {
            next_y = g_btm;
            XTri_OUT("z-bottom clip event");
            do_z_point = true;

            /// need next z_idx which is below next_y
        } else if(z_idx <= z_idx_top) { 
         //TODO not efficient: calls refine each time the test is taken
            xy = zline.algebraic_real_2(z_idx);
            _refine_xy(xy, SKELETONIZER_THRESHOLD);
            Rational lbnd = xy.lower_bound_y();
            if(lbnd <= next_y) {
                next_y = lbnd;
                do_z_point = true;
            }
        }

        Rational stepy = (next_y - cur_y) / max_step_y;
        // round up # of max steps to compute # of uniform steps
        int nsteps = static_cast< int >(std::ceil(CGAL::to_double(stepy)));
    // maximal distance between points over vertical line is 'max_step_y'
    // with the expection that event points must be separated by at least
    // 1 intermediate point

       // HACK HACK HACK
       // problem: when y-clip occurs and sy == 1, the number of steps
       // along y is also 1

//         if(nsteps <= 1)
//             nsteps = 2;

//         //FIXME no need to insert one more point if clip line is hit
        if(nsteps <= 0)
            nsteps = 1;
        if(nsteps == 1) {
            // TODO: check if this covers z-lines, since if 0 steps was
            // run fdi pointer won't be incremeted
            if((idx > idx_btm || z_idx > z_idx_btm) &&
                    (idx <= idx_top || z_idx <= z_idx_top))
                nsteps = 2;

            else if((idx == idx_btm || z_idx == z_idx_btm) &&
                    (clip_flag & (BOTTOM_CLIP | Z_BOTTOM_CLIP)))
                nsteps = 2;
        }

        XTri_OUT("original stepy: " << dbl(stepy) <<
            "; nsteps result: " << nsteps << "\n");

        stepy = (next_y - cur_y) / nsteps;

        Rational inc(0);
        // step to 'next_y' point itself in case if no clip event encountered
        if(next_y == g_top && !(clip_flag & (TOP_CLIP | Z_TOP_CLIP)))
            inc = stepy / 2;

        XTri_OUT("cur_y: " << CGAL::to_double(cur_y) << "; next_y: " <<
            CGAL::to_double(next_y) << "; stepy: " << CGAL::to_double(stepy) << " inc = " << CGAL::to_double(inc) << "\n");

        // in case of Z_BOTTOM_CLIP we shouldn't skip g_btm point
        if(cur_y != g_btm || (clip_flag & (BOTTOM_CLIP | Z_BOTTOM_CLIP)))
            // if not on bottom boundary, advance y by 1 step
            // to skip preceding event point y
            y = cur_y + stepy;

        Sil_point_2 pt;
        pt.type = INTERIOR_PT; // nothing special abt this point

        unsigned start_n = n_pts, idcel = -1u, n_zs = -1u;
        // ensure we make at least 1 iteration
        if(y < next_y + inc) {
        
            if(dry_run) {
                TM_START(tm_lift_in_face);
                idcel = surf_engine->locate_pt_in_face(x, y, dcel_vec);
                TM_STOP(tm_lift_in_face);
                
                features.push_back(Feature_desc(idcel));
            } else {
                idcel = fdi->idcel;
                fdi++;
            }
        }

        pt.idcel = idcel;
        for(; y < next_y + inc; y += stepy, n_pts++) {

            Double_vector zs;
            //TM_START(tm_lift_in_face);
            surf_engine->approximate_zs_in_face(poly_at_x, y,
                       SKELETONIZER_THRESHOLD, z_below, z_above, zs,
                        pt.z_skips_below, pt.z_skips_above);

            yd = CGAL::to_double(y);
            add_z_lift(xd, yd, zs, pt);

            //TM_STOP(tm_lift_in_face);
            stack.push_back(pt);
        }
        XTri_OUT("; actual steps run: " << (n_pts - start_n) << "\n");

        if(idx <= idx_top || do_z_point) { // process event point
            yd = CGAL::to_double(next_y);

            typename CA_line_1::Arc_pair ipair =
                (do_z_point ? zline.number_of_incident_branches(z_idx) :
                    sline.number_of_incident_branches(idx));

            if(z_sil_intersect) { // in case of z-sil intersection
                    // we need all incident branches
                    // z_idx-1 because already incremented
                typename CA_line_1::Arc_pair zpair =
                    zline.number_of_incident_branches(z_idx - 1);
                ipair.first += zpair.first;
                ipair.second += zpair.second;
            }

            sil_map_l.insert(sil_map_l.end(), ipair.first, n_pts);
            if(what_line & EVENT_LINE) // for event lines # of branches to the
                                        // left and right differs
                sil_map_r.insert(sil_map_r.end(), ipair.second, n_pts);

            pt.type = ON_SIL;
            Double_vector zs;

            if(do_z_point)
                XTri_OUT("&&&&&&&&&&&&&&&&& z-point: ");
            else
                XTri_OUT("&&&&&&&&&&&& point on sil: ");

            pt.z_skips_below = 0, pt.z_skips_above = 0;
            if(do_z_point) {
                TM_START(tm_lift_in_face);
                pt.type |= ON_Z_CURVE;
               //! NOTE NOTE NOTE: enable point location if smth goes wrong
                // it can happen that there are no regular point before
                // this z-point, hence we need point location here
                if(idcel == -1u)
                    idcel = surf_engine->locate_pt_in_face(x, next_y, dcel_vec);
                TM_STOP(tm_lift_in_face);
                
                unsigned _; // do not account for z-bounds
                surf_engine->approximate_zs_in_face(poly_at_x, next_y,
                            SKELETONIZER_THRESHOLD, 0, 0, zs, _, _, false);
            } else {
                //NOTE: dry_run for clip & event lines must always be set !!

                TM_START(tm_lift_on_sil);
                if(dry_run) {
                //! NOTE NOTE NOTE: enable point location if smth goes wrong
                    idcel = surf_engine->locate_pt_at_sil(x, idx, dcel_vec);    
                    features.push_back(Feature_desc(idcel));
                } else {
                    idcel = fdi->idcel;
                    fdi++;
                }
                
                // TODO: can skip z-lifts below/above z-boundaries
                surf_engine->approximate_zs_at_sil(x, idx, dcel_vec[idcel],
                     SKELETONIZER_THRESHOLD, zs);
                TM_STOP(tm_lift_on_sil);
            }

            pt.idcel = idcel;
            add_z_lift(xd, yd, zs, pt);

            stack.push_back(pt);
            cur_y = next_y, n_pts++;

            if(do_z_point) {
                y = next_y; // wouldn't be enough ??

                XTri_OUT("z-point done: y = " << dbl(y) << "; next_y = " <<
                    dbl(next_y) << "; g_top = " << dbl(g_top) << "\n");

                if(z_idx == z_idx_top && (clip_flag & Z_TOP_CLIP))
                    break;
                z_idx++;
            } else {
                if(idx == idx_top && (clip_flag & TOP_CLIP))
                    break;
                idx++;
            }
        }
    }

    sil_idx.push_back(sil_map_l);
    if(what_line & EVENT_LINE)   //! ATTENTION: 2 index vectors for event line
        sil_idx.push_back(sil_map_r);

    XTri_OUT("\n************ vertical line at: " << dbl(x)
        << "; DONE; stack size: " << stack.size() << "\n\n");
}

void create_skeletons(bool use_auto_bounds = true) {

    CA_2 ca_sil = surf_engine->silhouette_curve();
    CA_2 ca_z = surf_engine->z_curve(z_below, z_above);
//     CA_2 ca_sil = surf_engine->z_curve(Rational(0),Rational(0));

    {  typename Kernel_2::Decompose_2 decomp =
            kernel_2->decompose_2_object();

       std::vector< CA_2 > ca1, ca2, ca_common;
       if(decomp(ca_sil, ca_z, std::back_inserter(ca1),
           std::back_inserter(ca2), std::back_inserter(ca_common))) {

            if(ca2.size() > 0)
                ca_z = ca2[0];
            else
                ca_z = kernel_2->construct_curve_2_object()
                    (Polynomial_2(Integer(1)));
            std::cerr << "OUCH: z-curve & silhouette have a non-trivial"
                    " common part..\n";
        }
    }
    g_z_poly = CGAL::to_double(ca_z.polynomial_2());
    g_sil_poly = CGAL::to_double(ca_sil.polynomial_2());

    XTri_OUT("silhouette & z-curves:\n" << ca_sil.polynomial_2() << ",\n");
    XTri_OUT(ca_z.polynomial_2() << "\n\n");

    z_silhouette = z_below - (z_above - z_below) * 7 / 100;

    typename Kernel_2::Construct_curve_pair_2 ccp =
         kernel_2->construct_curve_pair_2_object();

    CPA_2 cpa = ccp(ca_sil, ca_z);
    int n_evts = cpa.number_of_status_lines_with_event();

    if(use_auto_bounds) {
#if 0    
    unsigned n_evts_sil = ca_sil.number_of_status_lines_with_event();
    if(n_evts_sil > 0) { // FIXME FIXME

        const X_coordinate_1& lx = ca_sil.status_line_at_event(0).x(),
            rx = ca_sil.status_line_at_event(n_evts_sil - 1).x();

        _refine_x(lx, SKELETONIZER_THRESHOLD);
        _refine_x(rx, SKELETONIZER_THRESHOLD);

        g_left = lbound_x(lx);
        g_right = ubound_x(rx);

        XTri_OUT("x-bounds are set by silhouette curve..\n");
    }
#endif    
    if(n_evts > 0) { // otherwise take events of CPA

        if(n_evts > 1) {
        //! HERE we take the 1st and the last event of curve pair analysis
        //! which are typically events of the silhouette ..?

            int ileft = 0, iright = n_evts - 1;
            CPA_line_1 line_l, line_r;
            while(ileft < iright) {
                 line_l = cpa.status_line_at_event(ileft);
                 if(line_l.number_of_events() > 0)
                     break;
                 ileft++;
            }
            while(ileft < iright) {
                 line_r = cpa.status_line_at_event(iright);
                 if(line_r.number_of_events() > 0)
                     break;
                 iright--;
            }
            const X_coordinate_1& lx = line_l.x(), rx = line_r.x();

            _refine_x(lx, SKELETONIZER_THRESHOLD);
            _refine_x(rx, SKELETONIZER_THRESHOLD);
            g_left = lbound_x(lx), g_right = ubound_x(rx);
        } else {
            // NOTE NOTE: always take events at intervals ?
            g_left = lbound_x(cpa.status_line_of_interval(0).x());
            g_right = ubound_x(cpa.status_line_of_interval(1).x());
        }
        XTri_OUT("x-bounds are set by CPA\n");

    } else { // no events: take 0th interval
        // always rational
        g_left = lbound_x(cpa.status_line_of_interval(0).x());
        g_right = g_left + Rational(1);
        g_left -= Rational(1);
    }
    Rational diff = g_right - g_left;
    g_left -= diff * enlarge_left, g_right += diff * enlarge_right;

    round_bfi(g_left, ROUND_BFI_PREC, 0);  // round down
    round_bfi(g_right, ROUND_BFI_PREC, 1); // round up
    
    } else { // use absolute bounds
        g_left = enlarge_left;
        g_right = enlarge_right;
    }

    int evt_idx = 0;
    CPA_line_1 line_r, line_l = cpa.status_line_for_x(X_coordinate_1(g_left),
                  CGAL::ZERO);

    line_r = cpa.status_line_for_x(X_coordinate_1(g_right),
                 CGAL::POSITIVE); //! positive: to stop after the last event
    evt_idx = line_l.index(); // starting from the first CPA event index
                // which corresponds to the first event of the silhouette
    // NOTE: n_evts: event after the last one
    n_evts = line_r.index(); // stop at the last event of sil

    XTri_OUT("index range: [" << evt_idx << "; " << (n_evts-1) << "]\n");

    XTri_OUT("g_left: " << g_left << "; " << dbl(g_left) << "; g_right: " <<
        g_right << "; " << dbl(g_right) << "\n\n");

    if(use_auto_bounds) {
    double y_min, y_max;
    bool is_set;

    TM_START(tm_y_range);
    compute_y_range(cpa, y_min, y_max, is_set);
    TM_STOP(tm_y_range);

    if(is_set) {
        XTri_OUT("y_range: " << y_min << "; " << y_max << "\n");
        g_btm = y_min; g_top = y_max;
    } else {
        g_btm = Rational(-1);
        g_top = Rational(1);
    }
    Rational diff = g_top - g_btm;
    g_btm -= diff * enlarge_btm, g_top += diff * enlarge_top;

    round_bfi(g_btm, ROUND_BFI_PREC, 0);  // round down
    round_bfi(g_top, ROUND_BFI_PREC, 1);  // round up
    
    } else { // use_auto_bounds
        g_btm = enlarge_btm;
        g_top = enlarge_top;
    }

    XTri_OUT("g_top: " << g_top << "; " << dbl(g_top) << "; g_btm: " <<
        g_btm << "; " <<  dbl(g_btm) << "\n\n");
    XTri_OUT("z_below: " << dbl(z_below) << "; z_above: " << dbl(z_above) <<
        "\n");

    Rational max_step_x = (g_right - g_left) / ((int)n_slices_x),
            max_step_y = (g_top - g_btm) / ((int)n_slices_y);

    sk_data.reserve(n_evts * 2); // at least twice the number of events
    event_idx.reserve(n_evts); // event index map

    // TODO: use cache for roots at y-boundary ??
    X_coord_vector btm_xs, top_xs, tb_xs, below_xs, above_xs, ba_xs, xs;
    surf_engine->roots_at_y_boundary(g_btm, SKELETONIZER_THRESHOLD, btm_xs,
         false);
    surf_engine->roots_at_y_boundary(g_top, SKELETONIZER_THRESHOLD, top_xs,
         false);

    std::set_union(btm_xs.begin(), btm_xs.end(), top_xs.begin(), top_xs.end(),
        std::back_inserter(tb_xs));

    Polynomial_2 z_poly = ca_z.polynomial_2();
    //TODO: check carefully whether z-polynomial is indeed square-free
    surf_engine->roots_at_y_boundary(g_btm, SKELETONIZER_THRESHOLD, below_xs,
        true, &z_poly);
    surf_engine->roots_at_y_boundary(g_top, SKELETONIZER_THRESHOLD, above_xs,
        true, &z_poly);

    std::set_union(below_xs.begin(), below_xs.end(), above_xs.begin(),
         above_xs.end(), std::back_inserter(ba_xs));

    // merge all clip points: note that there must be no duploicate roots !!!
    std::set_union(ba_xs.begin(), ba_xs.end(), tb_xs.begin(),
         tb_xs.end(), std::back_inserter(xs));

    XTri_OUT("clip points: ");
    for(unsigned jj = 0 ; jj < xs.size(); jj++) {
        XTri_OUT(dbl(xs[jj]) << "  ");
    }
    XTri_OUT("\n");

    typename X_coord_vector::const_iterator xit = xs.begin();
    Rational cur_x = g_left, x = cur_x, next_x;

    int line_idx = 0; // vertical line index (not *event* line index)
    int previous_line = 0;  // keeps a set of flags of the previous line

    while(x < g_right) {
        CPA_line_1 eline;

        XTri_OUT("<--------- line_idx = " << line_idx <<
            "; evt_idx = " << evt_idx << "; cur_x = " << dbl(cur_x) << "\n");

        //bool exact_hit = true;
        int what_line = INTERM_LINE;
        // need to find the next clip line which is > cur_x
        while(xit != xs.end() && *xit < cur_x)
            xit++;

        if(evt_idx < n_evts) {
            eline = cpa.status_line_at_event(evt_idx);
            const X_coordinate_1& eline_x = eline.x();
            _refine_x(eline_x, SKELETONIZER_THRESHOLD);
            next_x = lbound_x(eline_x);

            //! TODO: smth wrong here: we should not take events
            //! which are below the horizontal boundary !!!
            XTri_OUT("next_x: " << dbl(next_x) << "; CPA event_idx: " <<
                    evt_idx << "\n");

            what_line = EVENT_LINE;
            if(xit != xs.end()) {
                XTri_OUT("xit: " << dbl(*xit) << "; next_x: " <<
                        dbl(next_x) << "\n");
                       
                if(*xit == eline_x) { // advance otherwise can loop 4ever
                    what_line |= CLIP_LINE;
                    std::cerr << "clip event & line hit!!\n";
                    //exit(1);
                } else if(xit->compare_distinct(eline_x) == CGAL::SMALLER) {
                    XTri_OUT("&&&&&&& compare_distinct event !!\n");
                    what_line = CLIP_LINE;
                }
            }

        } else {
            if(xit != xs.end() && *xit <= g_right)
                what_line = CLIP_LINE;
            XTri_OUT("next x is set to g_right\n");
            next_x = g_right;
        }
        if(what_line & CLIP_LINE) {
            XTri_OUT("next x is set to lbound_x\n");
            next_x = lbound_x(*xit);
        }
        // FIXME FIXME: remove rounding if problems occur
//         round_bfi(next_x, ROUND_BFI_PREC, 0);

//         XTri_OUT("########### next_x: " << next_x << "\n");

        Rational stepx = (next_x - cur_x) / max_step_x;
        unsigned nsteps = static_cast< unsigned >(
            std::ceil(CGAL::to_double(stepx)));

        XTri_OUT("line_idx: " << line_idx << "; evt_idx: " << evt_idx <<
             "; # nsteps: " << nsteps << "; stepx = " <<
             dbl(stepx) << "; cur_x: " << dbl(cur_x)
             << "; next_x: " << dbl(next_x) << "\n");

        // to guarantee at least 1 intermediate line btw 2 event lines
        if(nsteps == 1) {
            if((what_line & EVENT_LINE))
                XTri_OUT("event line fired!\n");
            
            if((what_line & EVENT_LINE) != 0 &&  evt_idx > 0 &&
                evt_idx < n_evts && (previous_line != CLIP_LINE)) {
                // ensure that previous event line was not a clip line
                nsteps = 2;
            }
            // handling special case when previous line was clip & event
            // line at the same time while this line is a clip line
            if((what_line & CLIP_LINE) != 0 &&
                    (previous_line == (CLIP_LINE | EVENT_LINE))) {
                nsteps = 2;
            }
        }

        stepx = 0;
        if(nsteps != 0)
            stepx = (next_x - cur_x) / (int)nsteps;
        Rational inc(0);
        // in the lst step, iterate until x hits g_right (provided that
        // there is no event / clip line at g_right)
        if(next_x == g_right && (what_line & (EVENT_LINE|CLIP_LINE)) == 0)
            inc = stepx / 2;

        // if not on left boundary, advance x by 1 step
        // to skip preceding event line x. Special case might happen if
        // we start from event / clip line
        if(cur_x != g_left || ((previous_line & (EVENT_LINE|CLIP_LINE)) != 0))
            x = cur_x + stepx;

// NOTE: if we start from left boundary that happens
// to be an event line or clip line => no need for intermediate steps
//         if(cur_x == g_left && (what_line & (EVENT_LINE|CLIP_LINE)) != 0) {
//               x = next_x;
//               XTri_OUT("x = " << dbl(x) << " EVENT fired!\n");
//         }

        bool dry_run = true; // 1st intermediate line collects the information
            // about encountered faces/edges which is then reused by
            // subsequent lines
        for(; x < next_x + inc; x += stepx, line_idx++) {
            Point_stack stack;

            const CA_line_1& iline = ca_sil.status_line_at_exact_x(x),
                    zline = ca_z.status_line_at_exact_x(x);
            run_vertical_line(iline, zline, line_idx, max_step_y, stack,
                 dry_run);
            dry_run = false;
            sk_data.push_back(stack);
        }

        if((what_line & CLIP_LINE) || evt_idx < n_evts) {

            X_coordinate_1 xx = (what_line & EVENT_LINE ? eline.x() : *xit);

            const CA_line_1& iline =
                ca_sil.status_line_at_exact_x(xx),
                zline = ca_z.status_line_at_exact_x(xx);

            if(what_line & CLIP_LINE) {
                XTri_OUT("doing clip line at: " << dbl(*xit) <<
                    "\n");
            } else {
                XTri_OUT("doing event line at: " << dbl(iline.x())
                    << "; line_idx: " << line_idx << "\n");
            }

            bool intersect = (evt_idx >= n_evts ? false :
                  cpa.status_line_at_event(evt_idx).is_intersection());

            Point_stack stack;
            // dry_run is always true for event/clip lines
            run_vertical_line(iline, zline, line_idx, max_step_y, stack, true,
                 what_line, intersect);
            sk_data.push_back(stack);

            if(what_line & EVENT_LINE) {
                event_idx.push_back(line_idx);
                evt_idx++;
            }
            //TASK FIXME clip line indices are not necessary - only for debug
            if(what_line & CLIP_LINE) {
                xit++;
                clip_idx.push_back(line_idx);
            }
            cur_x = next_x, line_idx++;
        }
        // that is, clip line but not an event line
        previous_line = what_line;
    }
}

// catch me if you can))
// void catchme(const Point_index& i1, const Point_index& i2) {
// 
//     Sil_point_2 p1 = sk_data[i1.x][i1.y],
//         p2 = sk_data[i2.x][i2.y];
//         
//    const Point_3d& p1d = points[p1.istart],
//         p2d = points[p2.istart];
//     
//     XTri_OUT("catchme adjacency for " << i1 << " - " << i2
//         << ";\n p[i1]: [" << p1d[0] << "; " << p1d[1] << "]; " <<
//         ";\n p[i2]: [" << p2d[0] << "; " << p2d[1] << "]; \n");
//     try {
//         Adjacency_vector adj1;
//         surf_engine->adjacency_info(p1.idcel, p2.idcel, adj1);
//         std::cerr << "********* adjacency found: ";
//         for(unsigned j = 0; j < adj1.size(); j++)
//             std::cerr << "[" << adj1[j].first << "; " <<
//                 adj1[j].second << "] ";
//         std::cerr << "\n";
//     } catch(...) {
//         throw;
//     }
// }

void check_point(const Sil_point_2& pt1,
        const Sil_point_2& pt2, const Sil_point_2& pt3,
                 int i1, int i2, int i3) {
                 
    if(i1 < pt1.z_skips_below || i1 > pt1.n_zs + pt1.z_skips_below) {
        std::cerr << "************** i1 wrong: " <<
            pt1.z_skips_below << "; " << pt1.n_zs << "; i1 = " << i1 << "\n";
    }

    if(i2 < pt2.z_skips_below || i2 > pt2.n_zs + pt2.z_skips_below) {
        std::cerr << "************** i2 wrong: " <<
            pt2.z_skips_below << "; " << pt2.n_zs << "; i2 = " << i2 << "\n";
    }

    if(i3 < pt3.z_skips_below || i3 > pt3.n_zs + pt3.z_skips_below) {
        std::cerr << "************** i3 wrong: " <<
            pt3.z_skips_below << "; " << pt3.n_zs << "; i3 = " << i3 << "\n";
    }

}

// computes adjacecy btw 2 points and returns its index into adj_data array
int add_adjacency(const Point_index& i1, const Point_index& i2) {

//     XTri_OUT("========== computing adjacency for " << i1 << " - " << i2
//         << "\n");
//         
    int pos;
    Adjacency_vector adj;

    const Sil_point_2& pt1 = sk_data[i1.x][i1.y];
    const Dcel_feature_data& dcel1 = dcel_vec[pt1.idcel],
            dcel2 = dcel_vec[sk_data[i2.x][i2.y].idcel];

//     std::cerr << "dcel1: " << dcel1.first << "; " << dcel1.second << "\n";
//     std::cerr << "dcel2: " << dcel2.first << "; " << dcel2.second << "\n";

    typedef typename Feature_adjacency_map::iterator iterator;
    std::pair< iterator, bool > fi = feature_adj_map.insert(
            typename Feature_adjacency_map::value_type(
                   Addr_pair(dcel1.second, dcel2.second), -1u));

    if(!fi.second) {
        adj = adj_data[fi.first->second];
//         std::cerr << "adj info: [";
//         for(unsigned i = 0; i < adj.size(); i++)
//             std::cerr << "(" <<  adj[i].first << ", " <<  adj[i].second << ") ";
//         std::cerr << "]\n";
        return fi.first->second;
    }

    if(surf_engine->adjacency_info(dcel1, dcel2, adj)) {
        // i1 and i2 belong to the same feature
        unsigned n_zs = pt1.z_skips_below + pt1.n_zs +
                pt1.z_skips_above;
        if(n_zs > 1000) {
            std::cerr << "n_zs is corrupt: " << n_zs << "\n";
            exit(1);
        }
        adj.reserve(n_zs);
        // so that we do not skip other z-lifts
        for(unsigned i = 0; i < n_zs; i++)
            adj.push_back(Adj_pair(i, i));
    }

    std::cerr << "adj info: [";
    unsigned npops = 0;
    bool pops_done = false;
    for(unsigned i = 0; i < adj.size(); i++) {
        std::cerr << "(" <<  adj[i].first << ", " <<  adj[i].second << ") ";
        if(!pops_done && (adj[i].first == -1 || adj[i].second == -1))
            npops++;
        else {
            pops_done = true;
            if(npops != 0)
                adj[i - npops] = adj[i];
        }
    }
    std::cerr << "]\n";

    typename Adjacency_data::iterator it = std::find(adj_data.begin(),
        adj_data.end(), adj);
 
    pos = std::distance(adj_data.begin(), it);
    if(it == adj_data.end())
        adj_data.push_back(adj);

//     adj_map[idx] = pos;
    fi.first->second = pos;
    return pos;
}

void add_tri(Triangles_vector& vtris, const Tri_index& t) {

    vtris.push_back(t);
    
//     const Point_3d& p1 = points[t[0]], p2 = points[t[1]], p3 = points[t[2]];
// 
// //     if(std::abs(p1[2]) > 2.2 || std::abs(p2[2]) > 2.2 ||
// //             std::abs(p3[2]) > 2.2) {
//     std::cerr << "tri: [" << p1[0] << ", " << p1[1] << ", " <<
//         p1[2] << "] - [";
//     std::cerr << p2[0] << ", " << p2[1] << ", " << p2[2] << "] - [";
//     std::cerr << p3[0] << ", " << p3[1] << ", " << p3[2] << "]\n";
//     }
}

//! connects skeleton bones with joints
void connect_bones(const Point_index& i1,
        const Point_index& i2, const Point_index& i3) {

//TODO TASK FIXME etc.. there is an inherent problem with lower z-curves,
// namely that not all lifts starting from 0 must be considered but only
// a subset of them starting from predefined index..

/*    CGAL_assertion(i1.x <= i2.x);
    CGAL_assertion(i2.x <= i3.x);*/
    CGAL_assertion(i1.x != i2.x || i2.x != i3.x);

    const Sil_point_2& pt1 = sk_data[i1.x][i1.y],
        pt2 = sk_data[i2.x][i2.y], pt3 = sk_data[i3.x][i3.y];
 
/*
    printf("%d %d %d %d %d %d\n",
           i1.y, i2.y, i3.y, sk_data[i1.x].size(),
           sk_data[i2.x].size(), sk_data[i3.x].size())*/;

    CGAL_assertion(i3.x < sk_data.size());
    CGAL_assertion(i1.y < sk_data[i1.x].size());
    CGAL_assertion(i2.y < sk_data[i2.x].size());
    CGAL_assertion(i3.y < sk_data[i3.x].size());

    // NOTE NOTE: really strange thing: i2.y < sk_data[i2.x].size()
    //  km42 bug...

    // bool skip = false;

    if(pt1.type == INTERIOR_PT && pt1.n_zs == 0)
        return;
        
    if(pt2.type == INTERIOR_PT && pt2.n_zs == 0)
        return;

    if(pt3.n_zs == 0)
        return;

//     XTri_OUT("lift tri: " << i1 << " - " << i2 << " - " << i3 << "\n");
    Tri_index tri;

    // special points after all z-lifts are used for silhouette
    tri[0] = pt1.istart + pt1.n_zs;
    tri[1] = pt2.istart + pt2.n_zs;
    tri[2] = pt3.istart + pt3.n_zs;
    sil_triangles.push_back(tri);

/*    if(points[pt1.istart][0] > 0)
        return;*/

    int i12 = add_adjacency(i1, i2), i13 = add_adjacency(i1, i3),
        i23 = add_adjacency(i2, i3);

    const Adjacency_vector& adj12 = adj_data[i12], adj13 = adj_data[i13],
        adj23 = adj_data[i23];
        
    unsigned i1s = pt1.istart - pt1.z_skips_below,
             i2s = pt2.istart - pt2.z_skips_below,
             i3s = pt3.istart - pt3.z_skips_below;
    if(pt1.type == INTERIOR_PT) {

        unsigned z_skips = pt1.z_skips_below, zs = pt1.n_zs + z_skips;

//         XTri_OUT("p1.zs = " << pt1.n_zs << "; z_skips = "
//             << z_skips << "\n");

        zs = std::min(static_cast<unsigned>(adj12.size()), zs);
        zs = std::min(static_cast<unsigned>(adj13.size()), zs);

        if(tri_zlifts.size() < zs)
            tri_zlifts.resize(zs);

        for(int i = z_skips, flg = z_skips & 1; i < (int)zs;
                i++, flg ^= 1) {

            if(adj12[i].first != i || adj13[i].first != i) {
                XTri_OUT("FATAL: pt1 wrong triangles!!\n");
                continue;
            }
            
            CGAL_assertion(adj12[i].first == i);
            CGAL_assertion(adj13[i].first == i);

        // alternate CW/CCW vertex order
//             XTri_OUT("connect 1-12-13: " << i << " - " <<  adj12[i].second <<
//                  " - " << adj13[i].second << "; adj12 size: " << adj12.size() << "; adj13 size: " << adj13.size() << "\n");

//             check_point(pt1, pt2, pt3, i, adj12[i].second,
//                         adj13[i].second);

            tri[flg] = i1s + i;
            tri[flg ^ 1] = i2s + adj12[i].second;
            tri[2] = i3s + adj13[i].second;

//             tri_zlifts[i].push_back(tri);
            add_tri(tri_zlifts[i], tri);
        }
    } else if(pt2.type == INTERIOR_PT) {

        unsigned z_skips = pt2.z_skips_below, zs = pt2.n_zs + z_skips;
//         XTri_OUT("p2.zs = " << pt2.n_zs << "; z_skips = "
//             << z_skips << "\n");

        zs = std::min(static_cast<unsigned>(adj12.size()), zs);
        zs = std::min(static_cast<unsigned>(adj23.size()), zs);

        if(tri_zlifts.size() < zs)
            tri_zlifts.resize(zs);

        for(int i = z_skips, flg = z_skips & 1; i < (int)zs; i++,
                flg ^= 1) {

            if(adj12[i].second != i || adj23[i].first != i) {
                XTri_OUT("FATAL: pt2 wrong triangles!!\n");
                continue;
            }

            CGAL_assertion(adj12[i].second == i);
            CGAL_assertion(adj23[i].first == i);
//             CGAL_assertion(std::find(adj13.begin(), adj13.end(),
//                  std::make_pair(adj12[i].first, adj23[i].second))
//                     != adj13.end());

//             check_point(pt1, pt2, pt3, adj12[i].first, i,
//                         adj23[i].second);

            tri[flg] = i1s + adj12[i].first;
            tri[flg ^ 1] = i2s + i;
            tri[2] = i3s + adj23[i].second;

//             tri_zlifts[i].push_back(tri);
            add_tri(tri_zlifts[i], tri);
        }

    } else {
        CGAL_assertion(pt3.type == INTERIOR_PT);
        unsigned z_skips = pt3.z_skips_below, zs = pt3.n_zs + z_skips;

        zs = std::min(static_cast<unsigned>(adj13.size()), zs);
        zs = std::min(static_cast<unsigned>(adj23.size()), zs);

        if(tri_zlifts.size() < zs)
            tri_zlifts.resize(zs);
//         XTri_OUT("p3.zs = " << pt3.n_zs << "; z_skips = " << z_skips << "\n");
                
        for(int i = z_skips, flg = z_skips & 1; i < (int)zs; i++,
                flg ^= 1) {

            if(adj13[i].second != i || adj23[i].second != i) {
                XTri_OUT("FATAL: pt3 wrong triangles!!\n");
                continue;
            }

        // means that points 1 & 2 are adjacent to the ith lift of point 3
            CGAL_assertion(adj13[i].second == i);
            CGAL_assertion(adj23[i].second == i);

//             check_point(pt1, pt2, pt3, adj13[i].first, 
//                         adj23[i].first, i);

            tri[flg] = i1s + adj13[i].first;
            tri[flg ^ 1] = i2s + adj23[i].first;
            tri[2] = i3s + i;

//             tri_zlifts[i].push_back(tri);
            add_tri(tri_zlifts[i], tri);
        }
    }
}

//! skeletonizes region between x1 = \c stacki.x y1 = [\c stacki.y; \c lasti_y]
//! and x2 = \c stackj.x y2 = [\c stackj.y; \c lastj_y]
//! \c j_is_event indicates that j-th stack is an event line stack, i.e.,
//! one has to account for silhouette points with 0 outgoing branches
void skeletonize_region(const Point_index& stacki, unsigned lasti_y,
     const Point_index& stackj, unsigned lastj_y, bool j_is_event = false) {

    XTri_OUT("---->> region: { " << stacki.x << " [" << stacki.y << "; " <<
        lasti_y << "] } and { " << stackj.x <<  " [" << stackj.y << "; "
            << lastj_y << "] }\n");

    // for left stripe, stacki.x < stackj.x
    // for right stripe, stacki.x > stackj.x
    // jth stack is always event line stack by construction
    bool reversed = (stacki.x > stackj.x);
    if(stackj.y == lastj_y) { // trivial region -> only for jth range

        XTri_OUT("==== trivial region at: " << stackj << "\n");
        if(stacki.y == lasti_y) {
            XTri_OUT("FATAL: degenerate region..\n");
        }

        Point_index triv = stackj, pi = stacki,
            pip1(stacki.x, stacki.y + 1);
        for(; pi.y < lasti_y; pi.y++, pip1.y++) {

            if(!reversed) // stacki < stackj
                connect_bones(pi, pip1, triv); // points in CCW order
            else
                connect_bones(triv, pip1, pi); // points in CW order
        }
        XTri_OUT("==== trivial region done\n");
        return;
    }

    if(stacki.y == lasti_y) { // trivial region -> only for jth range

        XTri_OUT("==== trivial region at: " << stacki << "\n");

        Point_index triv = stacki, pj = stackj,
            pjp1(stackj.x, stackj.y + 1);
        for(; pj.y < lastj_y; pj.y++, pjp1.y++) {

            if(reversed) // stacki < stackj
                connect_bones(pj, pjp1, triv); // points in CCW order
            else
                connect_bones(triv, pjp1, pj); // points in CW order
        }
        XTri_OUT("==== trivial region done\n");
        return;
    }
    
    // non-trivial region
    const Point_stack& si = sk_data[stacki.x], sj = sk_data[stackj.x];
    Point_index pti = stacki, ptj = stackj, next;

    // (pti.y + 1)th points come directly after n_zs pti.y points
    Point_vec_3d::const_iterator yi = points.begin() + si[pti.y].istart,
            yj = points.begin() + sj[ptj.y].istart,
            yip1 = yi + si[pti.y].n_zs + 1, yjp1 = yj + sj[ptj.y].n_zs + 1;

    // rangei.y is the last index of i
    XTri_OUT("==== normal region: " << stacki.x << " - " << stackj.x << "\n");
    while(pti.y < lasti_y || ptj.y < lastj_y) {

        int inc_j = 1; // whether to advance i or j

        bool fired = false;
        if(ptj.y ==  lastj_y || (pti.y < lasti_y &&
                std::abs((*yip1)[1] - (*yj)[1]) <
                     std::abs((*yi)[1] - (*yjp1)[1]))) {
            inc_j = 0; // connect i+1 to j
        //NOTE alternatively can provide a flag 'rangej contains sil points
        // with 0 outgoing branches'
       // NOTE one is only allowed to step at the last point of the ith stack
       // if no silhouette points (with 0 outgoing
        // branches) are left in [ptj.y..rangej.y-1] interval
        // condition is to be checked only if rangei.y is a point on sil
            if(j_is_event && pti.y == lasti_y - 1 && si[lasti_y].type) {
                // scan over all points at jth stack
                unsigned j = ptj.y;
                while(!sj[j].type && j < lastj_y)
                    j++;
                if(j < lastj_y) {// point on sil was found-> reverse inc_j
                    inc_j = 1, fired = true;
                    XTri_OUT("j < rangej.y fired\n");
                }
            }
        }

        if(!fired) { // conditions are mutually exclusive
            unsigned inext = pti.y + (inc_j ^ 1), jnext = ptj.y + inc_j;
            if((inext < lasti_y || jnext < lastj_y) &&
                    si[inext].type && sj[jnext].type) {
                XTri_OUT("si[inext].type && sj[jnext].type\n");
                inc_j ^= 1;
            }
        }
        Point_index next = (inc_j ? ptj : pti);
        next.y++;

        if(!reversed)
            connect_bones(pti, next, ptj);
        else
            connect_bones(ptj, next, pti);

        if(inc_j) {
            ptj = next;
            yj = yjp1, yjp1 += sj[ptj.y].n_zs + 1;
        } else {
            pti = next;
            yi = yip1, yip1 += si[pti.y].n_zs + 1;
        }
    }
    XTri_OUT("==== normal region done\n");
}

void triangulate_silhouette() {

    unsigned line_idx = -1u;
    Index_vector::const_iterator git = event_idx.begin();
    std::vector< Index_vector >::const_iterator sil_map_it = sil_idx.begin();

    XTri_OUT("///////////////////////////// event line indices: ");
    for(; git != event_idx.end(); git++)
        XTri_OUT(*git << " ");

    XTri_OUT("\n\n clip line indices: ");
    for(git = clip_idx.begin(); git != clip_idx.end(); git++)
        XTri_OUT(*git << " ");

    git = event_idx.begin();

    for(; ;) {

        unsigned next_line = (git == event_idx.end() ? sk_data.size() :
                 *git);

        XTri_OUT("\n^^^^^^^^^^^^^doing intermediate triangles; line_idx: " <<
            (int)line_idx << "; next_line: " << next_line << "\n");

        for(unsigned li = line_idx + 1; li + 1 < next_line; li++) {

            unsigned s1_sz = sk_data[li].size(),
                s2_sz = sk_data[li + 1].size();

            const Index_vector& ln1_sil_map = *sil_map_it++,
                    ln2_sil_map = *sil_map_it;

            Index_vector::const_iterator ln1_sil_it = ln1_sil_map.begin(),
                ln2_sil_it = ln2_sil_map.begin();

/*            XTri_OUT("ln1_sil_map: [");
            for(Index_vector::const_iterator it1 = ln1_sil_map.begin();
                    it1 != ln1_sil_map.end(); it1++) {
                XTri_OUT(*it1 << " ");
            }

            XTri_OUT("]\nln2_sil_map: [");
            for(Index_vector::const_iterator it1 = ln1_sil_map.begin();
                    it1 != ln1_sil_map.end(); it1++) {
                XTri_OUT(*it1 << " ");
            }
            XTri_OUT("]\n")*/;

            unsigned i1 = 0, i2 = 0, next_i1 = i1, next_i2 = i2;

            if(ln1_sil_map.size() != 0 && *ln1_sil_it == 0 &&
                    curve_on_one_side(li, li + 1)) {
                XTri_OUT("ln1_sil_it++ fired\n");
                ln1_sil_it++;
            }
        //TODO shouldn't these conditions be mutually exclusive ??
            // if side_s has branches, [li + 1; 1] must be intermediate pt
            /*else*/ if(ln2_sil_map.size() != 0 && *ln2_sil_it == 0 &&
                    curve_on_one_side(li + 1, li)) {
                XTri_OUT("ln2_sil_it++ fired\n");
                ln2_sil_it++;
            }

//             XTri_OUT("cmp1: " << (ln1_sil_it == ln1_sil_map.end()) <<
//                     "; s1_sz: " << s1_sz << "\n");
// 
//             XTri_OUT("cmp2: " << (ln2_sil_it == ln2_sil_map.end()) <<
//                     "; s2_sz: " << s2_sz << "\n");

            while(i1 + 1 < s1_sz || i2 + 1 < s2_sz) {

                next_i1 = (ln1_sil_it == ln1_sil_map.end() ?
                    s1_sz - 1 : *ln1_sil_it++);

                next_i2 = (ln2_sil_it == ln2_sil_map.end() ?
                    s2_sz - 1 : *ln2_sil_it++);

                skeletonize_region(Point_index(li, i1), next_i1,
                    Point_index(li + 1, i2), next_i2);
                i1 = next_i1, i2 = next_i2;
            }
        }

        if(git == event_idx.end())
            break;

        line_idx = *git++;
        XTri_OUT("^^^^^^^^^^^^^^^^^^ line_idx: " << line_idx << "\n");

        unsigned s1_sz = sk_data[line_idx].size();

        for(unsigned which = 0; which < 2; which++) {

            Index_vector side_sil_map, event_sil_map;
            if(which == 0) {
                if(line_idx == 0)  { // no left lift here
                    sil_map_it++; // need to increment iterator
                    continue;
                }
                side_sil_map = *sil_map_it++;
                event_sil_map = *sil_map_it++;
                XTri_OUT("doing left line: \n");
            } else {
                if(sil_map_it ==  sil_idx.end())
                    break;
                event_sil_map = *sil_map_it++;

                if(sil_map_it ==  sil_idx.end())
                    break;
                side_sil_map = *sil_map_it;
                XTri_OUT("doing right line: \n");
            }

            Index_vector::const_iterator event_sil_it = event_sil_map.begin(),
                side_sil_it = side_sil_map.begin();

            unsigned side_line = line_idx + which*2 - 1,
                side_sz = sk_data[side_line].size();

            unsigned i1 = 0, i2 = 0, next_i1 = i1, next_i2 = i2;

//             if(side_sil_map.size() != 0)
//                 XTri_OUT("1st sil point on " << side_line << " is " <<
//                         *side_sil_it << "\n");

            if(side_sil_map.size() != 0 && *side_sil_it == 0 &&
                    curve_on_one_side(side_line, line_idx)) {
                side_sil_it++;
            // these conditions are exclusive
            } else if(event_sil_map.size() != 0 && *event_sil_it == 0 &&
                    curve_on_one_side(line_idx, side_line)) {
                event_sil_it++;
            }

            // this condition never fires
            while(i2 + 1 < side_sz || i1 + 1 < s1_sz) {

//                 XTri_OUT("i1: " << i1 << "; i2: " << i2 <<
//                     "; s1_sz: " << s1_sz << "; side_sz: " << side_sz << "\n");
                next_i1 = (event_sil_it == event_sil_map.end() ?
                    s1_sz - 1 : *event_sil_it++);

                next_i2 = (side_sil_it == side_sil_map.end() ?
                    side_sz - 1 : *side_sil_it++);

                // j_is_event is set
                skeletonize_region(Point_index(side_line, i2), next_i2,
                    Point_index(line_idx, i1), next_i1, true);

                i1 = next_i1, i2 = next_i2;
            }
        } // end of 'which' loop

    } // end of events loop

    // now we merge index ranges for different z-lifts in a single
    // array 'triangles' and also populating zlift_ranges array

    zlift_ranges.resize(tri_zlifts.size());
    for(unsigned i = 0, k = 0; i < tri_zlifts.size(); i++) {

        unsigned sz = tri_zlifts[i].size();
        triangles.resize(triangles.size() + sz);
        std::copy(tri_zlifts[i].begin(), tri_zlifts[i].end(),
                        triangles.begin() + k);
        k += sz;
        zlift_ranges[i] = k;
    }
}

void clear() {
    is_valid = false;

    tri_zlifts.clear();
    zlift_ranges.clear();
    triangles.clear();
    sil_triangles.clear();
    sk_data.clear();
    adj_data.clear();
    adj_map.clear();
    feature_adj_map.clear();
    dcel_vec.clear();
    points.clear();
    features.clear();
    event_idx.clear();
    clip_idx.clear();
    sil_idx.clear();

    g_start_idx = 0;
}

private:
    //! private data

    //! global skeleton data
    Skeleton_data sk_data; 
    //! adjacency info
    Adjacency_data adj_data;
    Adjacency_map adj_map;

    Feature_adjacency_map feature_adj_map;
    Dcel_feature_vec dcel_vec;

    //! surface triangles & silhouette triangles
    Triangles_vector triangles, sil_triangles;
    //! temporary storage to differentiate z-lifts
    std::vector< Triangles_vector > tri_zlifts;
    //! index ranges for different z-lifts
    std::vector< unsigned > zlift_ranges;
    
    //! points in 3D space (double precision)
    Point_vec_3d points;   //! set of lifted points in 3D
    //! supporting polynomials for z-curve and silhouette curve
    Poly_double_2 g_z_poly, g_sil_poly; 

    //! maps event stacks in Skeleton_data
    std::vector< unsigned > event_idx;
    //! indices of clip-points (needed for debug only)
    std::vector< unsigned > clip_idx;

    //! maps to points on silhouette: 1 event line is assigned to 2 index
    //! vectors - immediately to the left and to the right of it
    std::vector< Index_vector > sil_idx;
    std::vector< Feature_desc > features;

    //! # of slices in x / y direction
    unsigned n_slices_x, n_slices_y, g_start_idx;

    //! boundaries: x, y, z respectively
    Rational g_left, g_right, g_btm, g_top, z_below, z_above;
    Rational z_silhouette; // z-coordinate for silhouette plane

    Rational enlarge_btm, enlarge_top;
    Rational enlarge_left, enlarge_right;

    //! pointer to a surface engine
    Surface_analysis_interface *surf_engine;

    //! timers for debugging
#ifdef XTri_USE_OUT
    CGAL::Timer tm_lift_on_sil, tm_lift_in_face,
        tm_create_skeletons, tm_tri_silhouette, tm_load_surface,
            tm_y_range, tm_arcavoid_in_face;
#endif
    //! instance of kernel_2 for passed caching
    const Kernel_2 *kernel_2;
    //! indicates whether triangulation is valid
    bool is_valid;
}; // Surface_skeletonizer_impl

} // namespace CGAL

#endif // CGAL_SURFACE_SKELETONIZER_IMPL_H
