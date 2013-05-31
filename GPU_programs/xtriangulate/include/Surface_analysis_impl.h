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

#ifndef CGAL_SURFACE_ANALYSIS_IMPL_H
#define CGAL_SURFACE_ANALYSIS_IMPL_H

#include <CGAL/basic.h>

#define CGAL_USE_FILTERED_CKvA_2 0
#define CGAL_XTri_USE_ARCAVOID 1

#if CGAL_USE_FILTERED_CKvA_2
#include <CGAL/Curved_kernel_via_analysis_2/Filtered_curved_kernel_via_analysis_2_impl.h>
#else
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>
#endif

#include <CGAL/Algebraic_kernel_3/Algebraic_surface_3.h>
#include <CGAL/Algebraic_kernel_3/Algebraic_surface_3_z_at_xy_isolator_traits_base.h>
#include <CGAL/Algebraic_kernel_3/Algebraic_surface_3_z_at_xy_isolator_traits.h>
#include <CGAL/Algebraic_kernel_3/IO/Algebraic_surface_3_iostream.h>

#include <CGAL/Polynomial_type_generator.h>

#ifndef XTri_OUT
#define XTri_USE_OUT
#ifdef XTri_USE_OUT
#define XTri_OUT(x) std::cerr << x
#else
#define XTri_OUT(x) static_cast<void>(0)
#endif
#endif

namespace CGAL {


template < class AlgebraicKernel_2 >
class Surface_analysis_impl {

public:
    //! kernels
    typedef AlgebraicKernel_2 Kernel_2;
    typedef typename Kernel_2::Algebraic_kernel_d_1 Kernel_1;

    //! number types
    typedef typename Kernel_2::Arithmetic_kernel AT;
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;
    typedef typename AT::Bigfloat Bigfloat;
    typedef typename AT::Bigfloat_interval BFI;

    //! polynomials & traits
    typedef CGAL::Polynomial< Integer > Polynomial_1;
    typedef CGAL::Polynomial< Polynomial_1 > Polynomial_2;
    typedef CGAL::Polynomial< Polynomial_2 > Polynomial_3;
    typedef typename CGAL::Polynomial_type_generator< Rational, 3 >::Type
            Poly_rat_3;
    typedef CGAL::Polynomial_traits_d< Polynomial_3 > PT_3;
    
    //! analyses
    typedef typename Kernel_2::Curve_analysis_2 CA_2;

    //! coordinates
    typedef typename Kernel_2::Algebraic_real_1 X_coordinate_1;
    typedef typename Kernel_2::Algebraic_real_2 Xy_coordinate_2;

    typedef std::vector< double > Double_vector;
    typedef std::vector< X_coordinate_1 > X_coord_vector;

    //! adjacency stuff
    typedef CGAL::Adjacencies_3::Adjacency_vector Adjacency_vector;
    typedef CGAL::Adjacencies_3::Adjacency_pair Adj_pair;

private:
    //! private stuff:

    //! arrangement traits
#if CGAL_USE_FILTERED_CKvA_2
#warning "using filtered kernel"
    typedef CGAL::Curved_kernel_via_analysis_2< Kernel_2 > AFCK_2;
    typedef CGAL::Filtered_curved_kernel_via_analysis_2< AFCK_2 > CKvA_2;
#else
    typedef CGAL::Curved_kernel_via_analysis_2< Kernel_2 > CKvA_2;
#endif

    //! surface definitions and CAD traits
    typedef CGAL::Algebraic_surface_3< Integer > Surface_3;
    typedef CGAL::Algebraic_surface_3_z_at_xy_isolator_traits<
        CKvA_2, Surface_3 > Cad_traits;

    //! restricted CAD definitions
    typedef CGAL::Create_restricted_cad_3< Cad_traits >
        Restricted_cad_creator_3;
    typedef typename Restricted_cad_creator_3::Restricted_cad_3
        Restricted_cad_3;
    typedef CGAL::Restricted_cad_3_accessor< Restricted_cad_3 > RC3_accessor;
    typedef typename Restricted_cad_3::Z_at_xy_isolator Isolator;
    typedef typename Restricted_cad_3::Z_stack Z_stack;
    typedef typename Restricted_cad_3::Arrangement_traits_2 Arr_traits_2;

    //! point of an arrangement
    typedef typename Arr_traits_2::Point_2 Point_2;

    //! arrangement features (vertices, edges and faces handles)
    typedef typename Restricted_cad_3::Vertex_const_handle VCH;
    typedef typename RC3_accessor::Halfedge_const_handle ECH;
    typedef typename Restricted_cad_3::Face_const_handle FCH;

public:
    //! constructors

    Surface_analysis_impl(Kernel_2 *_kernel_2 = 0) :
        kernel_2(0), is_valid(false) {

        kernel_2 = &kernel_2_inst();
    }

public:
    //! interface

const Kernel_2& kernel_2_inst() const {
    return Arr_traits_2::instance().kernel();
}

bool surface_valid() {
    return is_valid;
}

void load_surface(const Polynomial_3& in) {

    is_valid = false;
    Polynomial_3 poly3 = in;

    // too big number is bad because leads to very lengthy/ coefficients
    unsigned n_xz = 2, n_xy = 1;
    unsigned what = 0; //! what ??!

    //int i = 1, j = 0; // check genericity conditions
    while(1) {
        bool make_sqfree = false;
Lbegin:        
        try {
            base_surf = Surface_3(poly3);
            std::cerr << "computing resultant.. \n";
            base_sil_poly = base_surf.resultant_f_fz();
        }
        catch(...) {
            if(make_sqfree) {
                // already checked for square-freeness
                std::cerr << "failed to make square-free.. \n";
                return;
            }
            make_sqfree = true;
            std::cerr << "Zero resultant..making square-free\n";
            poly3 = CGAL::make_square_free(poly3);
            goto Lbegin;
        }

        std::cerr << "making square-free.. \n";
        if(base_sil_poly.content().degree() == 0)
            break;

        unsigned n = (what == 0 ? n_xz++ : n_xy++);
        poly3 = affine_transform(poly3, n, what);
        what ^= 1;

        std::cerr << "o|o nontrivial content.. trying once again\n";
//         std::cerr << "content: " << base_sil_poly.content() << "\n";
//         poly = shear_surface(poly, i, j);
//         (i < j ? i++ : j++);
    }

    base_poly = base_surf.f();
    std::cerr << "creating cad for: " << base_poly << "\n";

    base_cad = Restricted_cad_creator_3()(base_surf);
    base_poly_swapped = typename PT_3::Swap()(base_poly, 0, 2); // swap x & z

    base_sil_poly = CGAL::make_square_free(base_sil_poly);
    std::cerr << "silhouette: " << base_sil_poly << " done\n";

    typename Kernel_2::Construct_curve_2 cc =
         kernel_2->construct_curve_2_object();
    base_ca_sil = cc(base_sil_poly);

    is_valid = true;
    std::cerr << "Surface loaded..\n";
}

/*bool overlay_surface(const char *filename) {

    Polynomial_2 sil;
    Surface_3 surf_bounds[2] = { base_surf };

    if(!_load_surface(filename, surf_bounds[1], sil))
        return false;

    base_cad = Restricted_cad_creator_3()(surf_bounds, surf_bounds + 2);

    Kernel_2::Construct_curve_2 cc =
         kernel_2_inst().construct_curve_2_object();
    base_ca_sil = cc(base_sil_poly * sil); // hack..

    std::cout << "Surface overlayed with: " << surf_bounds[1] << "\n; " <<
        "; resulting sil: " << base_ca_sil.polynomial_2() << "\n";

    return true;
}*/

CA_2 silhouette_curve() const {

    if(!is_valid)
        throw "Surface analysis is not valid";

    return base_ca_sil;
}

//! returns surface slice at z = z_below * z = z_above
CA_2 z_curve(Rational& z_below, Rational& z_above) {
    
    if(!is_valid)
        throw "ouch!";

    typedef CGAL::Fraction_traits< Rational > F_traits;
    typename F_traits::Numerator_type z_num;
    typename F_traits::Denominator_type z_den;
    typename F_traits::Decompose decompose;
    typename PT_3::Evaluate_homogeneous heval3;

    Rational dist = z_above - z_below;

    Polynomial_2 below_p, above_p, *ppoly = &below_p;
    Rational *z = &z_below, inc = -dist / 100;

    for(int i = 0; i < 2; i++) {
        while(1) {
            decompose(*z, z_num, z_den);
            *ppoly = heval3(base_poly, (z_num), (z_den));
            if(!CGAL::is_zero(*ppoly))
                break;
        
            std::cerr << "vanishing z-curve.. repeating\n";
            *z += inc;
        }
        *ppoly = CGAL::make_square_free(*ppoly);
        z = &z_above, inc = dist / 100;
        ppoly = &above_p;
    }

    typename Kernel_2::Decompose_2 decomp =
        kernel_2->decompose_2_object();

    typename Kernel_2::Construct_curve_2 cc =
         kernel_2->construct_curve_2_object();

    CA_2 c1 = cc(below_p), c2 = cc(above_p);
    if(c1.is_identical(c2)) {
        return c1;
    }

    std::vector< CA_2 > ca1, ca2, ca_common;

    if(decomp(c1, c2, std::back_inserter(ca1),
           std::back_inserter(ca2), std::back_inserter(ca_common))) {
        std::cerr << "z-curves have common part..ouch!\n";

        if(ca1.size() != 0) {
            return cc(ca1[0].polynomial_2() * above_p);
        }
        if(ca2.size() != 0) {
            return cc(below_p * ca2[0].polynomial_2());
        }
    }
    return cc(c1.polynomial_2() * c2.polynomial_2());
}


void compute_y_range(double& y_min, double& y_max, bool& is_set) {

    is_set = false;
    for(VCH vit = base_cad.vertices_begin();
            vit != base_cad.vertices_end(); vit++) {

        const std::pair< double, double>& xy =
            vit->point().xy().to_double();

        XTri_OUT("; xy: " <<
            xy.first << "; " << xy.second << "\n");
        if(!is_set || xy.second < y_min)
            y_min = xy.second;
        if(!is_set || xy.second > y_max) {
            y_max = xy.second;
            is_set = true;
        }
    }
}

void dump_arr() const {
        std::set<unsigned> idds;

        std::cerr << "EDGES:\n";

        for(typeof(base_cad.edges_begin()) cit =
                base_cad.edges_begin(); cit != base_cad.edges_end(); cit++) {

            if(
                idds.find(cit->curve().curve().id()) == idds.end()) {

                std::cerr << "SUPPORT@ " << cit->curve().curve().id() << "    " <<
                    cit->curve().curve().polynomial_2() << "\n\n";
                idds.insert(cit->curve().curve().id());
            }
            std::cerr << cit->curve() << "\n\n";
        }

        std::cerr << "VERTICES:\n";
        for(VCH vit =
                base_cad.vertices_begin(); vit != base_cad.vertices_end(); vit++) {

            std::cerr << vit->point();
            if(vit->is_isolated())
                std::cerr << " * isolated *\n\n";
            else
                std::cerr << "\n\n";
        }
}

// TODO: isolated points are removed from the arrangement..
unsigned locate_pt_at_sil(const X_coordinate_1& x,
        int arcno, Dcel_feature_vec& vec) {

    CGAL::Object obj;
    if(!is_valid)
        throw "Surface analysis is not valid";

    typename Arr_traits_2::Construct_point_2 cp_2 = Arr_traits_2::instance().
            construct_point_2_object();

//     XTri_OUT("locate pt on sil: " << CGAL::to_double(x) << "; "
//             << arcno << "\n" << silhouette_curve().polynomial_2() << "\n");

    Point_2 pt = cp_2(x, silhouette_curve(), arcno);
    obj = base_cad.locate(pt);

    VCH v;
    if(CGAL::assign(v, obj)) {
        vec.push_back(Dcel_feature_data(CGAL::VERTEX,
                static_cast< const void * >(v.ptr())));

    } else {
        // for isolated points point location fails..
        //std::cerr << "point location for: " << pt << "\n\n";

        ECH he;
        CGAL_assertion_code(bool check = )
        CGAL::assign(he, obj);
        CGAL_assertion(check);

        vec.push_back(Dcel_feature_data(CGAL::EDGE,
            static_cast< const void * >(he.ptr())));
    }
    return vec.size() - 1;
}

unsigned locate_pt_in_face(const X_coordinate_1& x,
        const Rational& y, Dcel_feature_vec& vec) {

    RC3_accessor acc(base_cad);
    Point_2 pt = acc.construct_point_with_rational_y(x, y);
    CGAL::Object obj = base_cad.locate(pt);

    FCH f;
    CGAL_assertion_code(bool check = )
    CGAL::assign(f, obj);
    CGAL_assertion(check);

    vec.push_back(Dcel_feature_data(CGAL::FACE,
            static_cast< const void * >(f.ptr())));

    return vec.size() - 1;
}

// CGAL::circulator_size(edge_it->face()->outer_ccb())
// dump_face(edge_it->face());
//  dump_face(edge_it->twin()->face());

void approximate_zs_at_sil(const X_coordinate_1& x, int arcno,
       const Dcel_feature_data& dcel, const double& threshold,
                  Double_vector& zs) {

    Isolator isl;

    if(dcel.first == CGAL::VERTEX) {
        VCH v = VCH(static_cast< typename VCH::pointer>
            (dcel.second));
        isl = base_cad.z_stack(v).isolator(base_surf);

    } else {
        CGAL_assertion(dcel.first == CGAL::EDGE);
        ECH he = ECH(static_cast< typename ECH::pointer >
            (dcel.second));

        typename Arr_traits_2::Construct_point_2 cp_2 =
                Arr_traits_2::instance().construct_point_2_object();
        Point_2 pt = cp_2(x, silhouette_curve(), arcno);
        isl = base_cad.isolator_at_point_on_halfedge(pt, base_surf, he);

//         CGAL::Nk nk = heh->data()->_nk(surface);
    }

    zs.reserve(isl.number_of_real_roots());
    for(int k = 0; k < isl.number_of_real_roots(); k++) {

        while(isl.right_bound(k) - isl.left_bound(k) > threshold)
            isl.refine_interval(k);

        double z = CGAL::to_double(isl.left_bound(k));
        zs.push_back(z);
    }
}

Polynomial_2 poly_at_rational_x(const Rational& x) {

    typedef CGAL::Fraction_traits< Rational > F_traits;
    typename F_traits::Numerator_type x_num;
    typename F_traits::Denominator_type x_den;
    typename F_traits::Decompose decompose;

    decompose(x, x_num, x_den);

    typename PT_3::Evaluate_homogeneous heval3;
    return heval3(base_poly_swapped, (x_num), (x_den));
}

#if 0
void exact_zs_in_face(const Xy_coordinate_2& xy,
        unsigned idcel,
        const Rational& threshold, const Rational& z_below,
        const Rational& z_above, Double_vector& zs) const {

    Arr_traits_2::Construct_point_2 cp_2 = Arr_traits_2::instance().
            construct_point_2_object();
    Point_2 pt = cp_2(xy.x(), xy.curve(), xy.arcno());

    CGAL_assertion(feature_vec[idcel].first == CGAL::FACE);
    FCH fh = FCH(static_cast< FCH::pointer >(feature_vec[idcel].second));

    Isolator isl = base_cad.isolator_at_point_in_face(pt, base_surf, fh);

    zs.reserve(isl.number_of_real_roots());
    for(int k = 0; k < isl.number_of_real_roots(); ) {

        double z = CGAL::to_double(isl.left_boundary(k));
        zs.push_back(z);
        k++;
        continue; // do u need this in fact ??

        Rational lbd = isl.left_boundary(k), rbd = isl.right_boundary(k);
        if(rbd < z_below || lbd > z_above) {
//             XTri_OUT("skipping root at lbd = " << dbl(lbd) << "; rbd = " <<
//                 dbl(rbd) << "\n");
            k++; continue;
        }

        if(rbd - lbd <= threshold &&
                ((lbd <= z_below && rbd >= z_below) ||
                    (lbd <= z_above && rbd >= z_above))) {
            double z = CGAL::to_double(isl.left_boundary(k));
            zs.push_back(z);
            k++;
        } else

        if(lbd >= z_below && rbd <= z_above) {

            while(isl.right_boundary(k) - isl.left_boundary(k) > threshold)
                isl.refine_interval(k);
            double z = CGAL::to_double(isl.left_boundary(k));
            zs.push_back(z);
            k++;
        } else {

            isl.refine_interval(k);
        }
    }
}
#endif

//! returns # of roots skipped on z_below
//! \c no_z_bounds - do not account for \c z-below/above bounds
void approximate_zs_in_face(const Polynomial_2& poly_at_x,
        const Rational& y, const Rational& threshold,
        const Rational& z_below, const Rational& z_above,
         Double_vector& zs, unsigned& n_skips_below,
          unsigned& n_skips_above, bool use_z_bounds) {

    if(!is_valid)
        throw "Surface analysis is not valid";

    typedef CGAL::Fraction_traits< Rational > F_traits;
    typename F_traits::Numerator_type y_num;
    typename F_traits::Denominator_type y_den;
    typename F_traits::Decompose decompose;
    decompose(y, y_num, y_den);

    typename CGAL::Polynomial_traits_d< Polynomial_2 >::Evaluate_homogeneous
        heval2;
    Polynomial_1 peval1 = heval2(poly_at_x, y_num, y_den);
    
//     XTri_OUT("poly: " << peval1 << "; " << CGAL::to_double(y) << "\n");

    n_skips_below = 0, n_skips_above = 0;
#if CGAL_XTri_USE_ARCAVOID

    typedef CGAL::internal::Bitstream_coefficient_kernel< Integer >
        IBCK;

    typedef CGAL::Arcavoid< IBCK, CGAL::Arcavoid_real_root_isolator_tag >
        Arcavoid;

    long prec = CGAL::set_precision (BFI(), 60);
    Arcavoid isolator(CGAL::Square_free_arcavoid_tag(), peval1);
//     typedef typename Arcavoid::Algebraic_complex_1_iterator
//         Algebraic_complex_1_iterator;

//     std::cerr << "approximate zs_in_face: " <<
//             CGAL::to_double(z_below) << "; " <<
//                 CGAL::to_double(z_above) << "; " << "\n";

    int n_roots = isolator.number_of_real_roots();

    zs.reserve(n_roots);
    for(int i = 0; i < n_roots; i++) {
        typename Arcavoid::Bound r = isolator[i].center().real();
//         std::cerr << i << ": " << CGAL::to_double(r) << "; " << Rational(r)
// << " (" << CGAL::to_double(isolator[i].radius()) << ") " <<
//                 n_skips << "\n";
        if(use_z_bounds) {
            if(r < z_below) {
                n_skips_below++;
                continue;
            } else if(r > z_above) {
                n_skips_above++;
                continue;
            }
        }
        zs.push_back(CGAL::to_double(r));
    }
    
    CGAL::set_precision (BFI(), prec);
#else

    X_coord_vector roots;
    unsigned n_roots = univariate_roots(peval1, threshold, roots);

    typename Kernel_1::Algebraic_real_traits::Lower_bound lbound_x;
    
    zs.reserve(n_roots);
    for(unsigned i = 0; i < n_roots; i++) {

        const X_coordinate_1& r = roots[i];
        if(use_z_bounds) {
            if(r < z_below) {
                n_skips_below++;
                continue;
            } else if(r > z_above) {
                n_skips_above++;
                continue;
            }
        }
        zs.push_back(CGAL::to_double(lbound_x(r)));
    }
#endif // CGAL_XTri_USE_ARCAVOID

#if 0
    Z_stack z_stack;
    //z_stack = base_cad.z_stack(f, base_surf).first;
    //z_stack = _m_cad.z_stack_at(pt).first;

//isl = base_cad.isolator_at_point_on_halfedge(pt, base_surf, he);
    Isolator isl = base_cad.isolator_at_point_in_face(pt, base_surf, f);
            //z_stack.isolator(base_surf);
    zs.reserve(isl.number_of_real_roots());
    for(int k = 0; k < isl.number_of_real_roots(); k++) {

        while(isl.right_bound(k) - isl.left_bound(k) > threshold)
            isl.refine_interval(k);

        double z = CGAL::to_double((isl.left_bound(k) +
             isl.right_bound(k)) / 2);
        zs.push_back(z);
    }
#endif

}

//! isolates real roots at rational boundary in y
void roots_at_y_boundary(const Rational& y_bnd,
       const Rational& threshold, X_coord_vector& roots,
        bool is_sq_free, Polynomial_2 *ppoly) {

    typedef CGAL::Fraction_traits< Rational > F_traits;
    typename F_traits::Numerator_type y_num;
    typename F_traits::Denominator_type y_den;
    typename F_traits::Decompose decompose;

    decompose(y_bnd, y_num, y_den);

    if(ppoly == NULL)
        ppoly = &base_sil_poly;

    typename CGAL::Polynomial_traits_d< Polynomial_2 >::Evaluate_homogeneous
        heval;
//     std::cout << "poly: " << *ppoly << "; y: " << y_bnd << "\n";
    Polynomial_1 eval = CGAL::make_square_free(heval(*ppoly, y_num, y_den));

    univariate_roots(eval, threshold, roots);
}

//! isolate polynomial roots
int univariate_roots(const Polynomial_1& poly, const Rational& threshold,
     X_coord_vector& roots, bool sq_free = true) {

    typedef CGAL::internal::Descartes< Polynomial_1, Rational > Isolator;
    CGAL::internal::Real_roots< X_coordinate_1, Isolator > real_roots;

    roots.reserve(poly.degree()); // at most 'degree' roots

    int n_roots;
    if(sq_free) {
        n_roots = real_roots(poly, std::back_inserter(roots));
    } else {
        n_roots = real_roots(poly, std::back_inserter(roots),
            CGAL::Emptyset_iterator()
            /*std::back_inserter(mults)*/);
    }

    typename AK_1::Algebraic_real_traits::Lower_bound lbd;
    typename AK_1::Algebraic_real_traits::Upper_bound ubd;

    for(typename X_coord_vector::iterator xit = roots.begin();
                xit != roots.end(); xit++) {
        const X_coordinate_1& x = *xit;
        while(ubd(x) - lbd(x) > threshold)
            x.refine();
    }
    return n_roots;
}

//! returns adjacency info of two points given by index coordinates \c i1 and
//! \c i2 , returns \c true if two points belong to features of the same time
bool adjacency_info(const Dcel_feature_data& dcel1,
        const Dcel_feature_data& dcel2, Adjacency_vector& adj) {
        
    if(!is_valid)
        throw "Surface analysis is not valid";

    VCH v1, v2;
    ECH e1, e2;
    FCH f1, f2;

    CGAL::Dcel_feature _1 = dcel1.first, _2 = dcel2.first;

    if(_1 == _2) // oops, surrogate bunnies !!
        return true;

    const void *bunny1 = dcel1.second, *bunny2 = dcel2.second;

    // dispatch fluffy bunnies..
    if(_1 == CGAL::VERTEX)  // bunny1 is a vertex
        v1 = VCH(static_cast< typename VCH::pointer >(bunny1));
    else if(_1 == CGAL::EDGE) // bunny1 is an edge
        e1 = ECH(static_cast< typename ECH::pointer >(bunny1));
    else  // bunny1 seems to be a face
        f1 = FCH(static_cast< typename FCH::pointer >(bunny1));

    CGAL::Adjacencies_3 sfx_adj;

    if(_2 == CGAL::VERTEX) { // bunny2 is a vertex
        v2 = VCH(static_cast< typename VCH::pointer >(bunny2));

/*        if(_1 == CGAL::EDGE)
            std::cerr << "e1-v2: " << e1->curve() << "; " << v2->point()
                << "\n";
        else
            std::cerr << "f1-v2: " << f1.ptr() << "; " <<
                v2->point()
                << "\n";*/
        
        sfx_adj = (_1 == CGAL::EDGE ? base_cad.adjacency(e1, base_surf, v2) :
                                      base_cad.adjacency(f1, base_surf, v2));

    } else if(_2 == CGAL::EDGE) { // bunny2 is an edge
        e2 = ECH(static_cast< typename ECH::pointer >(bunny2));

//         if(_1 == CGAL::FACE)
//             std::cerr << "f1-e2: " << f1.ptr() << "; " << e2->curve()
//                 << "\n";
//         else
//             std::cerr << "v1-e2: " << v1->point() << "; " << e2->curve()
//                 << "\n";

        sfx_adj = (_1 == CGAL::FACE ? base_cad.adjacency(f1, base_surf, e2) :
                                      base_cad.adjacency(v1, base_surf, e2));

    } else { // bunny2 seems to be a face
        f2 = FCH(static_cast< typename FCH::pointer >(bunny2));

//         if(_1 == CGAL::VERTEX)
//             std::cerr << "v1-f2: " << v1->point() << "; " << f2.ptr()
//                 << "\n";
//         else
//             std::cerr << "e1-f2: " << e1->curve() << "; " << f2.ptr()
//                 << "\n";

        sfx_adj = (_1 == CGAL::VERTEX ? base_cad.adjacency(v1, base_surf, f2) :
                                      base_cad.adjacency(e1, base_surf, f2));
    }

    adj.resize(sfx_adj.size());
    std::copy(sfx_adj.begin(), sfx_adj.end(), adj.begin());

    return false;
}

//! \c in and \c normals can point to the same container (in-place mod)
void compute_normals(const Point_vec_3d& in,
        Point_vec_3d& normals, bool is_inversed) {

    typedef CGAL::Polynomial_type_generator<double, 3>::Type Poly_double_3;

    if(!is_valid) // surface is not properly loaded
        return;

    // NOTE NOTE NOTE: no need to compute normals for silhouette points !!!

    Poly_double_3 base_poly_d = CGAL::to_double(base_poly);

    typedef CGAL::Polynomial_traits_d< Poly_double_3 > P_traits;
    typename P_traits::Differentiate der;
    typename P_traits::Substitute subst;

    // compute gradient vector
    Poly_double_3 pdx = der(base_poly_d, 0),
        pdy = der(base_poly_d, 1), pdz = der(base_poly_d, 2);

    double coords[3];
    normals.resize(in.size());
    for(unsigned i = 0; i < in.size(); i++) {

        coords[0] = in[i][0], coords[1] = in[i][1], coords[2] = in[i][2];
        double dx = subst(pdx, coords, coords+3),
            dy = subst(pdy, coords, coords+3),
            dz = subst(pdz, coords, coords+3),
            det = dx*dx + dy*dy + dz*dz;

        // some normals can be ill-posed
        det = (det >= 1e-10 ? 1.0 / std::sqrt(det) : 0.0);
        if(is_inversed)
            det = -det;

        normals[i][0] = dx * det;
        normals[i][1] = dy * det;
        normals[i][2] = dz * det;
    }
}

~Surface_analysis_impl() {
}

private:
    //! private stuff

Polynomial_3 shear_surface(Polynomial_3 p, int i, int j) {

    typename PT_3::Shift shift;
    Polynomial_3 x = shift(Polynomial_3(1), 1, 0),
        y = shift(Polynomial_3(1), 1, 1), z = shift(Polynomial_3(1), 1, 2);

    Polynomial_3 reps[] = { x + j*y, y + i*z, z }; // replacements for x, y, z
    Polynomial_3 res =  typename PT_3::Substitute()(p, reps, reps+3);

    return res;
}

//! what = 0: x-z rotation, what = 1: y-z rotation, what = 2: x-y rotation
static Polynomial_3 affine_transform(Polynomial_3 p, unsigned n,
        unsigned what) {

    int m = n + 1; // affine rotation with the angle close to 0
    int a = 2*m*n, b = m*m - n*n, c = m*m + n*n;
    Rational cosr(a, c), sinr(b, c);

    typename CGAL::Polynomial_traits_d< Poly_rat_3 >::Shift shift;
    Poly_rat_3 x = shift(Poly_rat_3(1), 1, 0),
        y = shift(Poly_rat_3(1), 1, 1), z = shift(Poly_rat_3(1), 1, 2);

    Poly_rat_3 reps[3];
    if(what == 0) {
        reps[0] =  x*cosr - z*sinr, reps[1] = y, reps[2] = x*sinr + z*cosr;
    } else if(what == 1) {
        reps[0] =  x*cosr - y*sinr, reps[1] = x*sinr + y*cosr, reps[2] = z;
    } else {
        reps[0] =  x, reps[1] = y*cosr - z*sinr, reps[2] = y*sinr + z*cosr;
    }

    Poly_rat_3 ratpoly =  typename PT_3::Substitute()(p, reps, reps+3);

    Polynomial_3 res;
    typedef CGAL::Fraction_traits< Poly_rat_3 > FTraits;
    typename FTraits::Denominator_type det(1);
    typename FTraits::Decompose decompose;
    decompose(ratpoly, res, det);

    return res;
}

private:
    //! private data

    //! defining polynomial and same with swapped x & z
    Polynomial_3 base_poly, base_poly_swapped; 
    //! defining surface
    Surface_3 base_surf;
    //! defining CAD
    Restricted_cad_3 base_cad;
    //! surface's silhouette defining polynomial
    Polynomial_2 base_sil_poly;
    //! surface's silhouette 
    CA_2 base_ca_sil;

    //! instance of kernel_2 for passed caching
    const Kernel_2 *kernel_2;
    //! indicates that the surface is valid
    bool is_valid;

}; // Surface_analysis_impl
    
} // namespace CGAL

#endif // CGAL_SURFACE_ANALYSIS_IMPL_H
