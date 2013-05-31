// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRcomANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : QdX
// File          : demos/xtriangulate/surf_analysis.C
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.6 $
// Revision_date : $Date: 2009-09-19 01:39:16 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

// #define NDEBUG 1

#include "include/surface_analysis_interface.h"
#include "include/Surface_analysis_impl.h"

#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

//! collects whatever expensive includes we need for surf analysis

#define PARSER_FLOAT_APPROX_BITS 32 // # of bits to approximate floating-point
#define PARSER_MAX_POLY_DEGREE 60 // maximal allowed polynomial degree

// extern Kernel_2 kernel_2_inst;

static CGAL::Surface_analysis_impl< Kernel_2 > surf_engine(0);
    //&kernel_2_inst);

template < class Poly_d_ >
struct Custom_parser_policy :
        public CGAL::Mixed_floating_point_parser_policy< Poly_d_ > {

    //! template argument type
    typedef Poly_d_ Poly_d;
    //! base class
    typedef CGAL::Mixed_floating_point_parser_policy< Poly_d > Base;
    //! type of polynomial coefficient
    typedef typename Base::Coeff Coeff;

    typedef typename CGAL::Get_arithmetic_kernel< Coeff >::
             Arithmetic_kernel AK;
    //! integer number type
    typedef typename AK::Integer Integer;
    //! rational number type
    typedef typename AK::Rational Rational;
    //! BFI type
    typedef typename AK::Bigfloat_interval BFI;
    //! BigFloat type
    typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;
    //! input coefficient types
    typedef typename Base::CoeffType CoeffType;

    virtual Coeff read_coeff_proxy(std::istream& is,
          CoeffType type) const {

        if(type == Base::COEFF_RATIONAL) {

            Integer num, denom;
            is >> CGAL::iformat(num); // read numerator
            is.get(); // skip '/'
            is >> CGAL::iformat(denom);
//             std::cout << "rational: " << num << "/" << denom << "\n";
            if(CGAL::is_zero(denom))
                throw CGAL::internal::Parser_exception("zero div error!");

            typedef CGAL::Fraction_traits< Rational > FT;
            return typename FT::Compose()(num, denom);

        } else if(type == Base::COEFF_FLOAT) {
            double ld;
            is >> CGAL::iformat(ld);
            BigFloat bf(ld);

            long prec = CGAL::get_precision(BFI());
            CGAL::set_precision(BFI(), prec);
            prec = CGAL::set_precision(BFI(), PARSER_FLOAT_APPROX_BITS);
            BFI bfi = CGAL::convert_to_bfi(bf);
            CGAL::set_precision(BFI(), prec);
            return CGAL::lower(bfi);

        } else
            return Base::read_coeff_proxy(is, type);
    }

    //! checking for degree overflow: can be used in real-time applications
    virtual bool exponent_check(unsigned e) const {
//         std::cout << "exponent_check: " << e << "\n";
        if(e > PARSER_MAX_POLY_DEGREE)
            return false;
        return true;
    }
};

bool XSurface_analysis::surface_valid() {
    return surf_engine.surface_valid();
}

const Kernel_2& XSurface_analysis::kernel_2_inst() {
    return surf_engine.kernel_2_inst();
}
   
Curve_analysis_2 XSurface_analysis::silhouette_curve() {
    return surf_engine.silhouette_curve();
}    

Polynomial_2 XSurface_analysis::poly_at_rational_x(const Rational& x) {
    return surf_engine.poly_at_rational_x(x);
}

unsigned XSurface_analysis::locate_pt_at_sil(const X_coordinate_1& x,
        int arcno, Dcel_feature_vec& vec) {
    return surf_engine.locate_pt_at_sil(x, arcno, vec);
}

unsigned XSurface_analysis::locate_pt_in_face(const X_coordinate_1& x,
        const Rational& y, Dcel_feature_vec& vec) {

    return surf_engine.locate_pt_in_face(x, y, vec);
}

void XSurface_analysis::approximate_zs_at_sil(const X_coordinate_1& x,
    int arcno, const Dcel_feature_data& dcel, const double& threshold,
                              Double_vector& zs) {

    surf_engine.approximate_zs_at_sil(x, arcno, dcel, threshold, zs);
}    

void XSurface_analysis::approximate_zs_in_face(const Polynomial_2& poly_at_x,
        const Rational& y, const Rational& threshold,
        const Rational& z_below, const Rational& z_above,
         Double_vector& zs, unsigned& n_skips_below,
          unsigned& n_skips_above, bool use_z_bounds) {

    surf_engine.approximate_zs_in_face(poly_at_x, y, threshold,
        z_below, z_above, zs, n_skips_below, n_skips_above, use_z_bounds);
}

int XSurface_analysis::univariate_roots(const Polynomial_1& poly, const
        Rational& threshold, X_coord_vector& zs, bool sq_free) {
    
    return surf_engine.univariate_roots(poly, threshold, zs, sq_free);
}
    
#if 0
    void exact_zs_in_face(const Xy_coordinate_2& xy,
        unsigned idcel, const Rational& threshold, const Rational& z_below,
        const Rational& z_above, Double_vector& zs) const;
#endif

void XSurface_analysis::roots_at_y_boundary(const Rational& y_bnd,
       const Rational& threshold, X_coord_vector& zs,
        bool is_sq_free, Polynomial_2 *ppoly) {
    
    surf_engine.roots_at_y_boundary(y_bnd, threshold, zs,
           is_sq_free, ppoly);
}

bool XSurface_analysis::adjacency_info(const Dcel_feature_data& dcel1,
             const Dcel_feature_data& dcel2, Adjacency_vector& adj) {
    
    return surf_engine.adjacency_info(dcel1, dcel2, adj);
}
    
void XSurface_analysis::compute_normals(const Point_vec_3d& in,
        Point_vec_3d& normals, bool is_inversed) {

    surf_engine.compute_normals(in, normals, is_inversed);
}    

Curve_analysis_2 XSurface_analysis::z_curve( Rational& z_below,
           Rational& z_above) {
    
    return surf_engine.z_curve(z_below, z_above);
}

void XSurface_analysis::compute_y_range(double& y_min, double& y_max,
         bool& is_set) {

    surf_engine.compute_y_range(y_min, y_max, is_set);
}
   
void XSurface_analysis::dump_arr() {
    surf_engine.dump_arr();
}

//! loads and analyses the surface
bool XSurface_analysis::load_surface(const char *filename) {

    typedef CGAL::Polynomial_type_generator< Rational, 3 >::Type
            Poly_rat_3;
    
    Polynomial_3 in;

    CGAL::Timer tm_load;
    tm_load.start();
    
    std::vector< Polynomial_3 > polys;
    if(!CGAL::read_file< AT >(filename, std::back_inserter(polys))) {
         std::cerr << "Could not read file: " << filename << std::endl;
         return false;
    }
    if(polys.size() != 0 && !CGAL::is_zero(polys[0])) {
        in = polys[0];
    } else {
        std::cerr << "Read failed, trying other format..\n";
    
        std::ifstream file(filename);
        std::string s;
        std::getline(file, s);

        CGAL::Polynomial_parser_d< Poly_rat_3,
           Custom_parser_policy< Poly_rat_3 > > parser;
        Poly_rat_3 ratpoly;

        std::cout << "str: " << s << "\n";

        if(!parser(s, ratpoly)) {
            std::cerr << "Poly parser failed..\n";
            return false;
        }
/// check with helix example
        std::cout << "ratpoly: " << ratpoly << "\n";
        typedef CGAL::Fraction_traits< Poly_rat_3 > FTraits;
        FTraits::Denominator_type det(1);
        FTraits::Decompose decompose;
        decompose(ratpoly, in, det);
    }

    CGAL::set_pretty_mode(std::cout);
    CGAL::set_pretty_mode(std::cerr);

    std::cout << "original surf = " << in << "\n";

    surf_engine.load_surface(in);

    tm_load.stop();
    std::cout << "\nLoad surface time elaplsed: " << tm_load.time () << "\n";
    
    return true;
}
