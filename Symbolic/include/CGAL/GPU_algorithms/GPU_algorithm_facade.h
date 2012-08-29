// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// this file is not part of any library ;-)
//
// ----------------------------------------------------------------------------
//
// Library       : CUDA MP
//
// File          : 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

//!\file GPU_algorithm_facade.h a facade to run algorithms on the GPU

#ifndef CGAL_GPU_ALGORITHM_FACADE_H
#define CGAL_GPU_ALGORITHM_FACADE_H

#include <CGAL/GPU_algorithms/resultant_algorithm.h>
#include <CGAL/GPU_algorithms/gcd_algorithm.h>
#include <CGAL/GPU_algorithms/GPU_algorithm_templates.h>

namespace CGAL {

struct GPU_algorithm_facade {

    template < class NT >
    static CGAL::Polynomial< NT > gcd(const CGAL::Polynomial< NT >& f_,
                const CGAL::Polynomial< NT >& g_, bool *pfailed = 0) {

        typedef CGAL::Polynomial< NT > Poly;
    
        if(CGAL::is_zero(f_)) {
            if(CGAL::is_zero(g_))
                return Poly(NT(1));
            else 
                return g_;
        }
        if(CGAL::is_zero(g_) || f_ == g_)
            return f_;

        Poly ggcd;
        // filter out constant polys
        if(f_.degree() != 0 && g_.degree() != 0) {

        Poly f = f_ / f_.content(), g = g_ / g_.content();
//         Poly f = CGAL::canonicalize(f_), g = CGAL::canonicalize(g_);
        internal::GPU_gcd& obj = internal::GPU_gcd::instance();

        internal::MPZ_vector_1 fv, gv, rv;
        internal::construct_mpz_vector_1(f, fv);
        internal::construct_mpz_vector_1(g, gv);

        unsigned deg_f = f.degree(), deg_g = g.degree(), bits;
        if(deg_f < deg_g) {
            std::swap(fv, gv);
            std::swap(deg_f, deg_g);
        }

        internal::compute_gcd_bitlength(f, g, bits);

        if(pfailed != 0)
            *pfailed = false;

        if(obj.internal_compute(fv, gv, rv, bits)) {
            ggcd = internal::construct_polynomial_from_mpz< NT >(rv);
//             std::cout << "GGCD succeeded\n";
            internal::dispose_mpz_vector(rv);
            
        } else  {
            std::cout <<  "GGCD failed\n";
            if(pfailed != 0)
                *pfailed = true;
            return Poly(NT(0));
        }
        ggcd /= ggcd.content();

        } else {
            ggcd = Poly(NT(1));
        }

        NT gcdcont = CGAL::gcd(f_.content(), g_.content());
        ggcd *= gcdcont;
        return ggcd;
    }

    template < class NT >
    static CGAL::Polynomial< NT > resultant(
        const CGAL::Polynomial< CGAL::Polynomial< NT > >& f,
        const CGAL::Polynomial< CGAL::Polynomial< NT > >& g,
                bool *pfailed = 0) {

         if(pfailed != 0)
           *pfailed = false;

        if(CGAL::is_zero(f) || CGAL::is_zero(g)) {
            return CGAL::Polynomial< NT >(NT(0));
        }

        if(f.degree() == 0) {
            return (CGAL::ipower(CGAL::canonicalize(f.lcoeff()), g.degree()));
        } else if(g.degree() == 0) {
            return (CGAL::ipower(CGAL::canonicalize(g.lcoeff()), f.degree()));
        }

        unsigned low_deg, high_deg, bits;

        internal::GPU_resultant& obj = internal::GPU_resultant::instance();
        internal::compute_resultant_bounds(f, g, low_deg, high_deg, bits);

        internal::MPZ_vector_2 fv, gv, *pfv = &fv, *pgv = &gv;
        unsigned deg_x1, deg_x2, deg_f = f.degree(), deg_g = g.degree();

        internal::construct_mpz_vector_2(f, fv, deg_x1);
        internal::construct_mpz_vector_2(g, gv, deg_x2);

        if(deg_f < deg_g) {
            std::swap(pfv, pgv);
            std::swap(deg_x1, deg_x2);
            std::swap(deg_f, deg_g);
        }

        CGAL::Polynomial< NT > res;
        bool failed = true;

        internal::MPZ_vector_1 r;

        if(obj.internal_compute(*pfv, *pgv, r, deg_f, deg_g, deg_x1, deg_x2,
                low_deg, high_deg, bits)) {

            res = internal::construct_polynomial_from_mpz< NT >(r);
            internal::dispose_mpz_vector(r);
//             std::cout <<  "GRES succeeded\n";
            failed = false;
        }

         // if the algorithm failed here: decompose one of the polynomials
         // as g * y^k, then the resultant is res(f, g/y^k) * f[0]^k
        if(failed) { // trying to remove trailing zero coeffs
            CGAL::Polynomial< CGAL::Polynomial< NT > >
                poly = f, reduced = g;

            if(CGAL::is_zero(poly[0])) {
                std::swap(poly, reduced);
            }

            // failed due to unkwnown reason or polynomials are not coprime
            if(CGAL::is_zero(poly[0]) || !CGAL::is_zero(reduced[0])) {
                std::cerr <<  "GRES failed; non-coprime ?\n";
                if(pfailed != 0)
                    *pfailed = true;

                return CGAL::Polynomial< NT >(0);
            }

            std::cout << "Trying factoring the resultant..\n";
            int nzeros = 0;
            while(1) {
               if(!CGAL::is_zero(reduced[0]))
                    break;
                nzeros++;
                reduced.divide_by_x();
            }
            // compose the resultant from pieces
            res = resultant(poly, reduced, &failed);
            if(failed) {
                if(pfailed != 0)
                    *pfailed = true;
                return CGAL::Polynomial< NT >(0);
            }
            res *= CGAL::ipower(poly[0], (int)nzeros);
        }
        if(pfailed != 0)
            *pfailed = false;
        return CGAL::canonicalize(res);
    }

}; // GPU_algorithm_facade

} // namespace CGAL

#endif // CGAL_GPU_ALGORITHM_FACADE_H
