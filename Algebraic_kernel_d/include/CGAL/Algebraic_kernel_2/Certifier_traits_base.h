// Copyright (c) 2010, 2011 Max-Planck-Institut fuer Informatik (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://asm@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Algebraic_kernel_2/include/CGAL/Certifier_cofactor_traits.h $
// $Id: Certifier_cofactor_traits.h 57108 2010-06-25 13:54:53Z eric $
//
//
// Author(s): Pavel Emeliyanenko <asm@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_TRAITS_BASE_H
#define CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_TRAITS_BASE_H 

/*! \file
 * Base class for certifier traits unifying the interface
 */

#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_2/Bi_solve_2_flags.h>

#include <CGAL/Algebraic_kernel_2/Bi_algebraic_real_2.h>

#include <CGAL/Algebraic_kernel_2/Root_box_1.h>
#include <CGAL/Algebraic_kernel_2/Type_conversion_internal_traits.h>

namespace CGAL {

namespace internal {

template < class AlgebraicKernel_d_1 >
class Certifier_traits_base {

public:

   //! type of ak1
   typedef AlgebraicKernel_d_1 Algebraic_kernel_d_1;
   
   //! type of Solution
   typedef Bi_algebraic_real_2< Algebraic_kernel_d_1 > Algebraic_real_2;

   //! type of univariate polynomial
   typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;

   //! type of bivariate polynomial
   typedef typename CGAL::Polynomial< Polynomial_1 > Polynomial_2;  

   //! type of Bound
   typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;

   //! type of Bound
   typedef typename Algebraic_kernel_d_1::Bound Bound;

   //! type of Multiplicity
   typedef typename Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;

    //! arithmetic kernel
   typedef typename CGAL::Get_arithmetic_kernel< Bound >::Arithmetic_kernel AK;

   //! bigfloat interval
   typedef typename AK::Bigfloat_interval BFI;

   //! our lovely bigfloats
   typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;

   //! BFI polynomials
   typedef typename CGAL::Polynomial_type_generator< BFI, 2 >::Type Poly_bfi_2;

   //! a root box
   typedef Root_box_1< Bound > Root_box;

    // no properties to save by default
    struct Candidate_props { 
    };

public:

    //!\name Constuctors
    //!@{

    Certifier_traits_base(Algebraic_kernel_d_1 *ak,
                          const Polynomial_2& f, const Polynomial_2& g,
                          const Polynomial_2& ft, const Polynomial_2& gt,
                          const Polynomial_1& res_x, const Polynomial_1& res_y) :
        _m_ak1(ak),
        _m_f(f), _m_g(g), _m_ft(ft), _m_gt(gt), _m_bfi_prec(0), 
        _m_res_x(res_x), _m_res_y(res_y) {
    }

    //! inits a fiber along which solutions are to be computed (in x or y)
    //! suffix '_x' is left for historical reasons
    void init_fiber(const Root_box& box_x, Multiplicity_type mul_x) {

        _m_box_x = box_x, _m_mul_x = mul_x;
    }

    //!@}

public:
    //! \name getters ?
    //!@{

    //! the algebraic kernel
    Algebraic_kernel_d_1* kernel() {
      return _m_ak1;
    }

    //! swaps x/y coordinates
    void swap_coordinates() {
        std::swap(_m_f, _m_ft);
        std::swap(_m_g, _m_gt);
        std::swap(_m_res_x, _m_res_y);
        // drop precision when swapping coordinates
        _m_bfi_prec = 0;
    }

    //! first polynomial
    const Polynomial_2& f() const {
        return _m_f;
    }

    //! second polynomial
    const Polynomial_2& g() const {
        return _m_g;
    }

    //! transposed f
    const Polynomial_2& ft() const {
        return _m_ft;
    }

    //! transposed g
    const Polynomial_2& gt() const {
        return _m_gt;
    }

    //! res(f, g, y)
    const Polynomial_1& res_in_x() const {
        return _m_res_x;
    }

    //! res(f, g, x)
    const Polynomial_1& res_in_y() const {
        return _m_res_y;
    }

    //! isolating box for x-root
    const Root_box& box_x() const {
        return _m_box_x;
    }

    //!@}

    //!\name main stuff
    //!@{

    //! norm test always succeeds by default
    bool norm_test(const Candidate_props& ps,
            const Bound& x_0, const Bound& y_0,
            const Root_box& box_y, Multiplicity_type mul_y, long& prec) {
        return true;
    }

protected:

    void _approximate_polynomials(long prec) {
        if(prec < _m_bfi_prec)
            return;

        long oldp = CGAL::set_precision(BFI(), prec);
        typename CGAL::Coercion_traits< Polynomial_2, Poly_bfi_2 >::Cast cast;
//TODO: use ia_bfi for caching ??
        _m_f_bfi = cast(_m_f);
        _m_g_bfi = cast(_m_g);

        CGAL::set_precision(BFI(), oldp);
        _m_bfi_prec = prec;
    }

    //!@}

protected:
    //!\name Data members
    //!@{

    //! the algebraic kernel
    Algebraic_kernel_d_1 *_m_ak1;

    //! f polynomial
    Polynomial_2 _m_f;

    //! g polynomial
    Polynomial_2 _m_g;

    //! f polynomial transposed
    Polynomial_2 _m_ft;

    //! g polynomial transposed
    Polynomial_2 _m_gt;

    //! approximation precision of polynomials
    long _m_bfi_prec;

    //! cached polynomials over BFIs
    Poly_bfi_2 _m_f_bfi, _m_g_bfi;

    //! the resultants
    Polynomial_1 _m_res_x, _m_res_y;

    //! multiplicity of a root in x
    Multiplicity_type _m_mul_x;
  
    //! root box for x
    Root_box _m_box_x;


public:

  typedef std::vector< std::pair< Polynomial_1, Multiplicity_type > >
  Squarefree_fact;
  Squarefree_fact *_m_res_in_x_factors_ptr, *_m_res_in_y_factors_ptr;

    //!@}

}; // Certifier_traits_base

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_TRAITS_BASE_H
// EOF

