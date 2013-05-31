// Copyright (c) 2010, 2011, 2012 Max-Planck-Institut fuer Informatik (Germany).
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
// $URL$
// $Id$
//
//
// Author(s): Eric Berberich    <eric.berberich@cgal.org>

#ifndef CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_COFACTOR_ARCAVOID_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_COFACTOR_ARCAVOID_TRAITS_H

/*! \file
 * Traits class implementing cofactors root certification with aracavoid active intervals
 */
#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_2/Bi_solve_2_flags.h>

#include <CGAL/Algebraic_kernel_2/Certifier_cofactor_traits.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel_at_alpha.h>

#include <CGAL/Algebraic_kernel_1/Arcavoid.h>

namespace CGAL {

namespace internal {

template < class AlgebraicKernel_d_1 >
class Certifier_cofactor_arcavoid_traits : public
            Certifier_cofactor_traits < AlgebraicKernel_d_1 > {

public:

    //! type of AK_1
    typedef AlgebraicKernel_d_1 Algebraic_kernel_d_1;

    //! base class
    typedef Certifier_cofactor_traits< Algebraic_kernel_d_1 > Base;

    //! type of solution
    typedef typename Base::Algebraic_real_2 Algebraic_real_2;

    //! type of univariate polynomial
    typedef typename Base::Polynomial_1 Polynomial_1;

    //! type of bivariate polynomial
    typedef typename Base::Polynomial_2 Polynomial_2;

    //! type of algebraic real
    typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;

    //! type of Multiplicity
    typedef typename Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;

    //! type of Bound
    typedef typename Base::Bound Bound;

    //! arithmetic kernel
    typedef typename CGAL::Get_arithmetic_kernel< Bound >::Arithmetic_kernel AK;

    //! bigfloat interval
    typedef typename AK::Bigfloat_interval BFI;

    //! our lovely bigfloats
    typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;

    //! a root box
    typedef typename Base::Root_box Root_box;


    //! range analysis with BFI using real arithmetic
    typedef typename Base::IA_real_bfi IA_real_bfi;

    //! bitstream coefficient kernel
    typedef CGAL::internal::Bitstream_coefficient_kernel_at_alpha< Algebraic_kernel_d_1 > BCK;

    //! Active intervals internal
    typedef CGAL::internal::Arcavoid_list< BCK> Arcavoid_list;

    //! Active intervals
    typedef Arcavoid_list Active_intervals_set;

    //! Active intervals iterator
    typedef typename Arcavoid_list::Cluster_iterator Active_intervals_set_iterator;

    //! const versions
    typedef typename Arcavoid_list::Cluster_const_iterator Active_intervals_set_const_iterator;

public:
    //!\name Constuctors
    //!@{

    //! default
    Certifier_cofactor_arcavoid_traits() :
        Base() {
    }

    //! standard
    Certifier_cofactor_arcavoid_traits
      (Algebraic_kernel_d_1 *ak,
       const Polynomial_2& f, const Polynomial_2& g,
       const Polynomial_2& ft, const Polynomial_2& gt,
       const Polynomial_1& res_x, const Polynomial_1& res_y) :
        Base(ak, f, g, ft, gt, res_x, res_y) {
    }

    //!@}
public:

    //!\name Active intervals
    //!@{

    // TODO cache
    template < class InputIterator >
    Active_intervals_set active_intervals_set_at(const Algebraic_real_1& x,
                                                 InputIterator first, InputIterator beyond,
                                                 long upper_bound_log2_abs = -1) { // -1 indicates local bound

      BCK bck(this->kernel(), x);

      return Arcavoid_list(bck, first, beyond);

    }

    inline
    Bound lower(const Active_intervals_set& ais, Active_intervals_set_const_iterator ait) const {
      //std::cout << "laitc: "<< ait->center().real() << std::endl;
      //std::cout << "laitr: "<< ait->radius() << std::endl;
      return ait->center().real() - ait->radius();
    }

    inline
    Bound upper(const Active_intervals_set& ais, Active_intervals_set_const_iterator ait) const {
      //std::cout << "uaitc: "<< ait->center().real() << std::endl;
      //std::cout << "uaitr: "<< ait->radius() << std::endl;
      return  ait->center().real() + ait->radius();
    }

    inline
    bool is_real(const Active_intervals_set& ais, Active_intervals_set_iterator ait) const {
      return ait->touch_real();
    }

    inline
    std::pair< Active_intervals_set_iterator, Active_intervals_set_iterator >
      refine(Active_intervals_set& ais, Active_intervals_set_iterator ait) {

      return ais.real_subdivide_cluster(ait);

    }

    inline
    Multiplicity_type upper_bound_on_multiplicity(const Active_intervals_set& ais,
                                                  Active_intervals_set_iterator ait) {
      return ait->multiplicity();
    }



    //!@}
protected:

}; // Certifier_cofactor_arcavoid_traits

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_COFACTOR_ARCAVOID_TRAITS_H
// EOF
