// Copyright (c) 2009, 2010 Max-Planck-Institut fuer Informatik (Germany).
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
// $URL: svn+ssh://asm@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Algebraic_kernel_2/include/CGAL/Subdivision_knot.h $
// $Id: Root_box_2.h 54940 2010-03-26 10:33:49Z eric $
//
//
// Author(s): Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_2_SUBDIVISION_KNOT_H
#define CGAL_ALGEBRAIC_KERNEL_2_SUBDIVISION_KNOT_H

/*! \file Subdivision_knot.h
 *  flyweight object representing a point for subdivision in C^2
 *  shares floating-point approximation of rational coordinates through rep
 */

#include <CGAL/config.h>
#include <boost/optional.hpp>

namespace CGAL {

namespace internal {

//! rep class ?? 
template < class Bound_, class BFI_ >
struct Subdivision_knot_rep {

    //! this instance's first template parameter
    typedef Bound_ Bound;
    //! this instance's second template parameter
    typedef BFI_ BFI;    

    //! our lovely bigfloats
    typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat; 
    
    typedef std::complex< Bound > C_bound;

    typedef std::complex< BigFloat > C_bigfloat;

    typedef std::complex< double > C_double;
    
    //!\name Constructors
    //!@{

      //! default constructor
    Subdivision_knot_rep() {
    }

    Subdivision_knot_rep(const Bound& re, const Bound& im) :
        z(C_bound(re, im)) {
    }

    //! non-default constructor 
    Subdivision_knot_rep(const C_bound& z_) :
        z(z_) {                
    }
            
    //!@}
    //!\name props
    //!@{

    C_bound z; // rational point
    // NOTE: these two things are mutually exclusive
    // bigfloat approximation (lo / hi)
    mutable boost::optional< std::pair< C_bigfloat, C_bigfloat > > z_bf; 
    // double approximation (lo / hi)
    mutable boost::optional< std::pair< C_double, C_double > > z_d;     

    //!@}
}; // class Subdivision_knot_rep

template < class Bound_, class BFI_ >
class Subdivision_knot: public CGAL::Handle_with_policy<
        Subdivision_knot_rep < Bound_, BFI_ > > {

public:
    //!\name seemingly typedefs
    //!@{

    //! this instance's first template parameter
    typedef Bound_ Bound;
    //! this instance's second template parameter
    typedef BFI_ BFI;    

    //! rep class
    typedef Subdivision_knot_rep < Bound, BFI > Rep;
    //! base class
    typedef CGAL::Handle_with_policy< Rep > Base;

    //! extract types from the rep class
    typedef typename Rep::BigFloat BigFloat; 
    //! complex value over Rationals
    typedef typename Rep::C_bound C_bound;
    //! complex value over bigfloats
    typedef typename Rep::C_bigfloat C_bigfloat;
    //! complex value over doubles
    typedef typename Rep::C_double C_double;

    typedef Interval_nt< true > Interval_d;

    //! myself
    typedef Subdivision_knot< Bound, BFI > Self;
  
    //!@}
public:
    //!\name Constructors
    //!@{

    //! default constructor
    Subdivision_knot() : 
        Base(Rep()) {
    }
  
    //! constructs a rational point in C: re + im*I
    Subdivision_knot(const Bound& re, const Bound& im) :
        Base(Rep(re, im)) {
    }

    //! constructs a rational point in C 
    explicit Subdivision_knot(const C_bound& z_) :
        Base(Rep(z_)) {
    }
  
    //!@}
    //!\name access point & approximations
    //!@{

    using Base::ptr;
    
    //! rational point
    const C_bound& z() const {
        return ptr()->z;
    }
  
    //! double approximation (lower bound)
    const C_double& zd_low() const {
        if(!ptr()->z_d) {
            _approximate_double();
        }
        return ptr()->z_d->first;
    }

    //! double approximation (upper bound)
    const C_double& zd_high() const {
        if(!ptr()->z_d) {
            _approximate_double();
        }
        return ptr()->z_d->second;
    }

    //! bigfloat approximation (lower bound)
    //! \c precision_change forces the bigfloat approximation being recomputed
    const C_bigfloat& zbf_low(bool precision_change = false) const {
        if(!ptr()->z_bf || precision_change) {
            _approximate_bigfloat();
        }
        return ptr()->z_bf->first;
    }

    //! bigfloat approximation (upper bound)    
    //! \c precision_change - the same as above
    const C_bigfloat& zbf_high(bool precision_change = false) const {
        if(!ptr()->z_bf || precision_change) {
            _approximate_bigfloat();
        }
        return ptr()->z_bf->second;
    }

    //!@}
protected:
    //!\name protected stuff: this is the property of Umbrella Corporation 
    //! (unauthorized access will be prosecuted with the full force of the law)
    //!@{

    void _approximate_double() const {
        
        Interval_d z_re = CGAL::to_interval(ptr()->z.real()),
                   z_im = CGAL::to_interval(ptr()->z.imag());
        // make pair (lower / upper)
        ptr()->z_d = std::make_pair(
              C_double(CGAL::lower(z_re), CGAL::lower(z_im)),
              C_double(CGAL::upper(z_re), CGAL::upper(z_im)));
    }

    void _approximate_bigfloat() const {
        BFI z_re = CGAL::convert_to_bfi(ptr()->z.real()),
            z_im = CGAL::convert_to_bfi(ptr()->z.imag());
        // make pair (lower / upper)
        ptr()->z_bf = std::make_pair(
              C_bigfloat(CGAL::lower(z_re), CGAL::lower(z_im)),
              C_bigfloat(CGAL::upper(z_re), CGAL::upper(z_im)));
    }

public:

    void write(std::ostream& os) const {
        os << "Knot: " << z() << "; bf (" <<  zbf_low() << ")\n";
    }
    //!@}
};

//! output
template < class Bound_, class BFI_ >
std::ostream& operator<< (std::ostream& os,
        const Subdivision_knot< Bound_, BFI_ >& knot) {
    
    knot.write(os);
    return os;
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_SUBDIVISION_KNOT_H
// EOF

