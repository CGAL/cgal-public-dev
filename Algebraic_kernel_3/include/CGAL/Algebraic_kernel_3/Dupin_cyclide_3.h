// Copyright (c) 2007-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $URL$
// $Id$
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_3_DUPIN_CYCLIDE_3_H
#define CGAL_ALGEBRAIC_KERNEL_3_DUPIN_CYCLIDE_3_H 1

/*!\file include/CGAL/Algebraic_kernel_3/Dupin_cyclide_3.h
 * \brief Class that defines a ring Dupin_cyclide in 3d.
 */

#include <CGAL/config.h>

#include <boost/array.hpp>

#include <CGAL/Algebraic_kernel_3/Algebraic_surface_3.h>

namespace CGAL {

// pre-declaration
template < class Coefficient_, class Unify_, class Rep_ >
class Dupin_cyclide_3;

namespace internal {

template < class Coefficient_, class Unify_ >
class Dupin_cyclide_3_rep : public 
internal::Algebraic_surface_3_rep< Coefficient_, Unify_> {

public:
  
  //! this instance' first template parameter
  typedef Coefficient_ Coefficient;
  
  //! this instance' second template parameter
  typedef Unify_ Unify;
  
  //! The class itself
  typedef Dupin_cyclide_3_rep< Coefficient, Unify > Self;
  
  //! the parent class
  typedef Algebraic_surface_3_rep< Coefficient, Unify > Base;
  
  //! type of univariate polynomial
  typedef typename Base::Polynomial_1 Polynomial_1;
  
  //! type of bivariate polynomial
  typedef typename Base::Polynomial_2 Polynomial_2;
  
  //! type of trivariate polynomial
  typedef typename Base::Polynomial_3 Polynomial_3;
  
  //! type of polynomial with 4 variables
  typedef typename 
  CGAL::Polynomial_type_generator<Coefficient,4>::Type Polynomial_4;
  
  //! Container for trivariate polynomials
  typedef std::vector< Polynomial_3 > Polynomial_3_container;
  
  //! point for pole
  typedef boost::array< Coefficient, 4 > Point_coeff_4;
  
  //! The linear algebra class
  typedef std::vector< Coefficient > Vector;
  typedef std::vector< Vector > Matrix;
  
  //!\name Constructors
  //!@{
  
  /*!\brief
   * Constructs rep from \c p
   */
  Dupin_cyclide_3_rep(const Polynomial_3& p) :
    Base(p) {
  }
  
  virtual ~Dupin_cyclide_3_rep() {
  }
  
  //!@}
  
private:
  
  void compute_parametrization() {
    
    Coefficient a  =  this->_m_radius_a;
    Coefficient b  =  this->_m_radius_b;
    // FINDING c IS POORLY IMPLEMENTED!!!
    Coefficient square_diff = a*a - b*b;
    Coefficient c = 0;
    while(c*c < square_diff) {
      c = c + Coefficient(1);
    }
    CGAL_assertion(a*a-b*b == c*c);
    Coefficient mu =  this->_m_mu;
    
    Polynomial_1 opus = 
      Polynomial_1(Coefficient(1),Coefficient(0),Coefficient(1));
    Polynomial_1 omus = 
      Polynomial_1(Coefficient(1),Coefficient(0),Coefficient(-1));
    Polynomial_2 opvs = 
      Polynomial_2( Polynomial_1(Coefficient(1)),
                    Polynomial_1(Coefficient(0)),
                    Polynomial_1(Coefficient(1)) );
    Polynomial_2 omvs = 
      Polynomial_2( Polynomial_1(Coefficient(1)),
                    Polynomial_1(Coefficient(0)),
                    Polynomial_1(Coefficient(-1)) );
    
    Polynomial_1 u = Polynomial_1(Coefficient(0),Coefficient(1));
    Polynomial_2 v = Polynomial_2(Polynomial_1(Coefficient(0)),
                                  Polynomial_1(Coefficient(1)));
    
    _m_x = (opvs*(opus*c) - omvs*(omus*a))*Polynomial_1(mu) +
      opvs*(omus*(b*b));
    
    _m_y = (opvs*Polynomial_1(a)-omvs*Polynomial_1(mu))*(u*(Coefficient(2)*b));
    
    _m_z = (v*(Polynomial_1(Coefficient(2)*b)))*(omus*c-opus*mu);
    
    _m_w = opvs*(opus*a) - omvs*(omus*c);
    _vector_transform(_m_x, _m_y, _m_z, _m_w, Polynomial_1());
/////////////////// tube ///////////////////////
    
    _m_tube[0] = (opus*c + omus*a)*mu - opus*(b*b);
    _m_tube[1] = Polynomial_1(0);
    _m_tube[2] = u*(-2*b*(c + mu));
    _m_tube[3] = opus*a + omus*c;
    
    _vector_transform(_m_tube[0], _m_tube[1], _m_tube[2], _m_tube[3], 
                      typename Polynomial_1::NT());
    
//////////////////// outer //////////////////////
    _m_outer[0] = (opus*c + omus*a)*mu + omus*(b*b);
    _m_outer[1] = u*(2*b*(a + mu));
    _m_outer[2] = Polynomial_1(0);
    _m_outer[3] = opus*a + omus*c;
    
    _vector_transform(_m_outer[0], _m_outer[1], _m_outer[2],
                      _m_outer[3], typename Polynomial_1::NT());
    
    _m_pole[0] = -mu*(a - c) - b*b;
    _m_pole[1] = Coefficient(0);
    _m_pole[2] = Coefficient(0);
    _m_pole[3] = a - c;
    
    _vector_transform(_m_pole[0], _m_pole[1], _m_pole[2], _m_pole[3],
                      Coefficient());
#if !NDEBUG
    CGAL::set_ascii_mode(std::cout);
    std::cout << "px:=" << _m_x << ";" << std::endl;
    std::cout << "py:=" << _m_y << ";" << std::endl;
    std::cout << "pz:=" << _m_z << ";" << std::endl;
    std::cout << "pw:=" << _m_w << ";" << std::endl;
    CGAL::set_pretty_mode(std::cout);
    std::cout << "px:=" << _m_x << ";" << std::endl;
    std::cout << "py:=" << _m_y << ";" << std::endl;
    std::cout << "pz:=" << _m_z << ";" << std::endl;
    std::cout << "pw:=" << _m_w << ";" << std::endl;
#endif
  }
  
  // obfuscate code, so enemies don't get what you've written 
  template <class _1, class _2>
  void _vector_transform(_1& _x, _1& _y, _1& _z, _1& _w, _2) {
    
    const Matrix& _3 = _m_base_plane_matrix;
    _1 _x_ = _x*_2(_3[0][0])+
      _y*_2(_3[0][1])+
      _z*_2(_3[0][2])+
      _w*_2(_m_center[0]),
      _y_ = 
      _x*_2(_3[1][0])+
      _y*_2(_3[1][1])+
      _z*_2(_3[1][2])+
      _w*_2(_m_center[1]),
      _z_ = 
      _x*_2(_3[2][0])+
      _y*_2(_3[2][1])+
      _z*_2(_3[2][2])+
      _w*_2(_m_center[2]);
    _x = _x_, _y = _y_, _z = _z_;
  }
  
  // members
  
  
  //! The greater radius of the base ellipse
  Coefficient _m_radius_a;
  
  //! The smaller radius of the base ellipse
  Coefficient _m_radius_b;
  
  //! The parameter defining the tube width
  Coefficient _m_mu;
  
  //! The (orthogonal) matrix defining the torus' base plane
  Matrix _m_base_plane_matrix;
  
  //! The torus center
  Vector _m_center;
  
  //! The parametrization variables
  Polynomial_2 _m_x, _m_y, _m_z, _m_w;
  
  //! the same but more complicated
  Polynomial_1 _m_tube[4];
  Polynomial_1 _m_outer[4];
  
  //! pole
  Point_coeff_4 _m_pole;
  
  // friends
  friend class Dupin_cyclide_3< Coefficient, Unify, Self >;
  
};

} // namespace internal

/*!\brief
 * Class template to represent a Dupin_cyclide in 3d
 */
template <
class Coefficient_,
class Unify_ = CGAL::Handle_policy_no_union,
class Rep_ = CGAL::internal::Dupin_cyclide_3_rep< Coefficient_, Unify_ >
 >
class Dupin_cyclide_3
    : public Algebraic_surface_3< Coefficient_, Unify_, Rep_ > {
    
public:
  
  //! this instance' first template parameter
  
  
  //! this instance' second template parameter
  typedef Coefficient_ Coefficient;
  
  //! this instance' third template parameter
  typedef Unify_ Unify;
  
  //! this instance' fourth template parameter
  typedef Rep_ Rep;
  
  //! The class itself
  typedef Dupin_cyclide_3< Coefficient, Unify, Rep > Self;
  
  //! base type
  typedef Algebraic_surface_3< Coefficient, Unify, Rep > Base;
  
  //! type of univariate polynomial
  typedef typename Rep::Polynomial_1 Polynomial_1;
  
  //! type of bivariate polynomial
  typedef typename Rep::Polynomial_2 Polynomial_2;
  
  //! type of trivariate polynomial
  typedef typename Rep::Polynomial_3 Polynomial_3;

  //! type of fourvariate polynomial
  typedef typename Rep::Polynomial_4 Polynomial_4;
  
  //! type of point for pole
  typedef typename Rep::Point_coeff_4 Point_coeff_4;
  
  //! Less_than type for surfaces
  typedef CGAL::Handle_id_less_than< Self > Surface_less_than;
  
  //! type of surface pair
  typedef std::pair< Self, Self > Two_surfaces;
  
  typedef typename Rep::Matrix Matrix;
  typedef typename Rep::Vector Vector;
  
  //! type of surface pair less than  
  typedef 
  CGAL::Pair_lexicographical_less_than< Self, Self, 
    Surface_less_than, Surface_less_than > 
  Surface_pair_less_than;
  
  //!\name Constructors
  //!@{
  
  /*!\brief default constructor
   *
   * A default-constructed surface supports no operation other than
   * having \c f().degree() return \c -1.
   */
  Dupin_cyclide_3() : 
    Base() {
  };
  
  
  /*!\brief
   * Constructs surface from polynomial \c p.
   */
  Dupin_cyclide_3(const Polynomial_3& p, bool z_regularity_check=false) :
    Base(p,z_regularity_check) {}
  
  /*! \brief Constructs a paramterizable Dupin cyclide from the parameters
   */
  Dupin_cyclide_3(Coefficient radius_a,
                  Coefficient radius_b,
                  Coefficient mu, 
                  Matrix basis, 
                  Vector center) 
    : Base() {
    this->ptr()->_m_radius_a  = radius_a;
    this->ptr()->_m_radius_b  = radius_b;
    this->ptr()->_m_mu        = mu;
    this->ptr()->_m_base_plane_matrix   = basis;
    this->ptr()->_m_center = center;
    this->ptr()->_m_f = Polynomial_3();
    this->ptr()->compute_parametrization();
    // TODO: Implicit equation
  }
  
  //!@}
  
public:
  
  Coefficient radius_a() const {
    return this->ptr()->_m_radius_a;
  }
  
  Coefficient radius_b() const {
    return this->ptr()->_m_radius_b;
  }
  
  Matrix base_plane_matrix() const {
    return this->ptr()->_m_base_plane_matrix;
  }
  
  Vector torus_center() const {
    return this->ptr()->_m_center;
  }
  
  Polynomial_2 x_param() const {
    return this->ptr()->_m_x;
  }
  
  Polynomial_2 y_param() const {
    return this->ptr()->_m_y;
  }
  
  Polynomial_2 z_param() const {
    return this->ptr()->_m_z;
  }
  
  Polynomial_2 w_param() const {
    return this->ptr()->_m_w;
  }
  
  const Polynomial_1 *tube_circle() const {
    return this->ptr()->_m_tube;
  }
  
  const Polynomial_1 *outer_circle() const {
    return this->ptr()->_m_outer;
  }
  
  const Point_coeff_4& pole() const {
    return this->ptr()->_m_pole;
  }
  
  
public:
  
  //!\name Cache
  //!@{
  
protected:
  
  /*!\brief
   * Canonicalizes trivariate polynomial.
   */
  class Canonicalizer {
  public:
    /*!\brief
     * Returns canonicalized version of \c p
     */
    Polynomial_3 operator()(const Polynomial_3& p) {
      return CGAL::canonicalize(p);
    }
  };
  
public:
  
  //! type of surface cache
  typedef CGAL::Cache< Polynomial_3, Self, 
                        CGAL::Creator_1< Polynomial_3,Self >,
                        Canonicalizer
    > Dupin_cyclide_cache;

/*\brief
 * Dupin cyclide cache instance
 */
static Dupin_cyclide_cache& surface_cache() {
  static Dupin_cyclide_cache cache;
  return cache;
}

//!@}
};

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_3_DUPIN_CYCLIDE_3_H
// EOF
