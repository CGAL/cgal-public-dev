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
// $URL: svn+ssh://algerbya@scm.gforge.inria.fr/svn/cgal/branches/features/Algebraic_kernel_d-RS_bivariate-nancy/Algebraic_kernel_d/include/CGAL/Algebraic_kernel_d/Algebraic_real_d_1.h $
// $Id: Algebraic_real_d_1.h 59002 2010-10-04 11:00:27Z lrineau $
// 
// Author: Luis Peñaranda <luis.penaranda@gmx.com>
// Author(s)     :  Yacine Bouzidi <bouzidi.yacine@gmail.com> 
//
// ============================================================================



#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_PURE_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_REAL_PURE_H


#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpfi.h>
#include <CGAL/Handle_for.h>
#include <boost/operators.hpp>


namespace CGAL{

  namespace NewAlg{
    
  template < class Coefficient_ > 
       class New_Algebraic_1;
  
  
  // bool operator<(const New_Algebraic_1&,const New_Algebraic_1&); 
  // bool operator==(const New_Algebraic_1&,const New_Algebraic_1&); 
    
// representation of algebraic numbers
template< class Coefficient_> 
class New_Algebraic_1_rep{
 
  typedef Coefficient_                             Coefficient;
  typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type Polynomial_1;
  
  
 public:
  mutable mpfi_t _mpfi;
  Polynomial_1 _poly;
  int _nr;
  int _mult;
  mpfi_ptr _prev,_next;
  mutable Sign _lefteval;
  
 New_Algebraic_1_rep():
  _poly(NULL),_nr(-1),_mult(-1),
    _prev(NULL),_next(NULL),_lefteval(ZERO){}
  ~New_Algebraic_1_rep(){}
  
 private:
  New_Algebraic_1_rep(const New_Algebraic_1_rep&);
  New_Algebraic_1_rep& operator=(const New_Algebraic_1_rep&);
 };
 
// The class of the algebraic numbers templated by the coefficient type.
 template< class Coefficient_ >
   class New_Algebraic_1 : public Handle_for< New_Algebraic_1_rep< Coefficient_ > >,
  boost::totally_ordered1< New_Algebraic_1< Coefficient_ > >{
 
 public:
  typedef Handle_for<New_Algebraic_1_rep< Coefficient_ > > Base;
  typedef New_Algebraic_1< Coefficient_ >  Self;
  typedef Coefficient_                             Coefficient;
  typedef typename CGAL::Polynomial_type_generator<Coefficient,1>::Type Polynomial_1;
  
 public:

    inline    
      New_Algebraic_1(const Self& i,mpfr_prec_t pl,mpfr_prec_t pr){
      if(pl>pr)
	mpfi_init2(mpfi(),pl);
      else
	mpfi_init2(mpfi(),pr);
        set_mpfi(i.mpfi());
        set_pol(i.pol());
        set_nr(i.nr());
        set_mult(i.mult());
        set_prev(i.prev());
        set_next(i.next());
        set_lefteval(i.lefteval());
}

    inline
      New_Algebraic_1(){
      mpfi_init(mpfi());
    }

inline
  New_Algebraic_1(int i){
        mpz_t temp;
        mpfi_init(mpfi());
        mpfi_set_si(mpfi(),(long int)i);
        mpz_init(temp);
        mpz_set_si(temp,(long int)i);
	mpz_neg(temp, temp);
	// le deuxieme coeff est un int au lieu d'un mpq, a voir si il n'est pas obligatoire de la caster...
        Polynomial_1 *p=new Polynomial_1(temp, 1);
        mpz_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
New_Algebraic_1(unsigned i){
        mpz_t temp;
        mpfi_init(mpfi());
        mpfi_set_ui(mpfi(),i);
        mpz_init(temp);
        mpz_set_ui(temp,(unsigned)i);
	mpz_neg(temp, temp);
        Polynomial_1 *p=new Polynomial_1(temp,1);
        mpz_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
New_Algebraic_1(long i){
        mpz_t temp;
        mpfi_init(mpfi());
        mpfi_set_si(mpfi(),i);
        mpz_init(temp);
        mpz_set_si(temp,i);
	mpz_neg(temp, temp);
        Polynomial_1 *p=new Polynomial_1(temp,1);
        mpz_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
New_Algebraic_1(unsigned long i){
        mpz_t temp;
        mpfi_init(mpfi());
        mpfi_set_ui(mpfi(),i);
        mpz_init(temp);
        mpz_set_ui(temp,i);
	mpz_neg(temp, temp);
        Polynomial_1 *p=new Polynomial_1(temp,1);
        mpz_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
New_Algebraic_1(double d){
        mpq_t temp;
	mpz_t num;
	mpz_t denum;
        mpfi_init(mpfi());
        mpfi_set_d(mpfi(),d);
        mpq_init(temp);
	mpz_init(num);
	mpz_init(denum);
        mpq_set_d(temp,d);
	mpq_get_num(num,temp);
	mpq_get_den(denum,temp);
	mpz_neg(num, num);
        Polynomial_1 *p=new Polynomial_1(num,denum);
        mpq_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
New_Algebraic_1(mpz_srcptr z){
        mpz_t temp;
        mpfi_init(mpfi());
        mpfi_set_z(mpfi(),z);
        mpz_init(temp);
        mpz_set(temp,z);
	mpz_neg(temp, temp);
        Polynomial_1 *p=new Polynomial_1(temp,1);
        mpz_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
New_Algebraic_1(mpq_srcptr q){
        mpfi_init(mpfi());
        mpfi_set_q(mpfi(),q);
	mpq_neg(q, q);
        Polynomial_1 *p=new Polynomial_1(q,1);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
New_Algebraic_1(mpfr_srcptr src){
        Gmpfr r(src);
        mpfi_init2(mpfi(),r.get_precision());
        mpfi_set_fr(mpfi(),r.fr());
        CGAL_assertion(mpfr_equal_p(r.fr(),&mpfi()->left)!=0);
        CGAL_assertion(mpfr_equal_p(r.fr(),&mpfi()->right)!=0);
        Polynomial_1 *rsp=new Polynomial_1(Gmpq(r).mpq());
        set_pol(*rsp);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
New_Algebraic_1(mpfi_srcptr i){
        mpfi_init(mpfi());
        set_mpfi(i);
}

inline
New_Algebraic_1(const Gmpz &z){
        mpq_t temp;
        mpfi_init(mpfi());
        mpfi_set_z(mpfi(),z.mpz());
        mpq_init(temp);
        mpq_set_z(temp,z.mpz());
	mpq_neg(temp, temp);
        Polynomial_1 *p=new Polynomial_1(temp,1);
        mpq_clear(temp);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
New_Algebraic_1(const Gmpq &q){
  mpq_t temp;
  mpfi_init(mpfi());
  mpfi_set_q(mpfi(),q.mpq());
  mpq_init(temp);
  temp = q.mpq;
  mpq_neg(temp, temp);
  Polynomial_1 *p=new Polynomial_1(temp, 1);
        set_pol(*p);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

inline
New_Algebraic_1(const Gmpfr &r){
        mpfi_init2(mpfi(),r.get_precision());
        mpfi_set_fr(mpfi(),r.fr());
        CGAL_assertion(mpfr_equal_p(r.fr(),&mpfi()->left)!=0);
        CGAL_assertion(mpfr_equal_p(r.fr(),&mpfi()->right)!=0);
	r = -r;
        Polynomial_1 *rsp=new Polynomial_1(Gmpq(r).mpq(), 1);
        set_pol(*rsp);
        set_nr(0);
        set_lefteval(CGAL::NEGATIVE);
}

// interesting constructor
inline
New_Algebraic_1(
                         mpfi_srcptr i,
                         const Polynomial_1 &p,
                         int n,
                         int m,
                         mpfi_ptr prevroot,
                         mpfi_ptr nextroot){
        mpfi_init(mpfi());
        set_mpfi_ptr(i);
        set_pol(p);
        set_nr(n);
        set_mult(m);
        set_prev(prevroot);
        set_next(nextroot);
        // we don't evaluate in the sf part of p, since p is sf
        // TODO: add assertion
	// il sagit de trouver le signe de l'evaluation du polynome sur un mpfr
        set_lefteval(p.sign_at(Gmpfr(&(i->left))));
}

// another interesting constructor, where we have already calculated the
// left evaluation
inline
New_Algebraic_1(mpfi_srcptr i,const Polynomial_1 &p,int n,int m,
                         mpfi_ptr prevroot,mpfi_ptr nextroot,CGAL::Sign s){
  mpfi_init(mpfi());
        set_mpfi_ptr(i);
        set_pol(p);
        set_nr(n);
        set_mult(m);
        set_prev(prevroot);
        set_next(nextroot);
        set_lefteval(s);
}
 
  
        /* New_Algebraic_1(const Self&,mpfr_prec_t,mpfr_prec_t); */
        /* New_Algebraic_1(); */
        /* New_Algebraic_1(int); */
        /* New_Algebraic_1(unsigned); */
        /* New_Algebraic_1(long); */
        /* New_Algebraic_1(unsigned long); */
        /* New_Algebraic_1(double); */
        /* New_Algebraic_1(mpz_srcptr); */
        /* New_Algebraic_1(mpq_srcptr); */
        /* New_Algebraic_1(mpfr_srcptr); */
        /* New_Algebraic_1(mpfi_srcptr); */
        /* New_Algebraic_1(const Gmpz&); */
        /* New_Algebraic_1(const Gmpq&); */
        /* New_Algebraic_1(const Gmpfr&); */

        /* // the only interesting constructor */
        /* New_Algebraic_1(mpfi_srcptr, */
        /*             const Polynomial_1&, */
        /*             int, */
        /*             int, */
        /*             mpfi_ptr, */
        /*             mpfi_ptr); */

        /* // the another interesting variant */
        /* New_Algebraic_1(mpfi_srcptr, */
        /*             const Polynomial_1&, */
        /*             int, */
        /*             int, */
        /*             mpfi_ptr, */
        /*             mpfi_ptr, */
        /*             CGAL::Sign); */

 

 inline
mpfi_srcptr mpfi()const{
   return Self::Ptr()->_mpfi;
}

inline
mpfi_ptr mpfi(){
  return Self::ptr()->_mpfi;
}

inline
Gmpfi interval()const{
        return Gmpfi(Self::Ptr()->_mpfi);
}

inline
Gmpfr inf()const{
        return Gmpfr(&(mpfi()->left));
}

inline
Gmpfr sup()const{
        return Gmpfr(&(mpfi()->right));
}

inline
  const Polynomial_1& pol()const{
        return Self::Ptr()->_poly;
}

inline
int nr()const{
        return Self::ptr()->_nr;
}

inline
int mult()const{
        return Self::ptr()->_mult;
}

inline
void set_mpfi_ptr(mpfi_srcptr x){
        // *mpfi()=*x;
        // mpfi_set(mpfi(),x);
        set_mpfi(x);
}

inline
void clear_pol(){
  ~(Self::ptr()->_poly)();
}

inline
void set_pol(const Polynomial_1 &p){
        Self::ptr()->_poly=p;
}

inline
void set_nr(int n){
        Self::ptr()->_nr=n;
}

inline
void set_mult(int m){
        Self::ptr()->_mult=m;
}

inline
void set_prec(mp_prec_t p){
        mpfi_round_prec(mpfi(),p);
}

inline
void set_prev(mpfi_ptr p){
        Self::ptr()->_prev=p;
}

inline
void set_next(mpfi_ptr n){
        Self::ptr()->_next=n;
}

inline
void set_lefteval(Sign s)const{
        Self::Ptr()->_lefteval=s;
}

inline
mp_prec_t get_prec()const{
        return mpfi_get_prec(mpfi());
}

inline
mpfr_srcptr left()const{
        return &(mpfi()->left);
}

inline
mpfr_srcptr right()const{
        return &(mpfi()->right);
}

inline
mpfi_ptr prev()const{
        return Self::ptr()->_prev;
}

inline
mpfi_ptr next()const{
        return Self::ptr()->_next;
}

inline
Sign lefteval()const{
        return Self::ptr()->_lefteval;
}

inline
bool is_consistent()const{
        return(&pol()==NULL?false:true);
}

inline
bool is_point()const{
        return(mpfr_equal_p(&(mpfi()->left),&(mpfi()->right))!=0);
}

inline
bool contains(int n)const{
        return(mpfr_cmp_si(&(mpfi()->left),n)<=0 &&
                        mpfr_cmp_si(&(mpfi()->right),n)>=0);
}

inline
bool contains(mpfr_srcptr n)const{
        return(mpfr_lessequal_p(&(mpfi()->left),n) &&
                        mpfr_greaterequal_p(&(mpfi()->right),n));
}

inline
bool contains(const Gmpz &n)const{
        return(mpfr_cmp_z(&(mpfi()->left),n.mpz())<=0 &&
                        mpfr_cmp_z(&(mpfi()->right),n.mpz())>=0);
}

inline
std::pair<double,double> to_interval()const{
        return std::make_pair(
                        mpfr_get_d(left(),GMP_RNDD),
                        mpfr_get_d(right(),GMP_RNDU));
}

inline
void set_mpfi(mpfi_srcptr x){
        mp_prec_t xp;
        xp=mpfi_get_prec(x);
        if(xp>mpfr_get_prec(left()) || xp>mpfr_get_prec(right()))
                mpfi_set_prec(mpfi(),xp);
        mpfi_set(mpfi(),x);
}

inline
bool overlaps(const Self &a)const{
        if(mpfr_lessequal_p(left(),a.left()))
                return (mpfr_lessequal_p(a.left(),right())!=0);
        else
                return (mpfr_lessequal_p(left(),a.right())!=0);
}

inline
bool is_valid()const{
        return (mpfi_nan_p(mpfi())==0);
}

inline
bool is_finite()const{
        return (mpfi_bounded_p(mpfi())!=0);
}

//template <class _Gcd_policy>
inline
  double to_double()const{
        /*typedef _Gcd_policy     Gcd;
        while(mpfr_get_d(left(),GMP_RNDU)!=mpfr_get_d(right(),GMP_RNDU))
                bisect_n<Gcd>(*this,33);*/
        refine_1(*this,100);
        CGAL_assertion(mpfr_get_d(left(),GMP_RNDD)==
                       mpfr_get_d(right(),GMP_RNDD));
        CGAL_assertion(mpfr_get_d(left(),GMP_RNDU)==
                       mpfr_get_d(right(),GMP_RNDU));
        CGAL_assertion(mpfr_get_d(left(),GMP_RNDN)==
                       mpfr_get_d(right(),GMP_RNDN));
        return mpfr_get_d(right(),GMP_RNDU);
}

inline
  Self sqrt()const{
  mpfi_t s;
  mpfi_init(s);
  mpfi_sqrt(s,mpfi());
  Self ret(s);
  return ret;
}

 
 
 
/* // functions related to the member data */
 
 Self operator+()const;
 Self operator-()const;
 
 /* mpfi_srcptr mpfi()const; */
 /* mpfi_ptr mpfi(); */
 /* Gmpfi interval()const; */
 /* Gmpfr inf()const; */
 /* Gmpfr sup()const; */
 /* const Polynomial_1& pol()const; */
 /* int nr()const; */
 /* int mult()const; */
 /* void set_mpfi(mpfi_srcptr); */
 /* void set_mpfi_ptr(mpfi_srcptr); */
 /* void clear_pol(); */
 /* void set_pol(const Polynomial_1 &); */
 /* void set_nr(int); */
 /* void set_mult(int); */
 /* void set_prec(mp_prec_t); */
 /* void set_prev(mpfi_ptr); */
 /* void set_next(mpfi_ptr); */
 /* void set_lefteval(Sign)const; */
 /* mp_prec_t get_prec()const; */
 /* mpfi_ptr prev()const; */
 /* mpfi_ptr next()const; */
 /* Sign lefteval()const; */
 /* mpfr_srcptr left()const; */
 /* mpfr_srcptr right()const; */
 /* bool is_consistent()const; */
 /* bool is_point()const;   // are endpoints equal? */
 /* bool overlaps(const New_Algebraic_1&)const; */
 /* bool contains(int n)const; */
 /* bool contains(mpfr_srcptr)const; */
 /* bool contains(const Gmpz&)const; */
 
 /* New_Algebraic_1 operator+()const; */
 /* New_Algebraic_1 operator-()const; */
 
 /* bool is_valid()const; */
 /* bool is_finite()const; */
 /* /\*template<class>*\/ double to_double()const; */
 /* std::pair<double,double> to_interval() const; */
 /* New_Algebraic_1 sqrt()const; */
 }; // class New_Algebraic_1
 
 
  } // namespace NewAlg
  
} // namespace CGAL

//#include <CGAL/RS/New_algebraic_1_constructors.h>
#include <CGAL/RS/New_algebraic_1_operators.h>
//#include <CGAL/RS/New_algebraic_1_member.h>       // other member functions
#include <CGAL/RS/New_algebraic_1_comparisons.h>

#include <CGAL/RS/New_algebraic_1_real_embeddable.h>       

#endif  // CGAL_RS_ALGEBRAIC_1_H









