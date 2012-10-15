// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#include <CGAL/basic.h>
#include <cassert>
#include <vector>

#include <CGAL/ipower.h>
/*#include <NiX/basic.h>
#include <NiX/number_type_utils.h>*/

#ifndef CGAL_TEST_REAL_ROOT_ISOLATOR_H
#define CGAL_TEST_REAL_ROOT_ISOLATOR_H

CGAL_BEGIN_NAMESPACE
  
namespace internal {
	
template <class RealRootIsolator, class Polynomial_>
int check_intervals_real_root_isolator(
    const Polynomial_& poly) {
  
  typedef RealRootIsolator Isolator;
  typedef typename Isolator::Bound Bound;
#ifdef CGAL_ISOLATOR_USES_INPUT_RANGE
  Isolator isolator(poly.begin(), poly.end());
#else
  Isolator isolator(poly);
#endif
  int n = isolator.number_of_real_roots();
  for(int i=0; i<n; i++) {
    Bound left  = isolator.left_bound(i);
    Bound right = isolator.right_bound(i);
    
    if (!isolator.is_exact_root(i)) {
      assert(left < right);
      //std::cout << " left = " << left << std::endl;
      //std::cout << " right = " << right << std::endl;
      //std::cout << " poly = " << poly << std::endl;
      assert(poly.sign_at(left) * poly.sign_at(right) == CGAL::NEGATIVE);
    } else {
      assert(left == right);
      assert(poly.sign_at(left) == CGAL::ZERO);
    }
  }
  return n;
};


// Not part of the concept 
/*
  template <class RealRootIsolator, class Polynomial, class Bound>
  int check_intervals_real_root_isolator( const Polynomial& p, 
  const Bound& a,
  const Bound& b) {
  
  typedef RealRootIsolator Isolator;
  Isolator isolator(p,a,b);
  int n = isolator.number_of_real_roots();
  for(int i=0; i<n; i++) {
  Bound left = isolator.left_bound(i);
  Bound right = isolator.right_bound(i);
  
  assert( left < right || isolator.is_exact_root(i));
  if(!isolator.is_exact_root(i)) {
  //std::cout << " left = " << left << std::endl;
  //std::cout << " right = " << right << std::endl;
  //std::cout << " p = " << p << std::endl;
  assert(CGAL::evaluate(p,left) * CGAL::evaluate(p,right) < 0);
  }
  }
  return n;
  };
*/


template <class RealRootIsolator>
void test_real_root_isolator() {
    typedef RealRootIsolator Isolator;
#ifdef CGAL_ISOLATOR_USES_INPUT_RANGE
    typedef typename Isolator::Coefficient NT;
    typedef CGAL::Polynomial< NT > Polynomial;
#else
    typedef typename Isolator::Polynomial    Polynomial;
    typedef typename Polynomial::NT          NT;
#endif
    typedef typename Isolator::Bound Bound;
     
    // just some Polynomials (not all are used)
    Polynomial p_00(NT(0));                   // zero polynomial
    Polynomial p_01(NT(1));                   // constant polynomial
    Polynomial p_1(NT(-1),NT(1));       //(x-1)
    Polynomial p_2(NT(-2),NT(1));       //(x-2)
    Polynomial p_3(NT(-3),NT(1));       //(x-3)
    Polynomial p_4(NT(-4),NT(1));       //(x-4)
    Polynomial p_12=p_1*p_2;    //(x-1)(x-2)
    Polynomial p_123=p_1*p_2*p_3;    //(x-1)(x-2)(x-3)
    Polynomial p_s2(NT(-2),NT(0),NT(1)); //(x^2-2)
    Polynomial p_s5(-NT(5),NT(0),NT(1)); //(x^2-5)
    Polynomial p_s10(-NT(10),NT(0),NT(1)); //(x^2-10)
    Polynomial p_s30(-NT(30),NT(0),NT(1)); //(x^2-30)
    Polynomial p_s2510= p_s2*p_s5*p_s10;
    Polynomial p_s530= p_s5*p_s30;
    
    // special cases: 
    /*
    {
        // from zero polynomial
#ifdef CGAL_ISOLATOR_USES_INPUT_RANGE
        Isolator isolator(p_00.begin(), p_00.end());
#else
        Isolator isolator(p_00);
#endif
        assert(isolator.number_of_real_roots() == -1);
#ifndef CGAL_ISOLATOR_USES_INPUT_RANGE
        assert(isolator.polynomial() == p_00); 
#endif
    }
    */
    
    {
        // from constant polynomial = 1 
        Polynomial poly(p_01);
#ifdef CGAL_ISOLATOR_USES_INPUT_RANGE
        Isolator isolator(poly.begin(), poly.end());
#else
        Isolator isolator(poly);
#endif
        assert(isolator.number_of_real_roots() == 0);
#ifndef CGAL_ISOLATOR_USES_INPUT_RANGE
        assert(isolator.polynomial() == p_01); 
#endif
    }{
        // copy constructor
#ifdef CGAL_ISOLATOR_USES_INPUT_RANGE
        Isolator isolator_1(p_123.begin(),p_123.end());
#else
        Isolator isolator_1(p_123);
#endif
        Isolator isolator_2(isolator_1);
        assert(isolator_1.number_of_real_roots() == 
                isolator_2.number_of_real_roots());
#ifndef CGAL_ISOLATOR_USES_INPUT_RANGE
        assert(isolator_1.polynomial() == isolator_2.polynomial()); 
#endif
    }{
        // assign
#ifdef CGAL_ISOLATOR_USES_INPUT_RANGE
        Isolator isolator_1(p_123.begin(),p_123.end());
#else
        Isolator isolator_1(p_123);
#endif
        Isolator isolator_2 = isolator_1;
        assert(isolator_1.number_of_real_roots() == 
                isolator_2.number_of_real_roots());
#ifndef CGAL_ISOLATOR_USES_INPUT_RANGE
        assert(isolator_1.polynomial() == isolator_2.polynomial()); 
#endif
    }
    {  
        int n = internal::check_intervals_real_root_isolator<Isolator>(p_123);
        assert( n == 3);
       
    }{  
        int n = internal::check_intervals_real_root_isolator<Isolator>(p_1);
        assert( n == 1);
       
    }{  
        int n = internal::check_intervals_real_root_isolator<Isolator>(p_123);
        assert( n == 3);
        
    }{  
        int n = internal::check_intervals_real_root_isolator<Isolator>(p_s2510);
        assert( n == 6);
       
    }{  // (x^2-2)*(x^2-3)
        std::vector<NT> vp(5);
        vp[0] = NT(6);
        vp[1] = NT(0);
        vp[2] = NT(-5);
        vp[3] = NT(0);
        vp[4] = NT(1);  
        Polynomial p(vp.begin(), vp.end()); 
        int n = internal::check_intervals_real_root_isolator<Isolator>(p);
        assert( n == 4);
       
    }
    {  // (x^2-2)*(x^2+2)*(x-1)
        std::vector<NT> vp(6);
        vp[0] = NT(4);
        vp[1] = NT(-4);
        vp[2] = NT(0);
        vp[3] = NT(0);
        vp[4] = NT(-1);  
        vp[5] = NT(1);  
        Polynomial p(vp.begin(), vp.end());
        int n = internal::check_intervals_real_root_isolator<Isolator>(p);
        assert( n == 3);
    }{       
        // std::cout << "Wilkinson Polynomial\n";
        int number_of_roots = 20;
        Polynomial p(1);
        for(int i=1; i<=number_of_roots ; i++) {
            p*=Polynomial(NT(-i),NT(1));
        }
#ifdef CGAL_ISOLATOR_USES_INPUT_RANGE
        Isolator isolator(p.begin(), p.end());
#else
        Isolator isolator(p);
#endif
        int n = internal::check_intervals_real_root_isolator<Isolator>(p);
        assert( n == number_of_roots);             
    }{
        //std::cout << "Kameny 3\n";
        // from http://www-sop.inria.fr/saga/POL/BASE/1.unipol
 
        NT c = CGAL::ipower(NT(10),12);
        Polynomial p(NT(-3),NT(0),c);
        p = p*p;   // (c^2x^2-3)^2
        Polynomial q(NT(0),NT(1));
        q = q*q; // x^2
        q = q*q; // x^4
        q = q*q; // x^8
        q = q*Polynomial(NT(0),c);//c^2x^9
        p = p+q;
        
        int n = internal::check_intervals_real_root_isolator<Isolator>(p);
        assert(n == 3);
    }{
        //std::cout << "Kameny 4\n";
        // from http://www-sop.inria.fr/saga/POL/BASE/1.unipol
    
        NT z(0);
        NT a = CGAL::ipower(NT(10),24); // a = 10^{24}
       
        Polynomial p(z,NT(4),CGAL::ipower(a,2),z,z,2*a,z,z,NT(1)); 
        // x^8+2*10^{24}*x^5+10^{48}*x^2+4*x  
        p = p * Polynomial(z,z,z,z,z,z,NT(1));
        // x^{14}+2*10^{24}*x^{11}+10^{48}*x^8+4*x^7
        p = p + Polynomial(NT(4),z,z,z,-4*a);
        // x^{14}+2*10^{24}*x^{11}+10^{48}*x^8+4*x^7-4*10^{24}*X^4+4
    
#ifdef CGAL_ISOLATOR_USES_INPUT_RANGE
        Isolator isolator(p.begin(), p.end());
#else
        Isolator isolator(p);
#endif
        int n = internal::check_intervals_real_root_isolator<Isolator>(p);
        assert( n == 4);
    }{
        //std::cout << "Polynomial with large and small clustered roots\n";
        // from http://www-sop.inria.fr/saga/POL/BASE/1.unipol
        // there seems to be some error or misunderstanding
        
        NT z(0);
        NT a = CGAL::ipower(NT(10),20); // a = 10^{20}
        
        Polynomial p(z,z,z,z,z,z,z,z,NT(1)); //x^8
        p = p*Polynomial(z,z,z,z,NT(1)); // x^{12}
        Polynomial r(NT(-1),a); // ax-1
        r = r*r;
        r = r*r; // (ax-1)^4
        p = p-r; // x^{12} - (ax-1)^4
        
        int n = internal::check_intervals_real_root_isolator<Isolator>(p);
        assert( n == 4);  
    }
};

} //namespace internal

CGAL_END_NAMESPACE

#endif // CGAL_TEST_REAL_ROOT_ISOLATOR_H
