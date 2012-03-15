//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/config.h>

#include <sstream>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_3/Algebraic_surface_3.h>

template < class AK >
void test_Algebraic_surface_3() {

  typedef typename AK::Integer Coefficient;
  
  typedef CGAL::Algebraic_surface_3< Coefficient > Algebraic_surface_3;
  
  typedef typename Algebraic_surface_3::Polynomial_1 Polynomial_1;
  typedef typename Algebraic_surface_3::Polynomial_2 Polynomial_2;
  typedef typename Algebraic_surface_3::Polynomial_3 Polynomial_3;
  
  typename Algebraic_surface_3::Surface_cache cache;
  
  Polynomial_3 t1, t2, t3;
  
  {
    Coefficient ra(10);
    Coefficient rb(5);
    
    Coefficient ra2 = ra*ra;
    Coefficient rb2 = rb*rb;
    
    Polynomial_1 c0(Coefficient(0));
    Polynomial_1 cp1(Coefficient(1));
    Polynomial_1 cm4a2(Coefficient(-4)*ra2);
    Polynomial_1 ca2(ra2);
    Polynomial_1 cb2(rb2);
    
    Polynomial_1 x2(Coefficient(0),Coefficient(0),Coefficient(1));
    
    Polynomial_2 x2y2(x2,c0,cp1);
    
    Polynomial_2 pt1_0 = x2y2 + Polynomial_2(ca2 - cb2);
    Polynomial_2 pt1_1(c0);
    Polynomial_2 pt1_2(cp1);
    Polynomial_2 pt2_0 = x2y2 * Polynomial_2(cm4a2);
    
    Polynomial_3 pt1(pt1_0, pt1_1, pt1_2);
    Polynomial_3 pt2(pt2_0);
    
    t1 = pt1*pt1 + pt2;
  }
  
  //::CGAL::set_pretty_mode(std::cout);
  //std::cout << "t1: " << t1 << std::endl;
  
  Algebraic_surface_3 sf1(cache(t1));
  
  assert(sf1.f() == t1);
  
  {
    Coefficient ra(12);
    Coefficient rb(5);
    
    Coefficient ra2 = ra*ra;
    Coefficient rb2 = rb*rb;
    
    Polynomial_1 c0(Coefficient(0));
    Polynomial_1 cp1(Coefficient(1));
    Polynomial_1 cm4a2(Coefficient(-4)*ra2);
    Polynomial_1 ca2(ra2);
    Polynomial_1 cb2(rb2);
    
    Polynomial_1 x2(Coefficient(0),Coefficient(0),Coefficient(1));
    
    Polynomial_2 x2y2(x2,c0,cp1);
    
    Polynomial_2 pt1_0 = x2y2 + Polynomial_2(ca2 - cb2);
    Polynomial_2 pt1_1(c0);
    Polynomial_2 pt1_2(cp1);
    Polynomial_2 pt2_0 = x2y2 * Polynomial_2(cm4a2);
    
    Polynomial_3 pt1(pt1_0, pt1_1, pt1_2);
    Polynomial_3 pt2(pt2_0);
    
    t2 = pt1*pt1 + pt2;
  }
  
  // cache
  //::CGAL::set_pretty_mode(std::cout);
  //std::cout << "t1: " << t1 << std::endl;
  
  Algebraic_surface_3 sf2(cache(t2));
  
  assert(sf2.f() == t2);
  
  
  Algebraic_surface_3 sf3(cache(t1));
  
  assert(sf3.id() == sf1.id());
  
  // operator<<
  
  std::stringstream strtest;
  std::stringstream strout;
  strtest << sf3.f() << std::endl;
  strout << sf3 << std::endl;
  assert(strout.str() == strtest.str());
  
  // TODO test for resultant/cofactors etc
}

int main() {
#if CGAL_USE_LEDA
    {
        typedef CGAL::LEDA_arithmetic_kernel AK;;
        test_Algebraic_surface_3< AK >();
    }
#endif
    {
        typedef CGAL::CORE_arithmetic_kernel AK;
        test_Algebraic_surface_3< AK >();
    }
    return EXIT_SUCCESS;
}
