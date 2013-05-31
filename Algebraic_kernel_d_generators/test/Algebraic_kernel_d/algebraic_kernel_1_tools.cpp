// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de> 
//
// ============================================================================

// Contains:  
// - Test compute_smallest_nonnegative_root
// - Test compare_smallest_nonnegative_roots


#include <CGAL/basic.h>
#include <CGAL/algebraic_kernel_1_tools.h>

template <class Algebraic_kernel_d_1>
void test_algebraic_kernel_1_tools(){
  typedef typename Algebraic_kernel_d_1::Coefficient Coefficient;
  typedef typename Algebraic_kernel_d_1::Bound Bound;
  
  typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Root;
  typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;
 
  typedef CGAL::Polynomial_traits_d<Polynomial_1> PT_1;
  
  Algebraic_kernel_d_1 ak; 
  
  Polynomial_1 x = typename PT_1::Shift()(Polynomial_1(1),1);
  Polynomial_1 p1 = x*x-2;
  Polynomial_1 p2 = x*x-3;

  typename Algebraic_kernel_d_1::Solve_1 solve_1 = ak.solve_1_object();
  
  std::vector<Root> roots1, roots2;
  solve_1(p1,std::back_inserter(roots1),false);
  solve_1(p2,std::back_inserter(roots2),false);
  
  assert(roots1.size() == 2);
  assert(roots2.size() == 2);
  
  Root r1 = roots1[1]; // sqrt(2)
  Root r2 = roots2[1]; // sqrt(3)
  
  assert(r1 == *CGAL::compute_smallest_nonnegative_root(ak,p1));
  assert(r2 == *CGAL::compute_smallest_nonnegative_root(ak,p2));
  
  assert(CGAL::compare_smallest_nonnegative_roots(ak,p1,p2) 
      == CGAL::compare(r1,r2));                 
  
}





int main(){

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  {
    typedef CGAL::GMP_arithmetic_kernel AK;
    typedef AK::Integer Coefficient;
    typedef AK::Rational Bound;
    
    typedef CGAL::Algebraic_kernel_d_1<Coefficient, Bound> AK_1;
    test_algebraic_kernel_1_tools<AK_1>();
  }
#else
  std::cerr << " GMP test skipped " << std::endl;
#endif  
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  {
    typedef CGAL::CORE_arithmetic_kernel AK;
    typedef AK::Integer Coefficient;
    typedef AK::Rational Bound;
    
    typedef CGAL::Algebraic_kernel_d_1<Coefficient, Bound> AK_1;
    test_algebraic_kernel_1_tools<AK_1>();
  }
#else
  std::cerr << " CORE test skipped " << std::endl;
#endif 
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  {
    typedef CGAL::LEDA_arithmetic_kernel AK;
    typedef AK::Integer Coefficient;
    typedef AK::Rational Bound;
    
    typedef CGAL::Algebraic_kernel_d_1<Coefficient, Bound> AK_1;
    test_algebraic_kernel_1_tools<AK_1>();
  }
#else
  std::cerr << " LEDA test skipped " << std::endl;
#endif  
  return 0;
}



