#include <CGAL/Bigfloat_traits.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/GMP_arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>
#include <CGAL/Cartesian_complex.h>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>
#include <iostream>

#include "../examples/formatters.h"

#include <boost/lexical_cast.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

/*
int main_1 () {
  using namespace std;

  typedef CGAL::GMP_arithmetic_kernel AK;
  typedef AK::Bigfloat T1;
  typedef AK::Integer T2;

  typedef CGAL::Cartesian_complex< T1 > Complex1;
  typedef CGAL::Cartesian_complex< T2 > Complex2;

  typedef CGAL::Coercion_traits< T1, T2 >::Type T;
  typedef CGAL::Coercion_traits< Complex1, Complex2 >::Type Complex;

  typedef CGAL::Complex_embeddable_traits< Complex1 > CET1;
  typedef CGAL::Complex_embeddable_traits< Complex2 > CET2;
  typedef CGAL::Complex_embeddable_traits< Complex >  CET;

  typedef CGAL::Algebraic_structure_traits< Complex1 > AST1;
  typedef CGAL::Algebraic_structure_traits< Complex2 > AST2;
  typedef CGAL::Algebraic_structure_traits< Complex >  AST;

  typedef complex< T > StdComplex;

  BOOST_STATIC_ASSERT
    ((boost::is_same< CGAL::Cartesian_complex< T >, Complex >::value));

  Complex1 w (16., 4.);
  Complex2 z (4, 2);
  Complex res;

#define TEST_ASSIGNMENT_OP(op)                                          \
  res = z;                                                              \
  res op w;                                                             \
  cout << "z " << #op << " w\t== " << res << endl
#define TEST_UNARY_OP(op)                                               \
  res = op z;                                                           \
  cout << #op << "z\t== " << res << endl
#define TEST_BINARY_OP(op)                                              \
  res = z op w;                                                         \
  cout << "z " << #op << " w\t== " << (StdComplex)res << endl
  // cout << "z " << #op << " w\t== " << res.real()                        \
  //      << " + " << res.imag() << "*I" << endl
  
  TEST_ASSIGNMENT_OP (=);
  TEST_ASSIGNMENT_OP (+=);
  TEST_ASSIGNMENT_OP (-=);
  TEST_ASSIGNMENT_OP (*=);
  TEST_ASSIGNMENT_OP (/=);

  TEST_UNARY_OP (+);
  TEST_UNARY_OP (-);
  TEST_UNARY_OP (~);

  TEST_BINARY_OP (+);
  TEST_BINARY_OP (-);
  TEST_BINARY_OP (*);
  TEST_BINARY_OP (/);

  if (Complex() == 0) cout << "equality" << endl;
  if (CET::Is_zero() (Complex())) cout << "is_zero" << endl;

  std::string s = "(5/3, -4/-2)";

  try {
    z = boost::lexical_cast< Complex2 > (s);
    cout << "Successfully parsed \"" << s << "\" as " << z << endl;
  } catch (boost::bad_lexical_cast &) {
    cerr << "Could not interpret \"" << s << "\" as a complex number" << endl;
  }
  
  cout.precision (53);
  cout << z << " to double is " << CET2::To_double() (z) << endl;
  cout << z << " to interval is " << CET2::To_interval() (z) << endl;
  // AST::Simplify() (z);
  cout << "simplified: " << z << endl;
  cout << "norm: " << CET2::Norm() (z) << endl;
  cout << "abs: " << z.abs() << endl;
  //  cout << "Unit part of " << w << " is " << AST1::Unit_part() (w) << endl;
  //  cout << "Unit part of " << z << " is " << AST2::Unit_part() (z) << endl;

  return 0;
}
*/

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)
#define PRINTME(x) cout << QUOTEME(x) << ":\t" << x << endl

int main () {
  using namespace std;

  typedef CGAL::GMP_arithmetic_kernel AK;
  typedef AK::Bigfloat BF;
  typedef AK::Bigfloat_interval BFI;

  typedef CGAL::Cartesian_complex< BF > CC;
  typedef CGAL::Cartesian_complex< BFI > CCI;

  typedef CGAL::Coercion_traits< CC, CC >::Type Complex;

  BOOST_STATIC_ASSERT ((boost::is_same< CGAL::Coercion_traits< BF, BFI >::Type, BFI >::value));
  BOOST_STATIC_ASSERT ((boost::is_same< CGAL::Coercion_traits< BFI, BF >::Type, BFI >::value));
  BOOST_STATIC_ASSERT ((boost::is_same< CGAL::Coercion_traits< CC, BF >::Type, CC >::value));
  BOOST_STATIC_ASSERT ((boost::is_same< CGAL::Coercion_traits< CC, BFI >::Type, CCI >::value));
  BOOST_STATIC_ASSERT ((boost::is_same< CGAL::Coercion_traits< CCI, BF >::Type, CCI >::value));
  BOOST_STATIC_ASSERT ((boost::is_same< CGAL::Coercion_traits< CCI, BFI >::Type, CCI >::value));

  CGAL::set_default_precision< BF > (800);
  CGAL::set_precision (BFI(), 1200);

  typedef CGAL::Polynomial_type_generator< BF, 1 >::Type BFPoly;
  typedef CGAL::Polynomial_traits_d< BFPoly > BFPT;
  typedef CGAL::Polynomial_type_generator< BFI, 1 >::Type BFIPoly;
  typedef CGAL::Polynomial_traits_d< BFIPoly > BFIPT;
  typedef CGAL::Polynomial_type_generator< CCI, 1 >::Type CCIPoly;
  typedef CGAL::Polynomial_traits_d< CCIPoly > CCIPT;
  // BOOST_STATIC_ASSERT ((boost::is_same< CGAL::Coercion_traits< CC, BFIPoly >::Type, CCIPoly >::value));
  /*
  BFPoly x = BFPT::Shift() (BFPoly(1), 1);
  BFPoly f = BF (2.) * x*x - x + 18;
  BFIPoly xI = BFIPT::Shift() (BFIPoly(1), 1);
  BFIPoly fI = BFI (2.) * xI*xI - xI + 18;

  PRINTME ((type_name (CGAL::Coercion_traits< CC, BFIPoly >::Type())));

  BFI I = 0;
  */
  BFI A (BF (-1.), BF (1.));
  BFI B (BF (-1.), BF (1.));
  BFI C;

  PRINTME (A);
  PRINTME (B);
  PRINTME (C);

  for (int i = 0; i < 100000; ++i) {
    A *= B;
    C = A+B;
  }

  PRINTME (A);
  PRINTME (B);
  PRINTME (C);

  /*
  CC z (17., 4.);

  typedef CGAL::Coercion_traits< CC, BFI > CT;
  CT::Cast cast;
  BOOST_STATIC_ASSERT ((boost::is_same< CT::Type, CCI >::value));

  PRINTME ((type_name (fI.evaluate (z))));
  CCI lcfI = cast (fI.lcoeff());

  PRINTME (lcfI);

  CCI fI_z = fI.evaluate (z);
  fI_z = fI.evaluate (lcfI);
  // error in cast Gmpfi -> Complex< Gmpfi >
  PRINTME(fI_z);
  PRINTME(fI_z.real().get_precision());
  PRINTME(fI_z.imag().get_precision());

  PRINTME ((type_name (f.evaluate (z))));
  CCI f_z = f.evaluate (z);
  PRINTME(f_z);
  PRINTME(f_z.real().get_precision());
  PRINTME(f_z.imag().get_precision());

  PRINTME(I);
  PRINTME(I.get_precision());
  PRINTME(I.get_precision());
  */
  return 0;
}
