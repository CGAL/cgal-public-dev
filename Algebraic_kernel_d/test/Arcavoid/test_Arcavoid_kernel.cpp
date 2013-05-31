#include <iostream>

#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Arcavoid_kernel.h>

int main () {
  using namespace std;
  using namespace CGAL;

#ifdef CGAL_USE_LEDA
  cerr << "Test LEDA" << endl
       << "================" << endl;
  {
    typedef leda_integer I;
    typedef leda_rational Q;
    typedef leda_bigfloat F;
    typedef leda_real R;

    typedef Polynomial_type_generator< I, 1 >::Type IPoly;
    typedef Polynomial_type_generator< Q, 1 >::Type QPoly;
    typedef Polynomial_type_generator< F, 1 >::Type FPoly;
    typedef Polynomial_type_generator< R, 1 >::Type RPoly;
    
    typedef Arcavoid_kernel< IPoly > IAK;
    typedef Arcavoid_kernel< QPoly > QAK;
    typedef Arcavoid_kernel< FPoly > FAK;
    typedef Arcavoid_kernel< RPoly > RAK;
    
    bool exact;
    F x;
    
    I Ifac20 = 1;
    for (int i = 2; i <= 20; ++i)
      Ifac20 *= i;

    IPoly If (Ifac20);
    QPoly Qf (Q (Ifac20, I(1) << 10));
    FPoly f;
    R Rfcs[3] = {
      leda::sqrt (leda_real (2)),
      - leda::root (leda_real (32), 4),
      R(1.),
    };
    RPoly Rf (Rfcs, Rfcs+3);
    
    cerr << "Test leda_integer -> leda_bigfloat" << endl;

    exact = IAK().convert_to_bigfloat_object() (Ifac20, 53, x);
    cerr << x << " = 20!;      exact?: " << exact << " should be 1" << endl;
    exact = IAK().convert_to_bigfloat_object() (Ifac20, 8, x);
    cerr << x << " = 20!;      exact?: " << exact << " should be 0" << endl;
    exact = IAK().convert_polynomial_to_bigfloat_object() (If, 53, f);
    cerr << f << " = 20!;      exact?: " << exact << " should be 1" << endl;
    exact = IAK().convert_polynomial_to_bigfloat_object() (If, 8, f);
    cerr << f << " = 20!;      exact?: " << exact << " should be 0" << endl;

    cerr << "Test leda_rational -> leda_bigfloat" << endl;
    
    exact = QAK().convert_to_bigfloat_object() (Q (I(1), Ifac20), 53, x);
    cerr << x << " = 1/20!;    exact?: " << exact << " should be 0" << endl;
    exact = QAK().convert_to_bigfloat_object() (Q (Ifac20, I(1) << 10), 53, x);
    cerr << x << " = 20!/2^10; exact?: " << exact << " should be 1" << endl;
    exact = QAK().convert_to_bigfloat_object() (Q (Ifac20, I(1) << 10), 8, x);
    cerr << x << " = 20!/2^10; exact?: " << exact << " should be 0" << endl;
    exact = QAK().convert_polynomial_to_bigfloat_object() (Qf, 53, f);
    cerr << f << " = 20!/2^10; exact?: " << exact << " should be 1" << endl;
    exact = QAK().convert_polynomial_to_bigfloat_object() (Qf, 8, f);
    cerr << f << " = 20!/2^10; exact?: " << exact << " should be 0" << endl;

    cerr << "Test leda_bigfloat -> leda_bigfloat" << endl;

    exact = FAK().convert_to_bigfloat_object() (leda::sqrt (leda_bigfloat (2.), 24), 53, x);
    cerr << x << " = sqrt(2);  exact?: " << exact << " should be 1" << endl;
    exact = FAK().convert_to_bigfloat_object() (leda::sqrt (leda_bigfloat (2.), 24), 12, x);
    cerr << x << " = sqrt(2);  exact?: " << exact << " should be 0" << endl;

    cerr << "Test leda_real -> leda_bigfloat" << endl;

    leda_real Rsqrt6 = leda::sqrt (leda_real (6));
    leda_real Rsqrt3 = leda::sqrt (leda_real (3));
    leda_real Rsqrt2 = leda::sqrt (leda_real (2));
    leda_real Rzero = Rsqrt6 / Rsqrt3 - Rsqrt2;

    exact = RAK().convert_to_bigfloat_object() (Rzero, 53, x);
    cerr << x << " = 0;        exact?: " << exact << " should be 1" << endl;
    long prec = 1000;
    // set to some huge number to see whether refinement is called twice
    // prec = 1000000;
    exact = RAK().convert_to_bigfloat_object() (Rsqrt2, prec, x);
    cerr << x << " = sqrt(2);  exact?: " << exact << " should be 0" << endl;
    exact = RAK().convert_to_bigfloat_object() (Rsqrt2, prec / 4 * 3, x);
    cerr << x << " = sqrt(2);  exact?: " << exact << " should be 0" << endl;
    exact = RAK().convert_polynomial_to_bigfloat_object() (Rf, 53, f);
    cerr << f << ";            exact?: " << exact << " should be 0" << endl;
  }
#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_GMP
#ifdef CGAL_USE_MPFR
  cerr << "Test GMP" << endl
       << "================" << endl;
  cerr << "Skipped: Arcavoid_kernel< Polynomial< GMPxx > > gives rubbish" << endl;
  if (false)
  {

    typedef Gmpz I;
    typedef Gmpq Q;
    typedef Gmpfr F;

    typedef Polynomial_type_generator< I, 1 >::Type IPoly;
    typedef Polynomial_type_generator< Q, 1 >::Type QPoly;
    
    typedef Arcavoid_kernel< IPoly > IAK;
    typedef Arcavoid_kernel< QPoly > QAK;

    bool exact;
    F x;
    
    I Ifac20 = 1;
    for (int i = 2; i <= 20; ++i)
      Ifac20 *= i;

    cerr << "Test Gmpz -> Gmpfr" << endl;
    cerr << Ifac20 << endl;
    
    exact = IAK().convert_to_bigfloat_object() (Ifac20, 53, x);
    cerr << x << " = 20!;      exact?: " << exact << " should be 1" << endl;
    exact = IAK().convert_to_bigfloat_object() (Ifac20, 8, x);
    cerr << x << " = 20!;      exact?: " << exact << " should be 0" << endl;
    
    cerr << "Test Gmpq -> Gmpfr" << endl;
    
    exact = QAK().convert_to_bigfloat_object() (Q (I(1), Ifac20), 53, x);
    cerr << Q (I(1), Ifac20) << endl;
    cerr << x << " = 1/20!;    exact?: " << exact << " should be 0" << endl;
    exact = QAK().convert_to_bigfloat_object() (Q (Ifac20, I(1) << 10), 53, x);
    cerr << x << " = 20!/2^10; exact?: " << exact << " should be 1" << endl;
    exact = QAK().convert_to_bigfloat_object() (Q (Ifac20, I(1) << 10), 8, x);
    cerr << x << " = 20!/2^10; exact?: " << exact << " should be 0" << endl;
  }
#endif // CGAL_USE_MPFR
#endif // CGAL_USE_GMP

  return 0;
}
