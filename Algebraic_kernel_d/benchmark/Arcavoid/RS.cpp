#include <vector>
#include <iostream>
#include <boost/lexical_cast.hpp>

#include <CGAL/basic.h>
#include <CGAL/GMP_arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>

// Includes for benchmarks
#include <CGAL/Timer.h>

// RS Solver
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpfi.h>
#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_d_1.h>

int main (int argc, char **argv) {
  using namespace std;

  typedef CGAL::GMP_arithmetic_kernel AK;
  typedef AK::Integer Integer;
  typedef CGAL::Polynomial_type_generator< Integer, 1 >::Type IPoly_1;
  typedef CGAL::Timer Timer;

  // READ POLYNOMIAL
  vector< string > args (argv, argv+argc);
  istream *in_ptr = &cin;
  if (args.size() > 1 && args.back() != "--")
    in_ptr = new ifstream (args.back().c_str());
  string line;
  while (line.empty() || line[0] == '#')
    getline (*in_ptr, line);

  IPoly_1 f_int = boost::lexical_cast< IPoly_1 > (line);
  
  ////////////////////////////////////////////////////////////////
  // The interesting part
  ////////////////////////////////////////////////////////////////

  typedef CGAL::Algebraic_kernel_rs_gmpz_d_1 AK_RS;
  typedef AK_RS::Polynomial_1                RS_Poly_1;
  typedef AK_RS::Algebraic_real_1            RS_AR_1;
  typedef AK_RS::Solve_1                     Solve_1;
  
  Timer t;
  t.start();

  vector< RS_AR_1 > rs_roots;
  AK_RS ak_rs;
  Solve_1 solve_1 = ak_rs.solve_1_object();
  
  set_rs_verbose(3);

  solve_1 (f_int, true, back_inserter (rs_roots));
  cerr << "[ " << rs_roots.size() << " isolating intervals found ]" << endl;
  for (int i = 0; i < rs_roots.size(); ++i)
    cerr << "r[" << i << "] = " << to_double (rs_roots[i]) << endl;

  t.stop();
  cerr << "Time: " << t.time() << " sec" << endl;

  // CLEAN UP
  if (in_ptr && in_ptr != &cin)
    delete in_ptr;

  return 0;
}
