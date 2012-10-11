#include <vector>
#include <iostream>
#include <boost/lexical_cast.hpp>

#include <CGAL/basic.h>
#include <CGAL/GMP_arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_type_generator.h>

// Includes for benchmarks
#include <CGAL/Timer.h>

// Descartes
#include <CGAL/Algebraic_kernel_d/Descartes.h>

// #include "include/RArcavoid.h"

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

  typedef CGAL::internal::Descartes
    < IPoly_1, AK::Rational >
    Isolator;

  // typedef CGAL::RArcavoid< AK > Isolator;
  
  Timer t;
  t.start();

  Isolator isolator (f_int);
  cerr << "[ " << isolator.number_of_real_roots() << " isolating intervals found ]" << endl;
  for (int i = 0; i < isolator.number_of_real_roots(); ++i)
    cout << "r[" << i << "] in (" << isolator.left_bound (i) << ", " << isolator.right_bound (i) << ")" << endl;  

  t.stop();
  cerr << "Time: " << t.time() << " sec" << endl;

  // CLEAN UP
  if (in_ptr && in_ptr != &cin)
    delete in_ptr;

  return 0;
}
