#undef CGAL_ISOLATE_USE_NT

#ifndef CGAL_ISOLATE_USE_GMP
#define CGAL_ISOLATE_USE_GMP 1 // 1 is default 
#endif
#if CGAL_ISOLATE_USE_GMP
#define CGAL_ISOLATE_USE_NT GMP
#warning Using Gmp number types
#endif

#ifndef CGAL_ISOLATE_USE_CORE
#define CGAL_ISOLATE_USE_CORE 0 // 0 is default 
#endif
#if CGAL_ISOLATE_USE_CORE
#ifdef CGAL_ISOLATE_USE_NT
#error Number type already chosen
#else
#define CGAL_ISOLATE_USE_NT CORE
#warning Using Core number types
#endif
#endif

#ifndef CGAL_ISOLATE_USE_RS
#define CGAL_ISOLATE_USE_RS 0
#endif
#if CGAL_ISOLATE_USE_RS
#if !CGAL_ISOLATE_USE_GMP
#error Rs needs Gmp
#endif
#warning Use Rs
#endif

#ifndef CGAL_ISOLATE_USE_DESCARTES
#define CGAL_ISOLATE_USE_DESCARTES 0
#endif
#if CGAL_ISOLATE_USE_DESCARTES
#warning Use descartes
#endif

#ifndef CGAL_ISOLATE_USE_BITSTREAM_DESCARTES
#define CGAL_ISOLATE_USE_BITSTREAM_DESCARTES 0
#endif
#if CGAL_ISOLATE_USE_BITSTREAM_DESCARTES
#ifndef CGAL_ACK_BITSTREAM_USES_E08_TREE
#define CGAL_ACK_BITSTREAM_USES_E08_TREE 1
#endif
#warning Use bitstream descartes
#endif

#ifndef CGAL_ISOLATE_USE_REAL_BITSTREAM_DESCARTES
#define CGAL_ISOLATE_USE_REAL_BITSTREAM_DESCARTES 0
#endif
#if CGAL_ISOLATE_USE_REAL_BITSTREAM_DESCARTES
#warning Use real bitstream descartes
#endif

#ifndef CGAL_ISOLATE_USE_DESCARTES_STAR
#define CGAL_ISOLATE_USE_DESCARTES_STAR 1
#endif
#if CGAL_ISOLATE_USE_DESCARTES_STAR
#warning Use descartes star
#endif


#ifndef CGAL_ISOLATE_SAVE_BENCHMARKS  
#define CGAL_ISOLATE_SAVE_BENCHMARKS 0 // 0 is default
#endif

#ifndef CGAL_ISOLATE_VERBOSE
#define CGAL_ISOLATE_VERBOSE 1
#endif

#ifndef CGAL_ISOLATE_DEBUG
#define CGAL_ISOLATE_DEBUG 1
#endif

#include <CGAL/config.h>

#include <algorithm>

#include <CGAL/Timer.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Polynomial.h>
#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>
#include <CGAL/Polynomial_type_generator.h>

#if CGAL_ISOLATE_USE_RS
#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>
#endif

#if CGAL_ISOLATE_USE_DESCARTES
#include <CGAL/Algebraic_kernel_d/Descartes.h>
#endif

#if CGAL_ISOLATE_USE_BITSTREAM_DESCARTES
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#endif

#if CGAL_ISOLATE_USE_REAL_BITSTREAM_DESCARTES
#include <CGAL/Algebraic_kernel_d/Real_bitstream_descartes.h>
#endif

#if CGAL_ISOLATE_USE_DESCARTES_STAR
#include <CGAL/Algebraic_kernel_d/Descartes_star.h>
#endif

CGAL::Timer t_total;

template < class Polynomial_, class Rational_ >
void read_from_file(const char *filename, std::vector< Polynomial_ >& polys,
		    Rational_) {

  typedef Polynomial_ Polynomial;
  typedef Rational_ Rational;

  const unsigned d = CGAL::internal::Dimension< Polynomial >::value;
  typedef typename CGAL::Polynomial_type_generator< Rational, d >::Type
    RPoly;

  CGAL::Polynomial_parser_d< RPoly,
			     CGAL::Mixed_rational_parser_policy< RPoly > > parser;
  std::ifstream in(filename);

  if (!in) {
    return;
  }
  while (!in.eof()) {
    Polynomial f;
    std::string s;
    std::getline(in, s);
    
    if (s.length() == 0) {
      continue;
    }

    RPoly rf;
    if (parser(s, rf)) {
      typedef typename CGAL::Fraction_traits< RPoly > FT;
      typename FT::Denominator_type _;
      typename FT::Decompose()(rf, f, _);
    } else {
      //std::cerr << "Parser error, trying another format..\n";
      try {
	std::stringstream ss(s);
	ss >> f;
      } catch(...) {
	std::cerr << "Invalid format of polynomial..skipping..\n";
	continue; 
      }
    }
    polys.push_back(f);
  }
}

template < class Coefficient_, class Rational_ >
int isolate(int argc, char **argv) {

  typedef Coefficient_ Coefficient;
  typedef Rational_ Rational;

  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);
  CGAL::set_pretty_mode(std::clog);
  
#if CGAL_ISOLATE_USE_RS
  typedef CGAL::Algebraic_kernel_rs_gmpz_d_1 Rs_ak_1;

  typedef typename Rs_ak_1::Polynomial_1 Polynomial_1;
#else
  typedef typename CGAL::Polynomial_type_generator< Coefficient, 1 >::Type Polynomial_1;
#endif

  typedef CGAL::Polynomial_traits_d< Polynomial_1 > PT_1;
  
  typedef typename PT_1::Coefficient_const_iterator_range 
    Coefficient_const_iterator_range;

  typedef typename PT_1::Construct_coefficient_const_iterator_range 
    Construct_coefficient_const_iterator_range;

  // we deal with two polynomials two bivariate polynomials
  Polynomial_1 f;
  
  if (argc < 2) {
    return 1;
  }
  
  std::ifstream in(argv[1]);
  std::vector< Polynomial_1 > polys;
  read_from_file(argv[1], polys, Rational());
  
  if (polys.size() == 0) {
    std::cerr << "Grr.. I need at least one polynomial..!!\n";
    return 1;
  }
  
#if CGAL_ISOLATE_SAVE_BENCHMARKS  
  std::ofstream os("benchmarks_out", std::ios::app);
  os << "#################################################################";
  os << "\n USE_GMP: " <<
    CGAL_ISOLATE_USE_GMP;                
  os << "\n USE_CORE: " <<
    CGAL_ISOLATE_USE_CORE;                
  //  os << "\n\n******************* POLYNOMIAL: " << argv[1] << "\n" << std::flush;
  CGAL::set_pretty_mode(os);
#else
  std::ostream& os = std::cout;
#endif    
  
  // TODO iterator over all polynomials in a file
  f = polys[0];
  f = CGAL::canonicalize(f);
  
  // TODO compile only on demand
  std::cerr << "Making square-free ... " << std::flush;
  f = CGAL::make_square_free(f);
  std::cerr << "done." << std::endl;

  //#if CGAL_ISOLATE_SAVE_BENCHMARKS  
  std::cout << "f:= " << f << ":\n";
  //#endif
  

#if CGAL_ISOLATE_USE_RS
  
  typedef typename Rs_ak_1::Algebraic_real_1 Algebraic_real_1;
  typedef typename Rs_ak_1::Multiplicity_type Multiplicity_type;

  std::list< std::pair< Algebraic_real_1, Multiplicity_type > > solutions;

  Rs_ak_1 ak1;
  
  t_total.start();
  ak1.solve_1_object()(f,std::back_inserter(solutions));
  t_total.stop();

#else

#if CGAL_ISOLATE_USE_DESCARTES
  typedef ::CGAL::internal::Descartes< Polynomial_1, Rational > Isolator;
#endif

#if CGAL_ISOLATE_USE_BITSTREAM_DESCARTES
  typedef CGAL::internal::Bitstream_descartes<
    CGAL::internal::Bitstream_descartes_rndl_tree_traits
    <CGAL::internal::Bitstream_coefficient_kernel< Coefficient > > > Isolator;
#endif

#if CGAL_ISOLATE_USE_REAL_BITSTREAM_DESCARTES
  typedef ::CGAL::internal::Real_bitstream_descartes< Coefficient, Rational > Isolator;
#endif

#if CGAL_ISOLATE_USE_DESCARTES_STAR
  typedef ::CGAL::internal::Descartes_star< Coefficient, Rational > Isolator;
#endif
  
#if CGAL_ISOLATE_USE_DESCARTES || CGAL_ISOLATE_USE_BITSTREAM_DESCARTES

  t_total.start();
  // TODO range constructor for classical (bitstream) descartes isolators?
  Isolator isolator(f);
  isolator.number_of_real_roots(); // ensure computation (avoid lazy construction)
  // TODO is #real_roots enough
  t_total.stop();

#else

  Coefficient_const_iterator_range range = 
    Construct_coefficient_const_iterator_range()(f);  
  
  t_total.start();
  Isolator isolator(range.first, range.second);
  isolator.number_of_real_roots(); // ensure computation (avoid lazy construction)
  // TODO is #real_roots enough
  t_total.stop();

#endif

#endif // !CGAL_ISOLATE_USE_RS
  
#if 0 && CGAL_ISOLATE_SAVE_BENCHMARKS  
  // TODO correct savings
  std::stringstream ss;
  CGAL::set_pretty_mode(ss);
  char *fname = basename(argv[1]);
  ss << "\n" << "f_" << fname << ":=" << f << ":\n\n";
  ss << "#####################################################\n";
  std::ofstream o("maple_out", std::ios::app);
  o << ss.str();
  o.close();
#endif
  
#if 0
  // TODO more outout
  os.precision(25);
  os << "\nSolutions: " << std::endl;
  for (typename Solution_list::const_iterator
	 it = solutions.begin(); it != solutions.end(); it++) {
#if CGAL_ISOLATE_USE_RS
    os << "x: " << CGAL::to_double(it->first.x()) << ", arcno: " <<
      it->first.arcno() << " (" << it->second << ")\n";
#else
    os << *it << std::endl;
#endif
  }
#endif

#if CGAL_ISOLATE_USE_RS
  os << "#SOLUTIONS: " << solutions.size() << "\n\n";
#else
  os << "#SOLUTIONS: " << isolator.number_of_real_roots() << "\n\n";
#endif  

  os.precision(7);
  os << "TIME total: " << t_total.time() << std::endl;
  
  os << "\n###########################################################\n";

  return 0;
}

int main(int argc, char **argv) {
#if CGAL_ISOLATE_USE_GMP
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  {
    std::cout << " Test ISOLATE using Gmp " << std::endl;
    typedef CGAL::GMP_arithmetic_kernel AK;
    isolate< AK::Integer, AK::Rational >(argc, argv);
  }
#else
#warning Gmp NOT AVAILABLE! Use other number type library
#endif
#endif
  
#if CGAL_ISOLATE_USE_CORE
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  {
    std::cout << " Test ISOLATE using Core " << std::endl;
    typedef CGAL::CORE_arithmetic_kernel AK;
    isolate< AK::Integer, AK::Rational >(argc, argv);
  }
#else
#warning Core NOT AVAILABLE! Use other number type library
#endif
#endif
  // #ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
  //   {
  //     std::cout << " Test ISOLATE using Leda " << std::endl;
  //     typedef CGAL::LEDA_arithmetic_kernel AK;
  //     isolate< AK::Integer, AK::Rational >(argc, argv);
  //   }
  // #endif
  
  return 0;
}
