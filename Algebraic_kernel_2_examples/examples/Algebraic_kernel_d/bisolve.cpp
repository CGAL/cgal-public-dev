#include <CGAL/config.h>

#ifndef CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING
#define CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING 0 // 0 is default
#endif

#define CGAL_BISOLVE_ENABLE_ARCAVOID 0 // default TODO?
#define CGAL_ACK_BITSTREAM_USES_E08_TREE 1 // do not change
#define CGAL_BISOLVE_USE_ADJUSTABLE_PRECISION 1 // 1 is default

#define CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER 1
#define CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE 1

#define CGAL_BISOLVE_DEBUG 0
#define CGAL_BISOLVE_VERBOSE 0

#define CGAL_BISOLVE_USE_RS_AK 0
#define CGAL_BISOLVE_USE_RS_ISOLATOR 1 // 1 is default

#define CGAL_ACK_DEBUG_FLAG 0
#define CGAL_ACK_DEBUG_PRINT std::cout

#define CGAL_BISOLVE_USE_RESULTANT_COFACTORS 1

#define CGAL_BISOLVE_USE_GMP 1
#define CGAL_BISOLVE_USE_CORE 0

#ifndef CGAL_BISOLVE_USE_GPU_RESULTANTS
#define CGAL_BISOLVE_USE_GPU_RESULTANTS 1 // default?
#define CGAL_BISOLVE_CHECK_GPU_RESULTANTS_SANITY 0 // default 0
#endif

#ifndef CGAL_BISOLVE_USE_GPU_GCDS
#define CGAL_BISOLVE_USE_GPU_GCDS 1  // default?
#define CGAL_BISOLVE_CHECK_GPU_GCDS_SANITY 0 // default 1
#endif

#ifndef CGAL_BISOLVE_USE_BIGCD
#define CGAL_BISOLVE_USE_BIGCD 1
#define CGAL_BIGCD_USE_SHIFT 0
#define CGAL_BIGCD_CHECK_SANITY 0
#endif

#ifndef CGAL_BISOLVE_USE_NTL
#define CGAL_BISOLVE_USE_NTL  1 // default 1 ??
#endif

#include <CGAL/symbolic_standalone.h>


#undef CGAL_BISOLVE_USE_NT

#ifndef CGAL_BISOLVE_USE_GMP
#define CGAL_BISOLVE_USE_GMP 1 // 1 is default
#endif
#if CGAL_BISOLVE_USE_GMP
#define CGAL_BISOLVE_USE_NT GMP
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING
#warning bisolve: Using Gmp number types
#endif
#endif

#ifndef CGAL_BISOLVE_USE_CORE
#define CGAL_BISOLVE_USE_CORE 0 // 0 is default
#endif
#if CGAL_BISOLVE_USE_CORE
#ifdef CGAL_BISOLVE_USE_NT
#error Number type already chosen
#else
#define CGAL_BISOLVE_USE_NT CORE
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING
#warning bisolve: Using Core number types
#endif
#endif
#endif


#ifndef CGAL_BISOLVE_USE_RS_AK
#define CGAL_BISOLVE_USE_RS_AK 0
#endif
#if CGAL_BISOLVE_USE_RS_AK
#if !CGAL_USE_GMP
#error Rs needs Gmp
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING
#warning bisolve: Using univariate algebraic kernel by Rs
#endif
#endif

#ifndef CGAL_BISOLVE_USE_RS_ISOLATOR
#define CGAL_BISOLVE_USE_RS_ISOLATOR 1 // 1 is default
#endif
#if CGAL_BISOLVE_USE_RS_ISOLATOR
#if !CGAL_USE_GMP
#error Rs needs Gmp
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING
#warning bisolve: Using Rs Solve_1
#endif
#endif



// global work
#ifndef CGAL_BISOLVE_ARRANGEMENTS
#define CGAL_BISOLVE_ARRANGEMENTS 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_ARRANGEMENTS
#warning bisolve: Arrangements are enabled
#endif

// multiply all curves into a combined one
#ifndef CGAL_BISOLVE_COMBINE_CURVES
#define CGAL_BISOLVE_COMBINE_CURVES 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_COMBINE_CURVES
#warning bisolve: Multiplying curves to combined one enabled
#endif


// processing several input curves: "single curves + f_y" or "all pairs"?
#ifndef CGAL_BISOLVE_SINGLE_CURVES
#define CGAL_BISOLVE_SINGLE_CURVES 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_SINGLE_CURVES
#warning bisolve: Single curves enabled
#endif

#ifndef CGAL_BISOLVE_CURVE_ANALYSES
#define CGAL_BISOLVE_CURVE_ANALYSES 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_CURVE_ANALYSES
#warning bisolve: Curve analyses are enabled
#endif

#ifndef CGAL_BISOLVE_USE_AK2 
#define CGAL_BISOLVE_USE_AK2 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_USE_AK2
#warning bisolve: Using AK2
#endif

#ifndef CGAL_BISOLVE_USE_RESULTANT_COFACTORS
#define CGAL_BISOLVE_USE_RESULTANT_COFACTORS 1 // 1 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && !CGAL_BISOLVE_USE_RESULTANT_COFACTORS
#warning bisolve: Using SUBDIVISION approach
#endif

#ifndef CGAL_BISOLVE_SHEAR_INPUT
#define CGAL_BISOLVE_SHEAR_INPUT 0 // 0 is default
#endif
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING && CGAL_BISOLVE_SHEAR_INPUT
#warning bisolve: Shearing input
#endif

// need to set/change flags of Algebraic_kernel_d 

#ifndef CGAL_BISOLVE_SAVE_BENCHMARKS
#define CGAL_BISOLVE_SAVE_BENCHMARKS 0 // 0 is default
#endif

#ifndef CGAL_BISOLVE_DEBUG
#define CGAL_BISOLVE_DEBUG 1 // 1 is default
#endif

#ifndef CGAL_BISOLVE_VERBOSE
#define CGAL_BISOLVE_VERBOSE 1 // 1 is default
#endif

#ifndef CGAL_BISOLVE_TELEMETRY
#define CGAL_BISOLVE_TELEMETRY 1 // should be 1 for this source
#endif

#ifndef CGAL_BISOLVE_WRITE_ONLY_RESULTANTS
#define CGAL_BISOLVE_WRITE_ONLY_RESULTANTS 0 // 0 is default
#endif

#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_2/Bi_solve_2_flags.h>
#include <iostream>
#include <fstream>

#include <CGAL/Timer.h>

CGAL::Timer t_cpa, t_msolve;

CGAL::Timer t_total,    // total running time
  // Project
  t_res[3],                // resultants
  t_sqfree[3],             // making square-free
  t_solve1[3],             // 1D solutions
  t_minprec[3],            // ensure minimal precision of isolating approximations

  // Seperate
  t_wellsep[3],            // well sepersting intervals
  t_refine_app,         // refining solutions
  t_filter,             // lense filter
  t_tshift,             // taylor shift
  t_ttest,              // T-test
  t_sort[3],               // sorting time

  // Validate 
  t_validate,            // lifting
  t_vl_detect,          // time to detect vertical lines
  t_ai_init,            // init active intervals
  t_arca_transform,     // transform polynomial
  t_ai_prep,            // prep active intervals
  t_ai_prep_subdiv,     // prep active intervals subdivision
  t_main_loop,          // main loop total
  t_combinatorial,      // combinatorial tests
  t_arca_cert,          // numerical certification
  t_exin,               // time for paper's inclusion and exclusion predicat
  t_exia2,              // exclusion time with 2D interval evaluation
  t_innorm,             // inclusion tiem with norm test
  t_approx,             // approximating solution to given precision
  t_subs,               // subsituting solutions to polynomial equations
  t_bounds,             // computing initial bounds (cofactors & resultant)
  t_ai_main_subdiv,     // main active intervals subdivison

  t_cel[3],                // create event line
  t_cel_vl,                // create event line, vertical line detect
  t_lcoeff,                // leading coeffiecient
  t_mult_rxy,              // mult res f_x f_y
  t_mult_xy_gcd,           // gcd(f_x,f_y)
  t_mult_res,              // mult cpu resultant
  t_mult_ntl_div,          // ntl div
  t_mult_sqfs,             // mult sqf factorization
  t_mult_sign,             // mult sign
  t_ffy_isol[3],           // ffyisolation
  t_ffyi_tes_isol[3],      // ffyisolator Teissier isolate
  t_ffyi_ak1_isol,         // exact ak1 isolate
  t_ffyi_bdv_isol,         // ffyisolator bisolve isolate
  t_bdv,                   // bidiffvanish 
  t_bdv_construct,         // construction
  t_bdv_bs0,               // bidiffvanish initial bisolve
  t_common,                // common factors of diffs
  t_bdv_bs,                // bidiffvanish further bisolves
  t_ffyi_bits_isol,        // ffyisolator bitstream isolate
  t_sllift[3],             // status line lifts
  t_sllift_init,           // event line init
  t_sllift_generic,        // event line generic
  t_sllift_bucket,         // event line bucket
  t_sllift_refine,         // event line refine
  t_sllift_assign,         // event line assign
  t_sllift_finish,         // event line finish
  t_ffy_refine, t_ffy_ais_subdivide,

  // misc                                                              
  t_sqrt,               // approximating square root
  t_timer, t_ais;              // general purpose timer


//CGAL::Timer tm_eval, tm_exact, tm_cvt;

int c_approx;
int c_ieval2_reject;
int c_t3over2_total;
int c_lensefilter_pass;

#if CGAL_BISOLVE_ENABLE_ARCAVOID
int c_arca_cert_total;
int c_arca_cert_success;
#endif

int c_status_lines_total[3];
int c_status_lines_tes[3];

long g_max_prec;


#if CGAL_BISOLVE_USE_GMP
#include <CGAL/GMP_arithmetic_kernel.h>
#endif

#if CGAL_BISOLVE_USE_CORE
#include <CGAL/CORE_arithmetic_kernel.h>
#endif

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Algebraic_kernel_2/Rounding_ak_d_1.h>

#if CGAL_BISOLVE_USE_RS_AK
#include <CGAL/Algebraic_kernel_d/Float_traits.h>
#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>
#else
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d_1_generator.h>
#endif

#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

#if CGAL_BISOLVE_SHEAR_INPUT
#include <CGAL/Algebraic_kernel_d/Shear_controller.h>
#include <CGAL/Algebraic_kernel_d/shear.h>
#endif

#if CGAL_BISOLVE_USE_RESULTANT_COFACTORS
#if CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
#include <CGAL/Algebraic_kernel_2/Certifier_cofactor_traits.h>
#else
#if CGAL_BISOLVE_ENABLE_ARCAVOID
#include <CGAL/Algebraic_kernel_2/Certifier_cofactor_arcavoid_traits.h>
#endif
// TODO DISABLE ACTIVE_INTERVAL FILTER
#if !CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER
#include <CGAL/Algebraic_kernel_2/Certifier_cofactor_bitstream_traits.h>
#endif
#endif
#else
#include <CGAL/Algebraic_kernel_2/Certifier_subdivision_traits.h>
#endif

#include <CGAL/Algebraic_kernel_2/Bi_solve_2.h>

#if CGAL_BISOLVE_USE_AK2 || CGAL_BISOLVE_CURVE_ANALYSES || CGAL_BISOLVE_ARRANGEMENTS
#include <CGAL/Algebraic_kernel_d_2.h>
#endif
#if CGAL_BISOLVE_CURVE_ANALYSES || CGAL_BISOLVE_ARRANGEMENTS
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>
#endif
#if CGAL_BISOLVE_ARRANGEMENTS
#include <CGAL/Arrangement_2.h>
#endif


template < class Poly_, class Rational_ >
void read_from_file(const char *filename, std::vector< Poly_ >& polys,
		    Rational_) {

  const unsigned d = CGAL::internal::Dimension< Poly_ >::value;
  typedef typename CGAL::Polynomial_type_generator< Rational_, d >::Type
    RPoly;

  CGAL::Polynomial_parser_d< RPoly,
			     CGAL::Mixed_rational_parser_policy< RPoly > > parser;
  std::ifstream in(filename);

  if(!in)
    return;
  while(!in.eof()) {
    Poly_ f;
    std::string s;
    std::getline(in, s);

    if(s.length() == 0)
      continue;

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
        // 	std::cerr << "Invalid format of polynomial..skipping..\n";
	    continue; 
      }
    }
    polys.push_back(f);
    CGAL::set_ascii_mode(std::cout);
    std::cout << f << std::endl;
    CGAL::set_pretty_mode(std::cout);
  }
}

template < class AK_1 >
int bisolve(int argc, char **argv) {

  c_approx = 0;
  g_max_prec = 2;
  c_ieval2_reject = 0;
  c_t3over2_total = 0;
  c_lensefilter_pass = 0;
  
#if CGAL_BISOLVE_ENABLE_ARCAVOID
  c_arca_cert_total = 0;
  c_arca_cert_success = 0;
#endif

  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);
  CGAL::set_pretty_mode(std::clog);
  
  typedef typename AK_1::Coefficient Coefficient;
  typedef typename AK_1::Bound Bound;
  typedef typename AK_1::Polynomial_1 Polynomial_1;

  //! type of Certifier traits
#if CGAL_BISOLVE_USE_RESULTANT_COFACTORS
#if CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
  typedef CGAL::internal::Certifier_cofactor_traits< AK_1 > Certifier_traits;
#else
#if !CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER
  typedef CGAL::internal::Certifier_cofactor_bitstream_traits< AK_1 > Certifier_traits;
#endif
#if CGAL_BISOLVE_ENABLE_ARCAVOID
  typedef CGAL::internal::Certifier_cofactor_arcavoid_traits< AK_1 > Certifier_traits;
#endif
#endif
#else
  typedef CGAL::internal::Certifier_subdivision_traits< AK_1 > Certifier_traits;
#endif 

  typedef CGAL::Bi_solve_2< Certifier_traits > Bi_solve_2;

  typedef typename Bi_solve_2::Polynomial_2 Polynomial_2;
  typedef typename Bi_solve_2::Multiplicity_type Multiplicity_type;

#if CGAL_BISOLVE_USE_AK2 || CGAL_BISOLVE_CURVE_ANALYSES || CGAL_BISOLVE_ARRANGEMENTS
  // TODO use Arr_algebraic_segment_traits< Coefficient > :: AK_1, if its "SLOW" default AK_1 has been replaced
  typedef CGAL::Algebraic_curve_kernel_2< AK_1 > AK_2;
#endif
#if CGAL_BISOLVE_CURVE_ANALYSES || CGAL_BISOLVE_ARRANGEMENTS
  typedef CGAL::Curved_kernel_via_analysis_2< AK_2 > Arr_traits_2;
  typedef typename Arr_traits_2::Curve_2 Curve_2;
#endif
#if !CGAL_BISOLVE_CURVE_ANALYSES && CGAL_BISOLVE_ARRANGEMENTS
  typedef CGAL::Arrangement_2< Arr_traits_2 > Arrangement_2;
#endif

  // we deal with two polynomials two bivariate polynomials
  Polynomial_2 f, g;

  if(argc < 2) {
    return 1;
  }

  std::ifstream in(argv[1]);
  std::vector< Polynomial_2 > polys;
  read_from_file(argv[1], polys, Bound());

  //std::cout << "#polynomials: " << polys.size() << std::endl;

  if (polys.size() == 0) {
    std::cerr << "Grr.. I need at least one polynomial..!!\n";
    return 1;
  }
  
#if 0 // convert to maple format
  std::ifstream in(argv[1]);
  std::vector< Polynomial_2 > polys;
  read_from_file(argv[1], polys, Rational());
  
  std::cout << "write file"<< std::flush;
  std::ofstream oo(strcat(argv[1], ".mpl")); // ofstream has buffer overflow so use
  //std::ofstream oo(argv[1]); // ofstream has buffer overflow so use
  std::stringstream ss;
  CGAL::set_pretty_mode(ss);
  for (int i = 0; i < polys.size(); i++) {
    ss << polys[i] << "\n";
  }
  oo << ss.str();            // strings first
  oo.close();
  std::cout << " done"  << std::endl;
  std::exit(0);
#endif

#if CGAL_BISOLVE_SAVE_BENCHMARKS  
  std::ofstream os("benchmarks_out", std::ios::app);
  os << "#################################################################";
  os << "\n USE_GMP: " <<
    CGAL_BISOLVE_USE_GMP;                
  os << "\n USE_CORE: " <<
    CGAL_BISOLVE_USE_CORE;                
  os << "\n USE_AK2: " <<
    CGAL_BISOLVE_USE_AK2;                
  os << "\n USE_RS_AK: " <<
    CGAL_BISOLVE_USE_RS_AK;              
  os << "\n USE_GPU_RESULTANTS: " <<
    CGAL_BISOLVE_USE_GPU_RESULTANTS;
  os << "\n USE_RS_ISOLATOR: " <<
    CGAL_BISOLVE_USE_RS_ISOLATOR;              
  os << "\n ENABLE_LENSE_FILTER_IN_T_TEST: " <<
    CGAL_BISOLVE_ENABLE_LENSE_FILTER_IN_T_TEST;
  os << "\n DISABLE_COMBINATORIAL_CERTIFICATION: " <<
    CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION;
  os << "\n DISABLE_BIDIRECTIONAL_CERTIFICATION: " <<
    CGAL_BISOLVE_DISABLE_BIDIRECTIONAL_CERTIFICATION;
  os << "\n\n******************* CURVE: " << argv[1] << "\n" << std::flush;
  CGAL::set_pretty_mode(os);
#else
  std::ostream& os = std::cout;
#endif    

#if CGAL_BISOLVE_COMBINE_CURVES
  Polynomial_2 combined = polys[0];
  for (int i = 1; i < polys.size(); i++) {
    combined *= polys[i];
  }
  polys.clear();
  polys.push_back(combined);
  std::cout << "Combined curve: " << combined << std::endl;
  std::cout << "Degree of combined curve: " << CGAL::degree(combined) << std::endl;
#endif
    
#if !CGAL_BISOLVE_SINGLE_CURVES
  if (polys.size() == 1) {
    polys.push_back(CGAL::differentiate(polys[0])); 
    std::cout << "Only one polynomial found; "
              << "use derivative as second polynomial." << std::endl;
  } 
#endif

#if CGAL_BISOLVE_DEBUG
  std::vector< Polynomial_1 > res_debug; 
  if(argc >= 3) { // try read out resultants from file
    read_from_file(argv[2], res_debug, Bound());
    if(res_debug.size() < 2) {
      std::cerr << "oops.. you have touched the very fabric of "
	"space-time causing the Universe to cease existence..\n";
    }
  }
  // view with extended precision because sometimes double-precision accuracy
  // is not enough
  std::clog.precision(30);
#endif

#if CGAL_BISOLVE_SHEAR_INPUT
  CGAL::internal::Shear_controller< Coefficient > shear_controller;
  Coefficient s(0);
  try {
    s = shear_controller.get_shear_factor();
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "Trying shear factor " 
                         << s << std::endl;
#endif
    
    std::vector< Polynomial_2> spolys;
    for (typename std::vector< Polynomial_2 >::iterator pit = polys.begin(); pit != polys.end(); pit++) {
        
      Polynomial_2 sh_pol = CGAL::internal::shear(*pit, s);
      //     if (CGAL::degree(CGAL::univariate_content(sh_pol)) > 0) {
        // TODO or if not a good shearing for another reason throw error!
      //  throw CGAL::internal::Non_generic_position_exception();
      //}
      
      spolys.push_back(sh_pol);
    }
    
    polys = spolys;
    
  } catch(CGAL::internal::Non_generic_position_exception /* err */) {
    shear_controller.report_failure(s);
#if CGAL_ACK_DEBUG_FLAG
    CGAL_ACK_DEBUG_PRINT << "Bad shear factor, retrying..." 
                         << std::endl;
#endif
  }
#endif

  for (typename std::vector< Polynomial_2 >::iterator pit = polys.begin(); pit != polys.end(); pit++) {
    *pit = CGAL::canonicalize(*pit);
    std::cout << *pit << std::endl;
  }
  
  int n_pairs = 0; // must be n(n-1) / 2

  // PREP END


#if CGAL_BISOLVE_ARRANGEMENTS
 
  Arrangement_2 arr;
  
  // Functor to create a curve from a Polynomial_2
  typename AK_2::Construct_curve_2 construct_curve = arr.geometry_traits()->kernel().construct_curve_2_object();

  std::list< Curve_2 > curves;

  for (typename std::vector< Polynomial_2 >::iterator pit = polys.begin(); pit != polys.end(); pit++) {
    curves.push_back(construct_curve(*pit));
  }
  
  t_total.start(); // TODO wrong timer as only "valid" for bisolve (used in substractions)
  CGAL::insert(arr, curves.begin(), curves.end());
  t_total.stop();


  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << " (" << arr.number_of_unbounded_faces() << ")"
            << std::endl;
  // condensed, for exp-run
  std::cout << "#SOLUTIONS: "
            << "V=" << arr.number_of_vertices()
            << ",E=" << arr.number_of_edges() 
            << ",F=" << arr.number_of_faces() << "(" << arr.number_of_unbounded_faces() << ")"
            << std::endl;
    
#else // no arrangements

#if CGAL_BISOLVE_USE_AK2

  std::cout << " Test BISOLVE using AK2 " << std::endl;

#endif

#if CGAL_BISOLVE_CURVE_ANALYSES

  std::cout << " Test BISOLVE using Curve Analyses " << std::endl;

#endif
  
#if CGAL_BISOLVE_USE_AK2 || CGAL_BISOLVE_CURVE_ANALYSES

  typedef typename AK_2::Algebraic_real_2 Algebraic_real_2;

  typedef std::list<
    std::pair< Algebraic_real_2, typename AK_2::Multiplicity_type >
    >
    Solution_list;
  
  Solution_list solutions;
  
  AK_2 ak2;

#else // CGAL_BISOLVE_USE_AK2

  typedef typename Bi_solve_2::Algebraic_real_2 Algebraic_real_2;
  typedef std::list< Algebraic_real_2 > Solution_list;

  Solution_list solutions;

#endif  // CGAL_BISOLVE_USE_AK2

  //////////////////
  // MAIN LOOP START

#if CGAL_BISOLVE_SINGLE_CURVES

  for (typename std::vector< Polynomial_2 >::iterator pit = polys.begin(); 
       pit != polys.end(); pit++) {
    {
      
      f = *pit;
      g = CGAL::differentiate(f); // choose diff
      
#else
      
  for (typename std::vector< Polynomial_2 >::iterator poit = polys.begin();  
    CGAL::cpp0x::next(poit) != polys.end(); poit++) {
    for (typename std::vector< Polynomial_2 >::iterator piit = CGAL::cpp0x::next(poit); 
         piit != polys.end(); piit++)
    {
        
      f = *poit;
      g = *piit; // choose next
      
#endif
      
      // std::cout << "f:  " << f << std::endl;
      // std::cout << "g:  " << g << std::endl;
      
      n_pairs++;

#if CGAL_BISOLVE_CURVE_ANALYSES

#if !CGAL_BISOLVE_SINGLE_CURVES && !CGAL_BISOLVE_COMBINE_CURVES
#error "Curve analyses need single curves"
#endif

      t_total.start();
      std::cout << "*************** Analysis starts ***************" << std::endl;
      
      AK_2 traits;
      
      Curve_2 curve = traits.construct_curve_2_object()(f);
      
      std::cout << "Now refine..." << std::flush;
      
      curve.refine_all(Bound(1,100));
      
      std::cout << "done" << std::endl;
      
      std::cout << "*************** Analysis finished ***************" << std::endl;
      
      for (typename Curve_2::Event_line_iterator eli = curve.event_begin(); eli != curve.event_end(); eli++) {
        for (int j = 0; j < eli->number_of_events(); j++) {
          if (eli->covers_line() || eli->is_event(j)) {
            solutions.push_back(std::make_pair(eli->algebraic_real_2(j),0));
          }
        }
      }

      t_total.stop();        

      std::cout << "Overall  time: " << t_total.time() << std::endl;

      std::cout << curve;
      
#else // CGAL_BISOLVE_CURVE_ANALYSES

#if CGAL_BISOLVE_USE_AK2
      t_total.start();
      ak2.solve_2_object()(f,g, std::back_inserter(solutions));
      t_total.stop();
#else  // CGAL_BISOLVE_USE_AK2
      
#if CGAL_BISOLVE_USE_GPU_RESULTANTS
      // now preallocate memory
      CGAL::internal::GPU_resultant& res_obj =
        CGAL::internal::GPU_resultant::instance();
      (void)res_obj;
#endif
#if CGAL_BISOLVE_USE_GPU_GCDS
//       // now preallocate memory
      CGAL::internal::GPU_gcd& gcd_obj =
          CGAL::internal::GPU_gcd::instance();
      (void)gcd_obj; // allocate gpu mem
#endif

      Bi_solve_2 bi_solve = (CGAL_BISOLVE_SINGLE_CURVES ? Bi_solve_2(f) : Bi_solve_2(f,g));

#if CGAL_BISOLVE_DEBUG
      if(res_debug.size() >= 2) { // load resultants from file
        // first: res_x = res(f, g, y)
        // second: res_y = res(f, g, x) = res(ft, gt, y)
        bi_solve.set_resultant_in_x(res_debug[0]);
        bi_solve.set_resultant_in_y(res_debug[1]);
      }
#endif  
      
#if 0
      Bound xmin(-1);
      Bound xmax( 1);
      
      Bound ymin(-1);
      Bound ymax( 1);
      
      typename Bi_solve_2::Box_2 box = CGAL::make_array(xmin, xmax, ymin, ymax);
      
      t_total.start();
      bi_solve(box, std::back_inserter(solutions));
      t_total.stop();
      
#else // box
      
      t_total.start();
      bi_solve(std::back_inserter(solutions));
      t_total.stop();
      
#endif // box
      
#endif // !CGAL_BISOLVE_USE_AK2
      
#endif // !CGAL_BISOLVE_CURVE_ANALYSES

    }
  }

  // MAIN LOOP END
  ////////////////

#endif // !CGAL_BISOLVE_ARRANGEMENTS

#if CGAL_BISOLVE_SAVE_BENCHMARKS  
//   std::stringstream ss;
//   CGAL::set_pretty_mode(ss);
//   char *fname = basename(argv[1]);
//   ss << "[[" << f << ", " << g << "], \"" <<  fname << "\"], \n\n";
//   std::ofstream o("maple_out", std::ios::app);
//   o << ss.str();
//   o.close();
//   return 1;

  std::stringstream ss;
  CGAL::set_pretty_mode(ss);
  char *fname = basename(argv[1]);

  ss << "[[";
  for (typeof(polys.begin()) poit = polys.begin(); poit != polys.end(); poit++) {
    ss << *poit;
    if(poit != polys.end()-1)
      ss << ',';
  }
  ss << "], \"" <<  fname << "\"], \n\n";
  std::ofstream o("maple_out", std::ios::app);
  o << ss.str();
  o.close();
  return 1;    
#endif

#if CGAL_BISOLVE_ARRANGEMENTS

  // TODO output for ARR

#else // CGAL_BISOLVE_ARRANGEMENTS
  
  os.precision(25);
  os << "\nSolutions: " << std::endl;
  for (typename Solution_list::const_iterator
	 it = solutions.begin(); it != solutions.end(); it++) {
#if CGAL_BISOLVE_USE_AK2 || CGAL_BISOLVE_CURVE_ANALYSES
    os << "x: " << CGAL::to_double(it->first.x()) << ", arcno: " <<
      it->first.arcno() << " (" << it->second << ")\n";
#else
    os << *it << std::endl;
#endif
  }
  os << "#SOLUTIONS: " << solutions.size() << "\n\n";
  
#endif // !CGAL_BISOLVE_ARRANGEMENTS

  os.precision(7);

#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING 
#warning Resetting nr of pairs to 1 before timer output for consistency with external timeouts
#endif
  n_pairs = 1;
  os << "Resetting nr of pairs to 1 before timer output for consistency with external timeouts"
     << std::endl;

  os << "TIME res          : " << t_res[0].time()/n_pairs << std::endl;
  os << "TIME   res_x      :   " << t_res[1].time()/n_pairs << std::endl;
  os << "TIME   res_y      :   " << t_res[2].time()/n_pairs << std::endl;
  os << "TIME sqfree       : " << t_sqfree[0].time()/n_pairs << std::endl;
  os << "TIME   sqfree_x   :   " << t_sqfree[1].time()/n_pairs << std::endl;
  os << "TIME   sqfree_y   :   " << t_sqfree[2].time()/n_pairs << std::endl;
  os << "TIME solve1       : " << t_solve1[0].time()/n_pairs << std::endl;
  os << "TIME   solve1_x   :   " << t_solve1[1].time()/n_pairs << std::endl;
  os << "TIME   solve1_y   :   " << t_solve1[2].time()/n_pairs << std::endl;
  os << "TIME minprec      : " << t_minprec[0].time()/n_pairs << std::endl;
  os << "TIME   minprec_x  :    " << t_minprec[1].time()/n_pairs << std::endl;
  os << "TIME   minprec_y  :    " << t_minprec[2].time()/n_pairs << std::endl;
  os << "TIME wellsep      : " << t_wellsep[0].time()/n_pairs << std::endl;
  os << "TIME   wellsep_x  :   " << t_wellsep[1].time()/n_pairs << std::endl;
  os << "TIME   wellsep_y  :   " << t_wellsep[2].time()/n_pairs << std::endl;
  os << "TIME   approx     :   " << t_refine_app.time()/n_pairs << std::endl;
  os << "TIME   lensef     :   " << t_filter.time()/n_pairs << " (successfull: "
     << c_lensefilter_pass << "/" << c_t3over2_total << ")" << std::endl;
  os << "TIME   tshift     :   " << t_tshift.time()/n_pairs << std::endl;
  os << "TIME   ttest      :   " << t_ttest.time()/n_pairs << std::endl;
  os << "TIME sort         : " << t_sort[0].time()/n_pairs << std::endl;
  os << "TIME   sort_x     :   " << t_sort[1].time()/n_pairs << std::endl;
  os << "TIME   sort_y     :   " << t_sort[2].time()/n_pairs << std::endl;

  os << "TIME validate     : " << t_validate.time()/n_pairs << std::endl;
  os << "TIME vert_line_det: " << t_vl_detect.time()/n_pairs << std::endl;
  os << "TIME aiinit       : " << t_ai_init.time()/n_pairs << std::endl;
#if CGAL_BISOLVE_ENABLE_ARCAVOID
  os << "TIME   transform  : " << t_arca_transform.time()/n_pairs << std::endl;
#else
  os << "TIME   transform  : " << "n/a" << std::endl;
#endif
  os << "TIME aiprep       : " << t_ai_prep.time()/n_pairs << std::endl;
  os << "TIME   aisubdiv   :   " << t_ai_prep_subdiv.time()/n_pairs << std::endl;
  os << "TIME main_loop    : " << t_main_loop.time()/n_pairs << std::endl;
  os << "TIME   combinat   :   " << t_combinatorial.time()/n_pairs << std::endl;
#if CGAL_BISOLVE_ENABLE_ARCAVOID
  os << "TIME     arca_cert:   " << t_arca_cert.time()/n_pairs << " (successfull: "
     << c_arca_cert_success << "/" << c_arca_cert_total << ")" << std::endl;
#else
  os << "TIME     arca_cert: " << "n/a" << std::endl;
#endif
#if CGAL_BISOLVE_USE_RESULTANT_COFACTORS
  os << "TIME   predicates :   " << t_exin.time()/n_pairs << std::endl;
  os << "TIME     approx   :     " << t_approx.time()/n_pairs << " (" << c_approx/n_pairs << ",maxprec=" << g_max_prec << ")" << std::endl;
  os << "TIME     exia2    :     " << t_exia2.time()/n_pairs <<
    " (successfull: " << c_ieval2_reject/n_pairs << ")" << std::endl;
  os << "TIME     innorm   :     " << t_innorm.time()/n_pairs << std::endl;
  os << "TIME       subs   :       " << t_subs.time()/n_pairs << std::endl;
  os << "TIME       bounds :       " << t_bounds.time()/n_pairs << std::endl;
#else
  os << "TIME   predicates :   " << "n/a" << std::endl;
  os << "TIME     approx   :     " << "n/a" << std::endl;
  os << "TIME     exia2    :     " << "n/a" << std::endl;
  os << "TIME     innorm   :     " << "n/a" << std::endl;
  os << "TIME       subs   :       " << "n/a" << std::endl;
  os << "TIME       bounds :       " << "n/a" << std::endl;
#endif
  os << "TIME   aisubdiv   :   " << t_ai_main_subdiv.time()/n_pairs << std::endl;
//   os << "TIME timer: " << t_timer.time() << "\n";


  os << std::endl;
  os << "CA/Arrangement Events"<< std::endl;
  os << "TIME CEL           :   " << t_cel[1].time()/n_pairs << std::endl;
  os << "TIME   vl          :     " << t_cel_vl.time()/n_pairs << std::endl;
  os << "TIME   lcoeff      :     " << t_lcoeff.time()/n_pairs << std::endl;
  os << "TIME   mult_r_fx_fx:     " << t_mult_rxy.time()/n_pairs << std::endl;
  os << "TIME     res       :       " << t_mult_res.time()/n_pairs << std::endl;
  os << "TIME     gcd       :       " << t_mult_xy_gcd.time()/n_pairs << std::endl;
  os << "TIME     ntl-div   :       " << t_mult_ntl_div.time()/n_pairs << std::endl;
  os << "TIME     sqfs      :       " << t_mult_sqfs.time()/n_pairs << std::endl;
  os << "TIME     signs     :       " << t_mult_sign.time()/n_pairs << std::endl;
  os << "TIME   ffyi        :     " << t_ffy_isol[1].time()/n_pairs << std::endl;
#if CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
  os << "TIME     ffyi-ctor :       " << t_ffyi_tes_isol[0].time()/n_pairs << std::endl;
  os << "TIME     ffyi-tes  :       " << t_ffyi_tes_isol[1].time()/n_pairs << " (successfull: "
     << c_status_lines_tes[1] << "/" << c_status_lines_total[1] << ")" << std::endl;
#else
  os << "TIME     ffyi-tes  :       " << "n/a" << std::endl;
#endif
  os << "TIME     ffyi-ak1  :       " << t_ffyi_ak1_isol.time()/n_pairs << std::endl;
  os << "TIME     ffyi-bdv  :       " << t_ffyi_bdv_isol.time()/n_pairs << std::endl;
  os << "TIME       bdv     :         " << t_bdv.time()/n_pairs << std::endl;
  os << "TIME         ctor  :           " << t_bdv_construct.time()/n_pairs << std::endl;
  os << "TIME         b0    :           " << t_bdv_bs0.time()/n_pairs << std::endl;
  os << "TIME         common:           " << t_common.time()/n_pairs << std::endl;
  os << "TIME         higher:           " << t_bdv_bs.time()/n_pairs << std::endl;
  os << "TIME     t_ffy_refine:         " << t_ffy_refine.time()/n_pairs <<
    std::endl;
  os << "TIME        t_ffy_ais_subdivide:  " << t_ffy_ais_subdivide.time()/n_pairs << std::endl;
  
  os << "TIME   sllift      :     " << t_sllift[1].time()/n_pairs << std::endl;
  os << "TIME     init      :       " << t_sllift_init.time()/n_pairs << std::endl;
  os << "TIME     generic   :       " << t_sllift_generic.time()/n_pairs << std::endl;
  os << "TIME     bucket    :       " << t_sllift_bucket.time()/n_pairs << std::endl;
  os << "TIME     refine    :       " << t_sllift_refine.time()/n_pairs << std::endl;
  os << "TIME     assign    :       " << t_sllift_assign.time()/n_pairs << std::endl;
  os << "TIME     finish    :       " << t_sllift_finish.time()/n_pairs << std::endl;
  os << std::endl;

  os << "CA/Arrangement Non-Events"<< std::endl;
  os << "TIME CNEL          :   " << t_cel[2].time()/n_pairs << std::endl;
  os << "TIME   ffyi        :     " << t_ffy_isol[2].time()/n_pairs << std::endl;
#if CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER
  os << "TIME     ffyi-tes  :       " << t_ffyi_tes_isol[2].time()/n_pairs << " (successfull: "
     << c_status_lines_tes[2] << "/" << c_status_lines_total[2] << ")" << std::endl;
#else
  os << "TIME     ffyi-tes  :       " << "n/a" << std::endl;
#endif
  os << "TIME     ffyi-bits :       " << t_ffyi_bits_isol.time()/n_pairs << std::endl;
  os << "TIME   sllift      :     " << t_sllift[2].time()/n_pairs << std::endl;
  os << "TIME t_ais   : " << t_ais.time() / n_pairs << "\n";
  os << std::endl;

 
  os << std::endl;

  double t_project_wo_res = t_sqfree[0].time() + t_solve1[0].time() + t_minprec[0].time();
  
  double t_total_cpu = t_total.time();
  double t_project_cpu = t_res[0].time() + t_project_wo_res;

  double p_project_cpu = (t_project_cpu*100)/t_total_cpu;
  double p_seperate = (t_wellsep[0].time())*100/t_total_cpu;
  double p_validate = (t_validate.time()*100)/t_total_cpu;
  double p_123_cpu = p_project_cpu + p_seperate + p_validate;

  os << "TIME CPU-PROJECT  : " << t_project_cpu/n_pairs << " ( " << p_project_cpu << "% )" << std::endl;
  os << "TIME SEPERATE     : " << t_wellsep[0].time()/n_pairs << " ( " << p_seperate << "% )" << std::endl;
  os << "TIME VALIDATE     : " << t_validate.time()/n_pairs << " ( " << p_validate << "% )" << std::endl;
  
  double t_misc = t_total_cpu - t_project_cpu - t_wellsep[0].time() - t_validate.time();
  
  os << "TIME misc         : " << t_misc/n_pairs << " ( " << (100-p_123_cpu) << "% )" << std::endl;

  os << std::endl;

  os << "TIME total-cpu    : " << (t_total_cpu)/n_pairs << std::endl;
  os << "TIME total w/o res: " << (t_total.time() - t_res[0].time())/n_pairs << std::endl;

#if CGAL_BISOLVE_USE_AK2
  os << "\n###########################################################\n";
  os << "TIME cpa         :   " << t_cpa.time()/n_pairs << std::endl;
  os << "TIME msolve      :   " << t_msolve.time()/n_pairs << std::endl;
#endif

  os << std::endl;
 
  return 0;
}

int main(int argc, char **argv) {

#if CGAL_BISOLVE_USE_GMP
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  std::cout << " Test BISOLVE using Gmp " << std::endl;
  typedef CGAL::GMP_arithmetic_kernel AK;
#else
#error Gmp NOT AVAILABLE! Use other number type library
#endif
#endif
  
#if CGAL_BISOLVE_USE_CORE
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
  std::cout << " Test BISOLVE using Core " << std::endl;
  typedef CGAL::CORE_arithmetic_kernel AK;
#else
#error Core NOT AVAILABLE! Use other number type library
#endif
#endif

  typedef AK::Integer Coefficient;
  // TODO enhance entire code to also work with BigFloat as Bound 
  typedef AK::Rational Rational;

#if CGAL_BISOLVE_USE_RS_AK

  typedef Algebraic_kernel_rs_gmpz_d_1 Algebraic_kernel_d_1;

#else

#if CGAL_BISOLVE_USE_RS_ISOLATOR
  typedef CGAL::Algebraic_kernel_d_1_generator< Coefficient, Rational >
    ::Algebraic_kernel_with_qir_and_rs_1 Actual_algebraic_kernel_d_1;
#else
  typedef CGAL::Algebraic_kernel_d_1_generator< Coefficient, Rational >
    ::Algebraic_kernel_with_qir_and_bitstream_1 Actual_algebraic_kernel_d_1;
  // TODO WE SHOULD ACTUALLY USE THIS ALTERNATIVE, but it is slow
  // typedef CGAL::Algebraic_kernel_d_1< Coefficient, Rational >
  // Actual_algebraic_kernel_d_1;

  // TODO compile with CORE!
  //typedef CGAL::Algebraic_kernel_d_1_generator< Coefficient, Rational >
  //  ::Algebraic_kernel_with_qir_and_bitstream_1 Actual_algebraic_kernel_d_1;

  // TODO test this kernel?
  // typedef CGAL::Algebraic_kernel_d_1_generator< Coefficient, Rational >
  //  ::Algebraic_kernel_with_qir_and_arcavoid_1 Actual_algebraic_kernel_d_1;
#endif

#endif

  // TODO where to move the rounding?
  //! the univariate algebraic kernel
  typedef CGAL::internal::Rounding_ak_d_1< Actual_algebraic_kernel_d_1 > Algebraic_kernel_d_1;
  
  bisolve< Algebraic_kernel_d_1 >(argc, argv);
  


  return 0;
}
