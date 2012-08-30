#ifndef CGAL_BI_SOLVE_2_FLAGS_H
#define CGAL_BI_SOLVE_2_FLAGS_H 1

// bisolve

#ifndef CGAL_BISOLVE_DISABLE_LOCALIZED_OPERATOR
#define CGAL_BISOLVE_DISABLE_LOCALIZED_OPERATOR 0 // 0 is default
#endif
#if CGAL_BISOLVE_DISABLE_LOCALIZED_OPERATOR
#warning Disabled localized operator for x
#endif

#ifndef CGAL_BISOLVE_ENABLE_LENSE_FILTER_IN_T_TEST
#define CGAL_BISOLVE_ENABLE_LENSE_FILTER_IN_T_TEST 0 // 0 is default
#endif
#if CGAL_BISOLVE_ENABLE_LENSE_FILTER_IN_T_TEST
#warning Lense filter in well-seperation is enabled
#endif
 
// certifier

#ifndef CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION
#define CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION 0 // 0 is default
#endif
#if CGAL_BISOLVE_DISABLE_COMBINATORIAL_CERTIFICATION
#warning Combinatorial certification is disabled
#endif

#ifndef CGAL_BISOLVE_DISABLE_BIDIRECTIONAL_CERTIFICATION
#define CGAL_BISOLVE_DISABLE_BIDIRECTIONAL_CERTIFICATION 0 // 0 is default
#endif
#if CGAL_BISOLVE_DISABLE_BIDIRECTIONAL_CERTIFICATION
#warning Certification in both directions is disabled
#endif

#ifndef CGAL_BISOLVE_ENABLE_ARCAVOID 
#define CGAL_BISOLVE_ENABLE_ARCAVOID 1 // default TODO?
#endif
#if CGAL_BISOLVE_ENABLE_ARCAVOID 
#undef CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER
#define CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER 1 // DISABLE BITSTREAM -- DO NOT CHANGE
#warning Enabled Arcavoid numerical solver
#endif

#ifndef CGAL_BISOLVE_ENABLE_NTL_FACTORIZE
#define CGAL_BISOLVE_ENABLE_NTL_FACTORIZE 0 // 0 is default
#endif
#if CGAL_BISOLVE_ENABLE_NTL_FACTORIZE
#warning NTL Factorize is enabled
#endif


#ifndef CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER
#define CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER 0 // 0 is default
#endif
#if CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER && !CGAL_BISOLVE_ENABLE_ARCAVOID
#warning Bitstreamfilter is disabled
#endif

#ifndef CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
#define CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS (CGAL_BISOLVE_DISABLE_BITSTREAM_FILTER && !CGAL_BISOLVE_ENABLE_ARCAVOID)
#endif
#if CGAL_BISOLVE_DISABLE_ACTIVE_INTERVALS
#warning Active intervals disabled
#endif

#ifndef CGAL_ACK_BITSTREAM_USES_E08_TREE
#define CGAL_ACK_BITSTREAM_USES_E08_TREE 1 // do not change
#endif

// set precision according to the results of the norm test
#ifndef CGAL_BISOLVE_USE_ADJUSTABLE_PRECISION
#define CGAL_BISOLVE_USE_ADJUSTABLE_PRECISION 1 // 1 is default
#endif

// verbosity

#ifndef CGAL_BISOLVE_DEBUG
#define CGAL_BISOLVE_DEBUG 0 // 0 is default
#endif

#ifndef CGAL_BISOLVE_VERBOSE
#define CGAL_BISOLVE_VERBOSE 0 // 0 is default
#endif

#if CGAL_BISOLVE_VERBOSE
#define Bisolve_out(x) std::cout << x;
#define dbl(x) CGAL::to_double(x)
#define bfi(x) CGAL::lower(CGAL::convert_to_bfi(x))
#define STILL_ALIVE std::clog << __LINE__ << "\n";
#else
#define Bisolve_out(x) static_cast< void >(0);
#define STILL_ALIVE static_cast< void >(0);
#endif

// telemetry

#ifndef CGAL_BISOLVE_TELEMETRY
#define CGAL_BISOLVE_TELEMETRY 0 // 0 is default
#endif

#if CGAL_BISOLVE_TELEMETRY
#define Bisolve_telemetry_code(x) x;
#else
#define Bisolve_telemetry_code(x) 
#endif


#endif // CGAL_BI_SOLVE_2_FLAGS_H
