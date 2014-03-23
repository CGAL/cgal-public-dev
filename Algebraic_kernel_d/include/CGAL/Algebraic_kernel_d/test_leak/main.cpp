#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_d/flags.h>

// Switches on/off tests for Sqrt-extension types
#if CGAL_ACK_WITH_ROTATIONS
#ifndef DO_SQRT_EXTENSION_TESTS
#define DO_SQRT_EXTENSION_TESTS 1
#endif
#endif

#include <CGAL/basic.h>

#include <sstream>

#if CGAL_ACK_USE_EXACUS
#include <AcX/Algebraic_curve_2.h>
#include <AcX/Algebraic_curve_pair_2.h>
#endif

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h>
#include <CGAL/Algebraic_kernel_d/Curve_analysis_2.h>
#include "Curve_pair_analysis_2.h"


typedef CGAL::CORE_arithmetic_kernel Arithmetic_kernel;
typedef Arithmetic_kernel::Rational Rational;
typedef Arithmetic_kernel::Integer Integer;

typedef Integer Coefficient;
typedef CGAL::Polynomial_type_generator<Coefficient,1>::Type Poly_1;
typedef CGAL::Polynomial_type_generator<Coefficient,2>::Type Poly_2;



template<typename Poly_> Poly_ from_string(const char* s) {
    std::stringstream ss(s);
    Poly_ f;
    ss >> f;
    return f;
}
int main()
{
        Poly_2 f=from_string<Poly_2>("P[4(0,P[4(3,-1)(4,2)])(2,P[1(1,1)])(4,P[0(0,1)])]");
        Poly_2 g=from_string<Poly_2>("P[4(0,P[4(4,1)])(1,P[2(2,1)])(3,P[0(0,-1)])(4,P[0(0,2)])]");
        Curve_analysis_2 c1=construct_curve_2(f);
        Curve_analysis_2 c2=construct_curve_2(g);
        Curve_pair_analysis_2 curve_pair=construct_curve_pair_2(c1,c2);

	return 0;
}
