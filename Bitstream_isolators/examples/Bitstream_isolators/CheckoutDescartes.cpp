#include <iostream>

#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_pure.h>
#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>


int main() {

	typedef CGAL::CORE_arithmetic_kernel::Integer Integer;
	typedef CGAL::internal::Bitstream_descartes<
            CGAL::internal::Bitstream_descartes_rndl_tree_traits
            <CGAL::internal::Bitstream_coefficient_kernel<Integer> > > Isolator;
	typedef CGAL::Polynomial_type_generator<Integer,1>::Type Polynomial;
  
	typedef CGAL::Polynomial_traits_d<Polynomial  > PT;
  
	PT::Construct_polynomial Construct_polynomial;

	CGAL::set_pretty_mode(std::cout);

	Integer coeff[4];
	coeff[0] = -6;
	coeff[1] = -7;
	coeff[2] = 0;
	coeff[3] = 1;
	Polynomial p = Construct_polynomial(coeff, coeff+4);
	std::cout << "Polynom: " << p << std::endl;
	Isolator isolator(p);
	std::cout << "Anzahl Nullstellen: " << isolator.number_of_real_roots() << std::endl;
	std::cout << "Erste Nullstelle: [" << isolator.left_bound(0) << "," << isolator.right_bound(0) << "]" << std::endl;
	std::cout << "Zweite Nullstelle: [" << isolator.left_bound(1) << "," << isolator.right_bound(1) << "]" << std::endl;
	std::cout << "Dritte Nullstelle: [" << isolator.left_bound(2) << "," << isolator.right_bound(2) << "]" << std::endl;
	
	return 0;

}
