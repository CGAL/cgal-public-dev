#include <iostream>

#include <CGAL/config.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Algebraic_kernel_1_generator.h>



int main() {

	CGAL::set_pretty_mode(std::cout);

	typedef CGAL::Polynomial_type_generator<double,1>::Type Poly;
  	typedef CGAL::Polynomial_traits_d<Poly> PT;
  	typedef PT::Innermost_coefficient_type ICoeff;
  	PT::Construct_polynomial Construct_polynomial;
  	PT::Degree Deg;
  	PT::Substitute Substitute;
    PT::Substitute_homogeneous hsubstitute;
  
	  //Polynom konstruieren mit Liste
  	std::list<ICoeff> a;
  	a.push_back(ICoeff(5.8));
  	a.push_back(ICoeff(1.5));
  	a.push_back(ICoeff(1.12));
  	//int arr[10];
  	//Poly p  = Construct_polynomial(arr, arr+10);

  	Poly p  = Construct_polynomial(a.begin(), a.end());
  	std::cout << "My Poly: " << p << std::endl;
  
  	//Polynom konstruieren mit Array
  	ICoeff b[3];
  	b[0] = ICoeff(2.5);
  	b[1] = ICoeff(3.5);
  	b[2] = ICoeff(1.12);
  	p = Construct_polynomial(b, b+3);
  	std::cout << "My second Poly: " << p << std::endl;
  
  	//Eintraege in Array schreiben  
  	double c [3];
  	for (int i = 0; i < 3; i++) {
    	c[i] = get_innermost_coefficient(p,i);
    	std::cout << "c[" << i << "] = " << c[i] << std::endl;
  	}
  
  	//Koeffizienten setzen
  	std::list<ICoeff> d;
  	for (int i = 0; i < 3; i++){
    	c[i] += 1;
    	d.push_back(ICoeff(c[i]));
  	}
  	p = Construct_polynomial(d.begin(), d.end());
  	std::cout << "My Shifted Poly second try: " << p << std::endl;

  	//Grad des Polynoms bestimmen
  	std::cout << "Degree: " << Deg(p) << std::endl;
  
  	std::cout << "p[0]: " << p[0] << std::endl;
  
  	std::cout << "p(1): " << p.evaluate(ICoeff(1)) << std::endl;
  
  	//scale_up und scale_down
  	double new_array[3];
  	new_array[0] = 1;
  	new_array[1] = 1;
  	new_array[2] = 1;
  	p = Construct_polynomial(new_array, new_array+3);
  	std::cout << "New Poly: " << p << std::endl;
  	Poly q = ::CGAL::scale_up(p,ICoeff(3));
  	std::cout << "scale_up(p,3): " << q << std::endl;
  	q = ::CGAL::scale_down(p,ICoeff(2));
  	std::cout << "scale_down(p,2): " << q << std::endl;
// 	p.scale_up(2);
// 	std::cout << "p.scale_up(2): " << p << std::endl;
  
	//translate
	q = ::CGAL::translate(p,ICoeff(2));
	std::cout << "translate(p,2): " << q << std::endl;
	q = ::CGAL::translate(p,ICoeff(3));
	std::cout << "translate(p,3): " << q << std::endl;
	
	//reversal
	q = reversal(q);
	std::cout << "reversal(q): " << q << std::endl;
    
    //evaluate	
	std::cout << "p.evaluate(2): " << p.evaluate(2) << std::endl;

    //substitute homogeneous
    ICoeff x_coeff[2];
    x_coeff[0] = ICoeff(3);
    x_coeff[1] = ICoeff(-2); 
    Poly x;
    x = Construct_polynomial(x_coeff, x_coeff+2);
    ICoeff y_coeff[2];
    y_coeff[0] = ICoeff(1);
    y_coeff[1] = ICoeff(1);
    Poly y;
    y = Construct_polynomial(y_coeff, y_coeff+2);
    std::list<Poly> replacements;
    replacements.push_back(x);
    replacements.push_back(y);
    q = hsubstitute(p, replacements.begin(), replacements.end());
    std::cout << "p: " << p << std::endl;
    std::cout << "x: " << x << std::endl;
    std::cout << "y: " << y << std::endl;
    std::cout << "hsubstitute(p,-2x+3,x+1): " << q << std::endl;

    //decompose
    typedef CGAL::CORE_arithmetic_kernel Arithmetic_kernel;
    typedef Arithmetic_kernel::Rational Rational;
    typedef CGAL::Fraction_traits<Rational> FT;
    FT::Decompose decompose;

    Rational test = Rational(1,5);
    FT::Numerator_type numerator;
    FT::Denominator_type denominator;

    decompose(test,numerator,denominator);
    std::cout << "numerator: " << numerator << std::endl;
    std::cout << "denominator: " << denominator << std::endl;

    //derivative
    std::cout << "p: " << p << std::endl;
    p = ::CGAL::diff(p);
    std::cout << "p': " << p << std::endl;
	
  	return 0;

}

