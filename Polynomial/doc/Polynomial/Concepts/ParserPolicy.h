
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

Parser policy concept defines the interface that any custom policy must provide.

\sa `CGAL::Polynomial<Coeff>`
\sa `CGAL::Polynomial_parser_d<Polynomial_d,ParserPolicy>`

\cgalHasModel `Default_parser_policy`
\cgalHasModel `Mixed_rational_parser_policy`
\cgalHasModel `Mixed_floating_point_parser_policy`

*/

class ParserPolicy {
public:

/// \name Types 
/// @{

/*!
  The type of the output polynomial to be constructed: must be a model of `Polynomial_d`
*/ 
typedef unspecified_type Polynomial_d; 

/*!
`Innermost_coefficient_type` of Polynomial_d
*/ 
typedef unspecified_type Innermost_coefficient_type; 

/// @}
/// \name Creation 
/// @{

/*!
Default constructor 
*/ 
ParserPolicy();
    
/// @}
/// \name Operations 
/// @{

/*!
Reads a polynomial coefficient from the input stream <i>is</i> and converts it to the output type `Innermost_coefficient_type`. <i>is</i> points to the first digit character of this coefficient. Throws an exception of type `Parser_exception` in case of error.
*/ 
Innermost_coefficient_type read_coefficient(std::istream& is) const;
    
/*!
Performs exponent check on the polynomial coefficients. 

The function returns <b>true</b> if the degree <i>deg</i> of a monomial being parsed lies within the accepted range and <b>false</b> otherwise. If no exponent check is required, the function should always return <b>true</b>.
*/ 
bool exponent_check(unsigned deg) const;

/*!
The function checks if character <i>ch</i> is a valid variable name and, if so, returns its index.

If <i>ch</i> is a not a valid polynomial variable name, the function returns <b>false</b>. Otherwise, the function computes a 0-based index of this polynomial variable (where 0 corresponds to the innermost variable) and returns it in <i>idx</i> reference. In other words, if <i>idx</i> is set to 0, then <i>ch</i> is the innermost variable of a polynomial. Otherwise, if <i>idx</i> is <i>d-1</i>, where <i>d</i> is the number of polynomial variables, then <i>ch</i> is the outermost variable. The function returns <b>true</b> if <i>ch</i> is a valid variable name.
*/ 
bool check_var_name(char ch, int& idx) const;

/// @}

}; /* end ParserPolicy */

