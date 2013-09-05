
namespace CGAL {
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

The class `Default_parser_policy` realizes `ParserPolicy` without explicit type conversions. In other words, this policy reads poynomial coefficients of type `Polynomial_d::Innermost_coefficient_type` directly from the input stream
(provided that the `operator >>` is available for `Polynomial_d::Innermost_coefficient_type`
) and constructs the resulting polynomial of type `Polynomial_d`.


\sa `CGAL::Polynomial_parser_d<Polynomial_d,ParserPolicy>`
\sa `Mixed_rational_parser_policy`
\sa `Mixed_floating_point_parser_policy`

\cgalModels `ParserPolicy`

*/

template < typename Polynomial_d >
class Default_parser_policy {
public:

/// \name Types 
/// @{

/*!
  The type of the output polynomial to be constructed: must be a model of `Polynomial_d`
*/ 
typedef unspecified_type Polynomial_d; 

/*!
 Innermost coefficient type of Polynomial_d
*/ 
typedef typename  CGAL::Polynomial_traits_d< Polynomial_d >::
        Innermost_coefficient_type Innermost_coefficient_type; 

/// @}
/// \name Creation 
/// @{

/*!
Standard constructor. 
Optionally takes a string of variable names <i>var_names</i>, where the first character is the name of the innermost polynomial variable, and the last one 
is the name of the outermost polynomial variable. If <i>var_names</i> is 
<b>NULL</b>, the default string of variable names 
\f$\mathsf{xyzwabcdfghijklmnopqrstuv}\f$ will be used. The character \f$\mathsf{e}\f$ 
is <b>not</b> used by default since it can be confused with a floating-point exponent.
*/ 
Default_parser_policy(const char *var_names = NULL);
    
/// @}
/// \name Operations 
/// @{
    
/*!
Reads a polynomial coefficient from the input stream <i>is</i>.

See `ParserPolicy` for details.
*/ 
Innermost_coefficient_type read_coefficient(std::istream& is) const;

/*!
Checks if character <i>ch</i> is a valid variable name and, if so, returns its index.

See `ParserPolicy` for details.
*/ 
bool check_var_name(char ch, int& idx) const;

/*!
Performs exponent check on the polynomial coefficients. 

See `ParserPolicy` for details. This function always 
returns <b>true</b> since `Default_parser_policy` does not perform exponent check.
*/ 
bool exponent_check(unsigned deg) const;
    
/// @}    

}; /* end Default_parser_policy */

} /* end namespace CGAL */