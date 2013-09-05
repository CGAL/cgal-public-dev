
namespace CGAL {

/*!
\ingroup PkgPolynomialClasses

The class `Polynomial_parser_d` provides functionality to read polynomials in a human-readable format from the input stream. In this sense, it can be regarded as an extension 
of `operator >>` from `Polynomial<Coeff>`. However, it should be noted that the `Polynomial_parser_d` is <B>not</B> a replacement of the latter input operator. Instead,
it is provided as an independent feature. 

The polynomial parser can parse algebraic expressions with `+`, `-`, `*` and `^` (raising to integer power) operations and arbitrary number of parentheses. Additionally, it is allowed to have a single equality operator `=` 
in an expression, such that \f$ f(x,y,\dots)=g(x,y,\dots)\f$ will be  
parsed as \f$ f(x,y,\dots)-g(x,y,\dots).\f$ 

We define a context-free grammar of the parser as follows. 
Non-terminating symbols:
<UL>
<LI> `E` - complete algebraic expression (starting symbol)
<LI> `sE` - algebraic subexpression 
<LI> `T` - single expression term
<LI> `F` - single term factor
</UL>
Terminating symbols:
<UL>
<LI> `var` - valid variable name (depends on the actual parser policy)
<LI> `C` - numeric constant (format depends on the actual parser policy)
<LI> `n` - non-negative integer constant (exponent)
<LI> <b>+</b>, <b>-</b>, <b>*</b>, <b>^</b> - binary operations
<LI> <b>=</b> - comparison operator
<LI> `(`, `)` - parantheses
</UL>

Parsing rules written in Backus Normal Form are given below: 
\f[
\begin{array}{rcl}
E  & ::= & sE\ |\ sE = sE\\
sE & ::= & T\  |\ sE + T\ |\ sE - T\\
T  & ::= & F\  |\ F^n |\ T * T\ |\ TT \\
F  & ::= & (E)\ |\ C\ |\ \mbox{var}
\end{array}
\f]
Remark that, the parser's grammar is not strict about the multiplication symbol (<b>*</b>)
which can be omitted. For instance, an expression of the form
\f[
34(x+yyyy)(x+z)^2*z
\f]
will be parsed as:
\f[
34*(x+y^4)*(x+z)^2*z
\f]
However, such ambiguity decreases readability of exressions, therefore we 
would always recommend to use the multiplication symbol even though it is optional.

As mentioned in the manual pages, the parser can deal with <i>any</i> sorts of 
space characters including `\f`, `\n`, `\r`, `\t`, `\v`. 
Although, an expression to be parsed can be multiline, the comments
inside the expression are <b>not</b> supported.

The parser's behavior is fully customizeable via the so-called <i>parser policy</i>
provided as a template parameter. 
In the first place, the policy is responsible for reading polynomial 
coefficients from the input stream and performing all necessary type conversions
for the output polynomial. In such a way, we "detach" the actual parsing algorithm
from type conversion routines, validity checks, etc. 
In particular, the parser policy is used to name polynomial variables,
perform maximal degree check, and detect various erroneous situations (such as division by zero). Additionally, with the help of the parser policy, one can mix integer, rational and floating-point numbers in a <i>single</i> algebraic expression 
or round coefficients to a certain precision since 
postprocessing of the polynomial coefficients is performed by the parser policy "on-the-fly" (for examples see the corresponding part of the CGAL users manual). 

Any user-defined policy must be a model of `ParserPolicy` concept.
If no parser policy is given to `Polynomial_parser_d`, the parser will use
`Default_parser_policy`.

*/
template< typename Polynomial_d, typename ParserPolicy >
class Polynomial_parser_d {
public:

/// \name Constants 
/// @{

/*!
The number of variables in the input polynomial (dimension)
*/ 
static const int d; 

/// @} 
    
/// \name Types 
/// @{

/*!
The type of the resulting polynomial to be constructed: must a model of `Polynomial_d`
*/ 
typedef unspecified_type Polynomial_d; 

/*!
  Parser policy type: must be a model of `ParserPolicy`
*/ 
typedef unspecified_type Policy;

/// @} 

/// \name Creation 
/// @{

/*!
Standard constructor: takes the parser policy object as an optional parameter
*/ 
Polynomial_parser_d(const Policy& policy = Policy());
    
/// @} 

/// \name Operations
/// @{    
    
/*!
Parses the input string <i>in</i> and constructs the resulting polynomial <i>poly</i>. Returns <b>true</b> if parsing was successful.
*/ 
bool operator()(const std::string& in, Polynomial_d& poly) const;
    
/// @}    

}; /* end Polynomial_parser_d */
} /* end namespace CGAL */
