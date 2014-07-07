#ifndef POLYNOMIAL_PARSER_H
#define POLYNOMIAL_PARSER_H
#include <iostream>
#include <string>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/bind.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <CGAL/Polynomial.h>

/**
Defines the static function Parse that turns user-input strings into CGAL
Polynomials.

When using this header, make sure to include this *before* any Qt GUI-related
headers. This is because the particular Boost libraries used here may
conflict with the auto-generated symbols generated elsewhere by Qt.
*/
template < class TPolynomial >
class PolynomialParser
{
  typedef TPolynomial Polynomial_2;

  // Internal struct used for parsing polynomials.
  struct PolynomialTerm
  {
    boost::optional<int> coefficient;
    boost::optional<int> x_exponent;
    boost::optional<int> y_exponent;
  };

  /**
  \param[in] begin and end specify the range of a collection of terms.
  The value type of Iterator is client::PolynomialTerm.
  */
  template < typename Iterator >
  static Polynomial_2 MakePolynomial( Iterator begin, Iterator end )
  {
    Polynomial_2 x = CGAL::shift( Polynomial_2( 1 ), 1, 0 );
    Polynomial_2 y = CGAL::shift( Polynomial_2( 1 ), 1, 1 );
    std::vector< Polynomial_2 > terms;
    for ( Iterator it = begin; it != end; ++it )
    {
      PolynomialTerm& tt = *it;
      int coefficient = (tt.coefficient)? *tt.coefficient : 1;
      int x_exponent = (tt.x_exponent)? *tt.x_exponent : 0;
      int y_exponent = (tt.y_exponent)? *tt.y_exponent : 0;
      Polynomial_2 term =
        coefficient * CGAL::ipower( x, x_exponent ) * CGAL::ipower( y, y_exponent );
      terms.push_back( term );
    }


    Polynomial_2 res;
    if ( ! terms.size( ) )
      return res;

    res = terms[0];
    for ( int i = 1; i < terms.size( ); ++i )
    {
      res = res + terms[i];
    }
    return res;
  }

  /**
  Defines the grammar for parsing bivariate polynomials supported in the
  Algebraic segment traits in the Arrangement_2 package.
  */
  template <typename Iterator>
  struct PolynomialGrammar : boost::spirit::qi::grammar<Iterator,
      std::vector<PolynomialTerm>(),
      boost::spirit::ascii::space_type>
  {
    typedef boost::spirit::ascii::space_type space_type;

    // parser rules
    boost::spirit::qi::rule<Iterator, int(), space_type> exponent;
    boost::spirit::qi::rule<Iterator, int(), space_type> x_exponent;
    boost::spirit::qi::rule<Iterator, int(), space_type> y_exponent;
    boost::spirit::qi::rule<Iterator, PolynomialTerm(), space_type> poly_term;
    boost::spirit::qi::rule<Iterator, PolynomialTerm(), space_type> negative_poly_term;
    boost::spirit::qi::rule<Iterator, std::vector<PolynomialTerm>(), space_type> start;

    PolynomialGrammar() : PolynomialGrammar::base_type(start)
    {
      namespace phx = boost::phoenix;
      namespace qi = boost::spirit::qi;
      namespace ascii = boost::spirit::ascii;

      using qi::int_;
      using qi::lit;
      using qi::double_;
      using qi::lexeme;
      using ascii::char_;
      using qi::_val;
      using qi::eps;

      exponent %= '^' >> int_;

      x_exponent %= 'x' >> ( exponent | eps[_val = 1] );

      y_exponent %= 'y' >> ( exponent | eps[_val = 1] );

      poly_term = eps[_val = PolynomialTerm()]
        >>
          -int_[phx::bind(&PolynomialTerm::coefficient, _val) = qi::_1]
        >>
          -x_exponent[phx::bind(&PolynomialTerm::x_exponent, _val) = qi::_1]
        >>
          -y_exponent[phx::bind(&PolynomialTerm::y_exponent, _val) = qi::_1]
      ;

      negative_poly_term = eps[_val = PolynomialTerm()]
        >> (
            int_[phx::bind(&PolynomialTerm::coefficient, _val) = -1 * qi::_1]
          |
            eps[phx::bind(&PolynomialTerm::coefficient, _val) = -1]
        )
        >>
          -x_exponent[phx::bind(&PolynomialTerm::x_exponent, _val) = qi::_1]
        >>
          -y_exponent[phx::bind(&PolynomialTerm::y_exponent, _val) = qi::_1]
      ;

      start = eps[_val = std::vector<PolynomialTerm>()]
        >>
          poly_term[phx::push_back(_val, qi::_1)]
        >>
          *(
              ('+' >> poly_term[phx::push_back(_val, qi::_1)])
            |
              ('-' >> negative_poly_term[phx::push_back(_val, qi::_1)])
          )
      ;
    }
  };

public:
  /**
  Parse the user-input string into a polynomial.
  \param[in] str the user input
  \param[out] res pointer to the polynomial output variable
  \return whether \a str was successfully parsed
  */
  static bool Parse( const std::string& str, Polynomial_2* res )
  {
    using boost::spirit::ascii::space;
    typedef std::string::const_iterator iterator_type;
    typedef PolynomialGrammar<iterator_type> ParserType;

    ParserType parser;
    std::vector<PolynomialTerm> poly;
    bool ok = phrase_parse(str.begin(), str.end(), parser, space, poly);
    if ( ok )
    {
      *res = MakePolynomial( poly.begin( ), poly.end( ) );
    }
    return ok;
  }
};
#endif // POLYNOMIAL_PARSER_H
