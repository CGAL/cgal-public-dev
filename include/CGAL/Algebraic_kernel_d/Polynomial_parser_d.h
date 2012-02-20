// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : ??
// File          : ../Polynomial_parser_d.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*! \file Polynomial_parser_d.h
 *  \brief Parser for d-variate polynomials
 *  
 *  Defines a parser for d-variate polynomials in Maple format
 */

#ifndef CGAL_POLYNOMIAL_PARSER_D_H
#define CGAL_POLYNOMIAL_PARSER_D_H

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Bigfloat_interval_traits.h>
// #include <boost/type_traits/is_same.hpp>

#include <vector>
#include <sstream>

namespace CGAL {

namespace internal {

class Parser_exception {

public:
    Parser_exception(std::string msg) :
        _m_error(msg) {
    }
        
    std::string get_message() const
    { return _m_error; }

private:
    std::string _m_error;
 
};

} // namespace internal

//! default policy: reads in a coefficient of type \c InputCoeff from the input
//! stream and converts it to \c Poly_d_ coefficients type using provided
//! type coercion
template < class Poly_d_, class InputCoeff =
        typename CGAL::Polynomial_traits_d< Poly_d_ >::
            Innermost_coefficient_type >
struct Default_parser_policy {

    //! first template argument type
    typedef Poly_d_ Poly_d;
    //! second template argument type
    typedef InputCoeff Input_coeff;
    //! coefficient type
    typedef typename  CGAL::Polynomial_traits_d< Poly_d >::
        Innermost_coefficient_type Coeff;

    enum CoeffType {
        COEFF_INTEGER,
        COEFF_RATIONAL,
        COEFF_FLOAT
    };

    //! default constructor
    Default_parser_policy() {
    }
    
    //! reads in the coefficient from the input stream and converts it
    //! to \c Coeff type using provided type coercion
    //! \c is points to the first digit character of the coefficient
    Coeff read_coefficient(std::istream& is) const {
        return read_coeff_proxy(is, COEFF_INTEGER);
    }

    //! checking for degree overflow: can be used in real-time applications
    virtual bool exponent_check(unsigned) const {
        return true;
    }

    virtual Coeff read_coeff_proxy(std::istream& is, CoeffType) const {
        return _read_coeff< Input_coeff >(is);
    }

    virtual ~Default_parser_policy() {
    }

    //! variable names listed in the order from innermost to the outermost
    //! variable as they appear in the resulting equation
    static const int n_var_names = 4;    
    static const char *var_names_lower;
    static const char *var_names_upper;

protected:

    template < class NT >
    Coeff _read_coeff(std::istream& is) const {
        NT cf;
        is >> CGAL::iformat(cf);
        typename CGAL::Coercion_traits< NT, Coeff >::Cast cvt;
        return cvt(cf);
    }

};

template < class Poly_d_, class InputCoeff >
const char * Default_parser_policy< Poly_d_, InputCoeff >::
        var_names_lower = "xyzw";

template < class Poly_d_, class InputCoeff >
const char * Default_parser_policy< Poly_d_, InputCoeff >::
        var_names_upper = "XYZW";

//! This parser policy allows to mix integer and rational coefficients in a
//! single equation. Appropriate type coercion with the \c Poly_d_ coefficient
//! type must be provided
template < class Poly_d_ >
struct Mixed_rational_parser_policy :
        public Default_parser_policy< Poly_d_ > {

    //! template argument type
    typedef Poly_d_ Poly_d;
    //! base class
    typedef Default_parser_policy< Poly_d > Base;
    //! type of polynomial coefficient
    typedef typename Base::Coeff Coeff;

    typedef typename CGAL::Get_arithmetic_kernel< Coeff >::
             Arithmetic_kernel AK;
    //! integer number type
    typedef typename AK::Integer Integer;
    //! rational number type
    typedef typename AK::Rational Rational;
    //! input coefficient types
    typedef typename Base::CoeffType CoeffType;

    //! default constructor
    Mixed_rational_parser_policy() {
    }

    //! reads in a rational or integer coefficient from the input stream
    //! depending on whether '/' has been met during parsing and converts
    //! it to the \c Coeff type using provided coercion
    //! \c is points to the first digit character of the coefficient
    //! throws \c Parser_exception if error occurred
    Coeff read_coefficient(std::istream& is) const {

        std::stringstream buf;

        bool read_rational; // recognised as rational coeff
        if(!_eat_characters(is, buf, read_rational))
            throw internal::Parser_exception("Wrong coefficient type");

        return read_coeff_proxy(buf, read_rational ?
            Base::COEFF_RATIONAL : Base::COEFF_INTEGER);
    }

    virtual Coeff read_coeff_proxy(std::istream& is, CoeffType type)
            const {

        if(type == Base::COEFF_INTEGER)
            return Base::template _read_coeff< Integer >(is);
        else 
            return Base::template _read_coeff< Rational >(is);
    }

protected:
    //! protected stuff

    bool _eat_characters(std::istream& is, std::ostream& buf,
            bool& read_rational) const {

        buf.put(is.get()); // the first character must be a digit
        read_rational = false; // recognised as rational coeff
        while(!is.eof()) {
            char ch = is.peek();
               
            if(isdigit(ch)) {
                buf.put(is.get());
            } else if(ch == '/') { // found separator

                if(read_rational)
                    return false;
                read_rational = true;
                buf.put(is.get());
                
                ch = is.peek();
                // character immediately after '/' must be a digit
                if(!isdigit(ch))
                    return false;
                buf.put(is.get());
            } else
                break;
        }
        return true;
    }
};

//! This parser policy allows to mix integer, rational and floating-point
//! coefficients in a single equation. Appropriate type coercion with the
//! \c Poly_d_ coefficient type must be provided
//! \c FP_rounding optional parameter can be used to round the floating-point
//! number to desired precision before returning it to the parser
template < class Poly_d_/*, class FP_rounding = CGAL::Identity<KeyType_>*/ >
struct Mixed_floating_point_parser_policy :
        public Mixed_rational_parser_policy< Poly_d_ > {

    //! template argument type
    typedef Poly_d_ Poly_d;
    //! base class
    typedef Mixed_rational_parser_policy< Poly_d > Base;
    //! type of polynomial coefficient
    typedef typename Base::Coeff Coeff;

    typedef typename CGAL::Get_arithmetic_kernel< Coeff >::
             Arithmetic_kernel AK;
    //! integer number type
    typedef typename AK::Integer Integer;
    //! rational number type
    typedef typename AK::Rational Rational;
    //! bigfloat number type
    typedef typename CGAL::Bigfloat_interval_traits<
            typename AK::Bigfloat_interval >::Bound BigFloat;
    //! input coefficient types
    typedef typename Base::CoeffType CoeffType;

    //! default constructor
    Mixed_floating_point_parser_policy() {
    }

    //! allows mixing of integer, rational and floating-point coefficients
    //! in a single equation
    //! \c is points to the first digit character of the coefficient
    //! throws \c Parser_exception if an error occurred
    Coeff read_coefficient(std::istream& is) const {

        std::stringstream buf;
        bool read_rational; // recognised as rational coeff
        // first try to find rational coeffs, so we do not bother ourselves
        // with them later 
        if(!Base::_eat_characters(is, buf, read_rational))
            throw internal::Parser_exception("Wrong coefficient type");

        if(read_rational) { // successfully read rational coeff => done
            return read_coeff_proxy(buf, Base::COEFF_RATIONAL);
        }

        // now we either hit a non-digit character, found the exponent
        // or fraction delimiter '.'; assume the first option
        int fp_type = -1; // -1 - integer, 0 - float, 1 - float scientific,

        while(!is.eof()) {
            char ch = is.peek();

            bool parse_exp = (ch == 'e' || ch == 'E');

// printf("ch = %c; fp_type: %d; parse_exp: %d: buf: %s\n", ch,
//                     fp_type, parse_exp, buf.str().c_str());

            if(ch == '.' || parse_exp) {

                if(parse_exp && (fp_type == -1 || fp_type == 0))
                   fp_type = 1; // found exponent => scientific format

               else if(fp_type != -1) // all other combinations: error
                   throw internal::Parser_exception("Wrong floating-point coefficient");

               buf.put(is.get()); // eat comma
               if(parse_exp) {
                   ch = is.peek();
                   if(ch == '-' || ch == '+') {
                       buf.put(is.get());
                       ch = is.peek();
                   }
                   if(!isdigit(ch))
                       throw internal::Parser_exception("Wrong fp exponent");
               }

               if(fp_type == -1)
                   fp_type = 0; // recognised as floating-point

            } else if(isdigit(ch)) {
                buf.put(is.get()); // eat this
            } else 
                break;
        }
        return read_coeff_proxy(buf, (fp_type == -1 ? Base::COEFF_INTEGER :
                Base::COEFF_FLOAT));
    }

    virtual Coeff read_coeff_proxy(std::istream& is, CoeffType type)
            const {

        if(type == Base::COEFF_INTEGER)
            return Base::template _read_coeff< Integer >(is);
        else if(type == Base::COEFF_RATIONAL)
            return Base::template _read_coeff< Rational >(is);
        else // floating-point
            return Base::template _read_coeff< BigFloat >(is);
    }

};

/*! \brief parsing d-variate polynomials in Maple format
 * 
 * input format (y is outermost variable), e.g.:
 * (y-1)^4 + (-1)*y^3 + (x + z)^2 - (2123234523*x^2 - 2*y*y*x + 3*x*132123)^3
 * (y + x - z + w - 3)^3 = x + y - 123/12312 + 1.00001z*x^2
 */
template <class Poly_d_, class ParserPolicy =
                    Default_parser_policy< Poly_d_ > >
struct Polynomial_parser_d
{
    //!\name public typedefs
    //!@{

    //! this instance's template argument
    typedef Poly_d_ Poly_d;
    //! this instance's second template argument
    typedef ParserPolicy Policy;
    
    //! polynomial policy
    typedef CGAL::Polynomial_traits_d< Poly_d > PT;

    //! polynomial coefficient type
    typedef typename Policy::Coeff NT;
    //! multivariate dimension
    static const int d = CGAL::internal::Dimension< Poly_d >::value;

    //!@}
public: 
    /// \name Public methods
    //!@{

    //! default constructor
    Polynomial_parser_d() {
    }

    //! \brief functor invokation operator
    bool operator()(const std::string& in, Poly_d& poly) {
        
        try {
            // remove white spaces from the string: look for all possible
            // whitespace characters
            std::string s(in);
            const char *whitespace = " \f\n\r\t\v";
            std::size_t cnt = s.find_first_of(whitespace,0);
            while(cnt != std::string::npos){
                s.erase(cnt, 1);
                cnt = s.find_first_of(whitespace, cnt);
            }
                        
            //to handle the case when there is one '=' sign
            //we compute sides separetely and then subtract one from another
            std::size_t loc;
            std::string rhs;
            if((loc = s.find('=', 0)) != std::string::npos) {
                rhs = s.substr(loc + 1); // right-hand side
                s.erase(loc);
            }
            poly = get_poly(s);

            if(loc != std::string::npos) {
                Poly_d tmp = get_poly(rhs);
                poly -= tmp;
            }
        } 
        catch(internal::Parser_exception ex) {
            std::cerr << "Parser error: " << ex.get_message() << std::endl;
            return false;
        }
//         catch(...) {
//             std::cerr << "Unhandled exception during parsing" << std::endl;
//             return false;
//         }
        return true;
    }
    
    //!@}
protected: 
    /// \name protected methods
    //!@{ 
    
    //! given a string \c cstr of length \c len starting with an open
    //! parentheses returns the place marking the end of the corresponding
    //! closing parentheses
    std::size_t match_parenth(std::istream& is)
    {   
        int count = 0;
        std::size_t pos = 0, start = is.tellg(); // save pos of the first '('
        do {
            if(is.eof())
            return -1; // illegal number of parentheses
            char ch = is.peek();
            if(ch == '(')
                count++;
            else if(ch == ')') 
                count--;
            is.seekg(1, std::ios_base::cur);
            pos++;
        } while(count != 0); //j is one more than the matching ')'
        is.seekg(start);
        // pos is the number of chars including ( and )
        return pos - 1u; // return position of the last ')'
    }

    //! constructs {x/y}^exp and returns the result as bivariate polynomial 
    //! \c res \c var encodes a variable (x or y) 
    Poly_d construct_monomial(int idx, unsigned int exp) {
        static Poly_d one(NT(1));
        if(exp == 0) // just return 1
            return one;
        typename PT::Shift shift;
        Poly_d monom = shift(one, exp, idx);
        return monom;
    }

    bool check_var_names(char ch, int& idx) {

        int i;
        for(i = 0; i < Policy::n_var_names; i++) {
            if(Policy::var_names_lower[i] == ch ||
                    Policy::var_names_upper[i] == ch)
                break;
        }
        if(i == Policy::n_var_names || i >= d)
            return false;
        idx = i;
        return true;
    }

    void get_basic_term(std::istringstream& is, Poly_d& res)
    {
        char var = 'x', ch = is.peek();
        int idx = -1;

        Poly_d tmp;
        NT coeff;
        unsigned int which_case = 0, power = 1;

        if(isdigit(ch)) {
            coeff = policy.read_coefficient(is);
            which_case = 0; // number
        
        } else if(check_var_names(ch, idx)) {
            which_case = 1;
            var = is.get();
            
        } else if(ch =='(') {
        
            // pos is the index of closing parentheses relative 
            // to the opening ')'
            std::size_t pos = match_parenth(is);
            if(pos == -1u)
                throw internal::Parser_exception(
                    "Parentheses do not match in basic term");
        
            std::size_t beg = (std::size_t)is.tellg() + 1u;
            std::string t = is.str().substr(beg, pos-1);
            is.seekg(pos+1, std::ios_base::cur);
            tmp = get_poly(t);
            which_case = 2;

        } else {
            printf("############ch: %c\n", ch);
            throw internal::Parser_exception("Error in parsing basic term");
        }

        // adjust i to be pointed to the next char
        if(is.peek() == '^') {
            is.ignore(1); // ignore one character
            if(!isdigit(is.peek()))
                throw internal::Parser_exception("Incorrect power for basic term");
            is >> CGAL::iformat(power);
            if(!policy.exponent_check(power))
                throw internal::Parser_exception("Power is too large for basic term");
        }

        switch(which_case) {
        case 0:
            coeff = CGAL::ipower(coeff, static_cast< long >(power));
            tmp = Poly_d(coeff);
            break;
        case 1:
            tmp = construct_monomial(idx, power);
            break;
        case 2: // control degree overflow
            int degree = CGAL::total_degree(tmp);
            if(!policy.exponent_check(degree * power))
                throw internal::Parser_exception(
                    "Power is too large for polynomial in basic term ");
            tmp = CGAL::ipower(tmp, static_cast< long >(power));
        }
        res = tmp;
    }
    
    void get_term(std::istringstream& is, Poly_d& res)
    {
        if(is.eof()) {
            res = Poly_d(NT(0));
            return;
        }
        Poly_d mul;
        get_basic_term(is, mul);
        
        char ch = is.peek();
        //printf("getterm next char to be read: %c\n", ch);
        //while(ind < len && cstr[ind] != '+' && cstr[ind] != '-') {
        while(!is.eof() && ch != '+' && ch != '-') {
            //Recursively get the basic terms till we reach the end or see
            // a '+' or '-' sign.
            if(ch == '*') 
                is.ignore(1);
            Poly_d tmp;
            get_basic_term(is, tmp);
            mul *= tmp;
            ch = is.peek();
            //printf("getterm next char to be read: %c\n", ch);
        }
        res = mul;
        //std::cout << "getterm result: " << res << "\n";
    }
    
    Poly_d get_poly(std::string &s)
    {
       //std::cout << "getpoly: " << s << "\n";
        std::size_t len = s.length();
        if(len == 0) // zero polynomial
            return Poly_d(NT(0));
               
        len = s.length();
        std::istringstream is(s);
        Poly_d res;
        // res will be the polynomial in which we accumulate the
        // sum and difference of the different terms.
        if(is.peek() == '-') {
            is.ignore(1);
            get_term(is, res);
            res *= Poly_d(NT(-1)); // negate polynomial
        } else 
            get_term(is, res);
        
        //  here ind either points to +/- or to the end of string        
        while(!is.eof()) {
            
            Poly_d tmp;
            char ch = is.get(); // must be +/-
            get_term(is, tmp);
            
            if(ch == '+')
                res += tmp;
            else if(ch == '-')
                res -= tmp;
            else
                throw internal::Parser_exception(
                    "Illegal character while parsing polynomial");
        }
        return res;
    }
    
    //!@}
protected:
    //! parser policy object
    Policy policy;
      
};  // Polynomial_parser_d

} // namespace CGAL

#endif // CGAL_POLYNOMIAL_PARSER_D_H
