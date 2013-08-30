// Copyright (c) 2009-2013 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://eric@scm.gforge.inria.fr/svn/cgal/branches/features/Polynomial-Polynomial_parser-asm/Polynomial/include/CGAL/Polynomial/Polynomial_type.h $
// $Id: Polynomial_type.h 67240 2012-01-18 09:57:46Z afabri $
//
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-inf.mpg.de>
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

//! default policy: reads in a polynomial of type \c Polynomial_d_ from the
//! input stream without any conversions
template < class Polynomial_d_ >
struct Default_parser_policy {

    //! first template argument type
    typedef Polynomial_d_ Polynomial_d;
    //! coefficient type
    typedef typename  CGAL::Polynomial_traits_d< Polynomial_d >::
        Innermost_coefficient_type Innermost_coefficient_type;

    enum CoefficientTypeID {
        COEFF_INTEGER,
        COEFF_RATIONAL,
        COEFF_FLOAT
    };

    //! default constructor which optionally initilizes variable names
    Default_parser_policy(const char *var_names_ = 0) {
        _set_var_names(var_names_);
    }
        
    //! reads in the coefficient from the input stream and converts it
    //! to \c Innermost_coefficient_type type using provided type coercion
    //! \c is points to the first digit character of the coefficient
    Innermost_coefficient_type read_coefficient(std::istream& is) const {
        return read_coeff_proxy(is, COEFF_INTEGER);
    }

    //! checking for degree overflow: can be used in real-time applications
    //! checks if the total degree \c deg of a monomial lies within given boundaries
    //! and returns \c true in this case
    virtual bool exponent_check(unsigned deg) const {
        static_cast< void > (deg);
        return true;
    }
    
    //! computes a zero-based index of a polynomial variable name \c ch where
    //! 0 corresponds to an innermost variable while d-1 corresponds to an
    //! outermost variable (d is the number of polynomial variables)
    //! if \c ch is not a valid variable name, returns \c false
    virtual bool check_var_names(char ch, int& idx) const {

        unsigned i;
        ch = tolower(ch);
        for(i = 0; i < var_names.size() && var_names[i] != ch; i++);
        if(i == var_names.size())
            return false;
        idx = (int)i;
        return true;
    }

    virtual Innermost_coefficient_type read_coeff_proxy(std::istream& is, CoefficientTypeID) const {
        return _read_coeff< Innermost_coefficient_type >(is);
    }

    virtual ~Default_parser_policy() {
    }
    
protected:

    template < class NT >
    Innermost_coefficient_type _read_coeff(std::istream& is) const {
        NT cf;
        is >> CGAL::iformat(cf);
        typename CGAL::Coercion_traits< NT, Innermost_coefficient_type >::Cast cvt;
        return cvt(cf);
    }
    
    void _set_var_names(const char *str) {
        
        if(str == 0) { //! set default variable names
            str = "xyzw";
        }    
        int n = strlen(str);
        var_names = std::vector< char >(n);
        for(int i = 0; i < n; i++) {
            var_names[i] = tolower(str[i]);
        }
    }
    
    //! variable names listed in the order from innermost to the outermost variable
    std::vector< char > var_names;
};

//! This parser policy allows to mix integer and rational coefficients in a
//! single equation. Appropriate type coercion with the \c Polynomial_d_ coefficient
//! type must be provided
template < class Polynomial_d_ >
struct Mixed_rational_parser_policy :
        public Default_parser_policy< Polynomial_d_ > {

    //! template argument type
    typedef Polynomial_d_ Polynomial_d;
    //! base class
    typedef Default_parser_policy< Polynomial_d > Base;
    //! type of polynomial coefficient
    typedef typename Base::Innermost_coefficient_type Innermost_coefficient_type;

    typedef typename CGAL::Get_arithmetic_kernel< Innermost_coefficient_type >::
             Arithmetic_kernel AK;
    //! integer number type
    typedef typename AK::Integer Integer;
    //! rational number type
    typedef typename AK::Rational Rational;
    //! input coefficient types
    typedef typename Base::CoefficientTypeID CoefficientTypeID;

     //! default constructor which optionally initializes variable names 
    Mixed_rational_parser_policy(const char *var_names_ = 0) :
        Base(var_names_) {
    }
    
    //! reads in a rational or integer coefficient from the input stream
    //! depending on whether '/' has been met during parsing and converts
    //! it to the \c Innermost_coefficient_type type using provided coercion
    //! \c is points to the first digit character of the coefficient
    //! throws \c Parser_exception if error occurred
    Innermost_coefficient_type read_coefficient(std::istream& is) const {

        std::stringstream buf;

        bool read_rational; // recognised as rational coeff
        if(!_eat_characters(is, buf, read_rational))
            throw internal::Parser_exception("Wrong coefficient type");

        return read_coeff_proxy(buf, read_rational ?
            Base::COEFF_RATIONAL : Base::COEFF_INTEGER);
    }

    virtual Innermost_coefficient_type read_coeff_proxy(std::istream& is, CoefficientTypeID type)
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
//! \c Polynomial_d_ coefficient type must be provided
//! \c FP_rounding optional parameter can be used to round the floating-point
//! number to desired precision before returning it to the parser
template < class Polynomial_d_/*, class FP_rounding = CGAL::Identity<KeyType_>*/ >
struct Mixed_floating_point_parser_policy :
        public Mixed_rational_parser_policy< Polynomial_d_ > {

    //! template argument type
    typedef Polynomial_d_ Polynomial_d;
    //! base class
    typedef Mixed_rational_parser_policy< Polynomial_d > Base;
    //! type of polynomial coefficient
    typedef typename Base::Innermost_coefficient_type Innermost_coefficient_type;

    typedef typename CGAL::Get_arithmetic_kernel< Innermost_coefficient_type >::
             Arithmetic_kernel AK;
    //! integer number type
    typedef typename AK::Integer Integer;
    //! rational number type
    typedef typename AK::Rational Rational;
    //! bigfloat number type
    typedef typename CGAL::Bigfloat_interval_traits<
            typename AK::Bigfloat_interval >::Bound BigFloat;
    //! input coefficient types
    typedef typename Base::CoefficientTypeID CoefficientTypeID;

    //! default constructor which optionally initializes variable names 
    Mixed_floating_point_parser_policy(const char *var_names_ = 0) :
        Base(var_names_) {
    }

    //! allows mixing of integer, rational and floating-point coefficients
    //! in a single equation
    //! \c is points to the first digit character of the coefficient
    //! throws \c Parser_exception if an error occurred
    Innermost_coefficient_type read_coefficient(std::istream& is) const {

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

    virtual Innermost_coefficient_type read_coeff_proxy(std::istream& is, CoefficientTypeID type)
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
template <class Polynomial_d_, class ParserPolicy =
                    Default_parser_policy< Polynomial_d_ > >
struct Polynomial_parser_d
{
    //!\name public typedefs
    //!@{

    //! this instance's template argument
    typedef Polynomial_d_ Polynomial_d;
    //! this instance's second template argument
    typedef ParserPolicy Policy;
    
    //! polynomial policy
    typedef CGAL::Polynomial_traits_d< Polynomial_d > PT;

    //! polynomial coefficient type
    typedef typename Policy::Innermost_coefficient_type NT;
    //! multivariate dimension
    static const int d = CGAL::internal::Dimension< Polynomial_d >::value;

    //!@}
public: 
    /// \name Public methods
    //!@{

    //! default constructor
    Polynomial_parser_d(const Policy& policy_ = Policy()) :
            policy(policy_) {
    }

    //! \brief functor invokation operator
    bool operator()(const std::string& in, Polynomial_d& poly) const {
        
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
                Polynomial_d tmp = get_poly(rhs);
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
    std::size_t match_parenth(std::istream& is) const
    {   
        int count = 0;
        std::size_t pos = 0, start = is.tellg(); // save pos of the first '('
        do {
            if(is.eof() || is.tellg() == (std::streampos)-1) {
                return (std::size_t)-1u; // illegal number of parentheses
            }
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
    Polynomial_d construct_monomial(int idx, unsigned int exp) const {
        static Polynomial_d one(NT(1));
        if(exp == 0) // just return 1
            return one;
        typename PT::Shift shift;
        Polynomial_d monom = shift(one, exp, idx);
        return monom;
    }

    void get_basic_term(std::istringstream& is, Polynomial_d& res) const
    {
        char ch = is.peek();
        int idx = -1;

        Polynomial_d tmp;
        NT coeff;
        unsigned int which_case = 0, power = 1;

        if(isdigit(ch)) {
            coeff = policy.read_coefficient(is);
            which_case = 0; // number
        
        } else if(policy.check_var_names(ch, idx) && idx < d) {
            which_case = 1;
            /* char var = */ is.get();
            
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
            printf("############ ch: %c\n", ch);
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
            tmp = Polynomial_d(coeff);
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
    
    void get_term(std::istringstream& is, Polynomial_d& res) const
    {
        if(is.eof()) {
            res = Polynomial_d(NT(0));
            return;
        }
        Polynomial_d mul;
        get_basic_term(is, mul);
        
        char ch = is.peek();
        //printf("getterm next char to be read: %c\n", ch);
        //while(ind < len && cstr[ind] != '+' && cstr[ind] != '-') {
        while(!is.eof() && ch != '+' && ch != '-') {
            //Recursively get the basic terms till we reach the end or see
            // a '+' or '-' sign.
            if(ch == '*') 
                is.ignore(1);
            Polynomial_d tmp;
            get_basic_term(is, tmp);
            mul *= tmp;
            ch = is.peek();
            //printf("getterm next char to be read: %c\n", ch);
        }
        res = mul;
        //std::cout << "getterm result: " << res << "\n";
    }
    
    Polynomial_d get_poly(std::string &s) const
    {
       //std::cout << "getpoly: " << s << "\n";
        std::size_t len = s.length();
        if(len == 0) // zero polynomial
            return Polynomial_d(NT(0));
               
        len = s.length();
        std::istringstream is(s);
        Polynomial_d res;
        // res will be the polynomial in which we accumulate the
        // sum and difference of the different terms.
        if(is.peek() == '-') {
            is.ignore(1);
            get_term(is, res);
            res *= Polynomial_d(NT(-1)); // negate polynomial
        } else 
            get_term(is, res);
        
        //  here ind either points to +/- or to the end of string        
        while(!is.eof()) {
            
            Polynomial_d tmp;
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
