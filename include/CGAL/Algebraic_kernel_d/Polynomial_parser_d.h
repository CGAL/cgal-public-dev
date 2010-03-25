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
// Library       : NumeriX
// File          : include/NiX/Polynomial_parser_d.h
// SoX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*! \file Polynomial_parser_d.h
 *  \brief Defines ... ?!
 *  
 *  Parser for .. polynomials
 */

#ifndef CGAL_POLYNOMIAL_PARSER_D_H
#define CGAL_POLYNOMIAL_PARSER_D_H

#include <CGAL/config.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial_traits_d.h>

#include <vector>
#include <sstream>

namespace {

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

// coefficient validity check
template <class _Poly_d >
struct _Default_checker {

    // coefficient type
    typedef typename CGAL::Polynomial_traits_d< _Poly_d >::
            Innermost_coefficient_type NT;

    bool operator()(const NT& x) const {
        return true;
    }
};

} // anonymous namespace 

namespace CGAL {

/*! \brief this functor implements parsing of bivariate polynomials
 * 
 * input format (y is outermost variable), e.g.:
 * (y-1)^4 + (-1)*y^3 + (x + y)^2 - (2123234523*x^2 - 2*y*y*x + 3*x*132123)^3
 * (y + x - 3)^3 = x + y - 123 + y*x^2
 */
template <class _Poly_d, class ValidyChecker = _Default_checker< _Poly_d > >
struct Polynomial_parser_d
{
    //!\name public typedefs
    //!@{

    //! this instance's template argument
    typedef _Poly_d Poly_d;
    //! this instance's second template argument
    typedef ValidyChecker Validy_checker;
    
    //! an instance of univariate polynomial
    //typedef typename Poly_3::NT Poly_d;

    //! polynomial traits
    typedef CGAL::Polynomial_traits_d< Poly_d > PT;

    //! coefficient type
    typedef typename PT::Innermost_coefficient_type NT;

protected:

#ifdef CGAL_POLYNOMIAL_PARSE_FLOATING_POINT
    typedef typename CGAL::Get_arithmetic_kernel< NT >::
            Arithmetic_kernel Arithmetic_kernel;
    typedef typename Arithmetic_kernel::Integer Integer;

    //! multi-precision float NT
    typedef typename Arithmetic_kernel::Bigfloat_interval BFI;
    typedef typename CGAL::Bigfloat_interval_traits<BFI>::Bound 
        Bigfloat; 
#endif // CGAL_POLYNOMIAL_PARSE_FLOATING_POINT

    //! multivariate dimension
    static const int d = CGAL::internal::Dimension< Poly_d >::value;

    //! variable names
    static const int n_vars = 4;
    static const char var_lower[n_vars];
    static const char var_upper[n_vars];

    //!@}
public: 
    /// \name Public methods
    //!@{

    //! if \c max_exp_ > 0 it defines the maximal allowed exponent 
    //! for input polynomial
    //! \c max_fp_bits_ - amount of bits to approximate floating point number
    Polynomial_parser_d(unsigned max_fp_bits_ = -1,
            unsigned max_exp_ = -1) {

        max_exp = max_exp_;
        max_fp_bits = max_fp_bits_;
    }

    //! \brief functor invokation operator
    //!

    bool operator()(const std::string& in, Poly_d& poly) {
        
        try {
            // remove white spaces from the string: look for all possible
            // whitespace characters
            std::string s(in);
            const char *whitespace = " \f\n\r\t\v";
            unsigned int cnt = s.find_first_of(whitespace,0);
            while(cnt != std::string::npos){
                s.erase(cnt, 1);
                cnt = s.find_first_of(whitespace, cnt);
            }
                        
            //to handle the case when there is one '=' sign
            //we compute sides separetely and then subtract one from another
            unsigned int loc;
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
        catch(Parser_exception ex) {
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
    unsigned int match_parenth(std::istream& is)
    {   
        int count = 0;
        unsigned int pos = 0, start = is.tellg(); // save pos of the first '('
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
    Poly_d construct_monomial(int idx, unsigned int exp)
    {
        static Poly_d one(NT(1));
        if(exp == 0) // just return 1
            return one;
        typename PT::Shift shift;
        Poly_d monom = shift(one, exp, idx);
        return monom;
    }

    bool check_var_names(char ch, int& idx) {

        int i;
        for(i = 0; i < n_vars; i++)
            if(ch == var_lower[i] || ch == var_upper[i])
                break;

        if(i == n_vars || i >= d)
            return false;
        idx = i;
        return true;
    }

    void get_basic_term(std::istringstream& is, Poly_d& res)
    {
        char var = 'x', ch = is.peek();
        int idx = -1;
        //std::cout << "getbasicterm: " << is.tellg() << std::endl;
    
        Poly_d tmp;
        NT coeff;
        unsigned int which_case = 0, power = 1;

        if(isdigit(ch)) {
#ifdef CGAL_POLYNOMIAL_PARSE_FLOATING_POINT
            int which_type = -1; // 0 - float, 1 - float scientific, 
                                 // 1 - rational
            std::string ttmp;
            do {
                ttmp.push_back(is.get());
                ch = is.peek();
            
                bool parse_exp = (ch == 'e' || ch == 'E'); 
                if(ch == '.' || ch == '/' || parse_exp) {

                    if(parse_exp && (which_type == -1 || which_type == 0))
                        which_type = 1;
                    else if(which_type != -1)
                        throw Parser_exception("Wrong coefficient");

                    if(parse_exp) {
                        ttmp.push_back(is.get());
                        ch = is.peek();
                        if(ch == '-' || ch == '+') {
                            ttmp.push_back(is.get());
                            ch = is.peek();
                        }
                        if(!isdigit(ch))
                            throw Parser_exception("Wrong coefficient");
                    } 
                    if(which_type == -1)
                        which_type = (ch == '/' ? 2 : 0);
                } else if(!isdigit(ch))
                    break;
            } while(1);
        
            if(which_type < 2) { 
                Bigfloat bf(ttmp.c_str());
                //std::cout << "bf = " << NT(bf) << "; str: " << ttmp << "\n";
                if(max_fp_bits != -1u) { 
                    bf = CGAL::round(bf, max_fp_bits);
                    //std::cout << "bf2 = " << NT(bf) << "\n";
                }
                coeff = static_cast<NT>(bf);
            } else
                coeff = static_cast<NT>(ttmp.c_str());
#else        
            is >> CGAL::iformat(coeff);
#endif // CGAL_POLYNOMIAL_PARSE_FLOATING_POINT

            if(!checker(coeff))
                throw Parser_exception("Coefficient validity check failed");
            which_case = 0; // number        
        
        } else if(check_var_names(ch, idx)) {
            which_case = 1;
            var = is.get();
            
        } else if(ch =='(') {
        
            // pos is the index of closing parentheses relative 
            // to the opening ')'
            unsigned int pos = match_parenth(is);
            if(pos == -1u)
                throw Parser_exception(
                    "Parentheses do not match in basic term");
        
            unsigned int beg = (unsigned int)is.tellg() + 1u;
            std::string t = is.str().substr(beg, pos-1);
            is.seekg(pos+1, std::ios_base::cur);
            tmp = get_poly(t);
            which_case = 2;

        } else 
            throw Parser_exception("Error in parsing basic term");
        
        // adjust i to be pointed to the next char
        if(is.peek() == '^') {
            is.ignore(1); // ignore one character
            if(!isdigit(is.peek()))
                throw Parser_exception("Incorrect power for basic term");
            is >> CGAL::iformat(power);
            if(power >= max_exp)
                throw Parser_exception("Power is too large for basic term");
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
            if(degree * power >= max_exp) 
                throw Parser_exception(
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
        unsigned int len = s.length();
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
                throw Parser_exception(
                    "Illegal character while parsing polynomial");
        }
        return res;
    }

protected:

    Validy_checker checker; 
         
    //! a maximal exponent allowed for input polynomial
    unsigned max_exp;
    unsigned max_fp_bits;
};  // Polynomial_parser_d

template <class _Poly_d, class ValidyChecker >
const char Polynomial_parser_d< _Poly_d, ValidyChecker >::var_lower[] =
        {'x', 'y', 'z', 'w'};

template <class _Poly_d, class ValidyChecker >
const char Polynomial_parser_d< _Poly_d, ValidyChecker >::var_upper[] =
        {'X', 'Y', 'Z', 'W'};

} // namespace CGAL

#endif // CGAL_POLYNOMIAL_PARSER_D_H


