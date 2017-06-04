// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s)     : Tianyu Zhou <zhoutianyu16tue@gmail.com>

#include "AlgebraicCurveExpressionParser.h"

AlgebraicCurveExpressionParser::AlgebraicCurveExpressionParser(std::string&poly_expr):
    poly_expr(poly_expr)
{}

void AlgebraicCurveExpressionParser::pre_hanlde_poly_expr( std::string& str )
{
    std::vector<int> neg_sign_index;
    
    for ( int i = 0; i < str.size(); i++ )
    {
        if ( str[i] == '+' )
        {
            str.replace(i, 1, ",");
        }
        
        if ( str[i] == '-')
        {
            neg_sign_index.push_back(i);
        }

        if ( str[i] == '*')
        {
            str.replace(i, 1, " ");
        }
    }
    
    for ( int i = int(neg_sign_index.size()-1); i >=0 ; i-- )
    {
        str.insert(neg_sign_index[i], "," );
    }
}

int AlgebraicCurveExpressionParser::extract_poly_coefficient(std::string& poly_expr, struct term& term)
{
    int ret_val = 0;
    
    if ( poly_expr[0] == '-'
        && (poly_expr[1]=='x' || poly_expr[1]=='y') )
    {
        term.coefficient = -1;
        ret_val = 1;
    }
    else if (poly_expr[0]=='x' || poly_expr[0]=='y')
    {
        term.coefficient = 1;
        ret_val = 0;
    }
    else{
        term.coefficient = std::stoi(poly_expr);
        std::string tmp = std::to_string(term.coefficient);
        ret_val = int(tmp.size());
    }
    return ret_val;
}

void AlgebraicCurveExpressionParser::extract_poly_exponents(std::string& sub_poly_expr, struct term& term)
{
    term.x_exponent = 0;
    term.y_exponent = 0;
    
    std::size_t found_x = sub_poly_expr.find('x');
    if (found_x != std::string::npos)
    {
        if(sub_poly_expr[found_x+1] == '^')
        {
            term.x_exponent = std::stoi(sub_poly_expr.substr(found_x+2));
        }
        else
        {
            term.x_exponent = 1;
        }
    }
    
    std::size_t found_y = sub_poly_expr.find('y');
    if (found_y != std::string::npos)
    {
        if(sub_poly_expr[found_y+1] == '^')
        {
            term.y_exponent = std::stoi(sub_poly_expr.substr(found_y+2));
        }
        else
        {
            term.y_exponent = 1;
        }
    }
}

void AlgebraicCurveExpressionParser::print_poly_term( const term& term )
{
    std::cout << "term\n";
    std::cout << "  coefficient: " << term.coefficient << "\n"
    << "  x_exponent: " << term.x_exponent << "\n"
    << "  y_exponent: " << term.y_exponent << "\n";
}

void AlgebraicCurveExpressionParser::extract_poly_components(std::string& poly_expr, struct term& term)
{
    // Extract the coefficient
    int starting_index = extract_poly_coefficient(poly_expr, term);
        
    // Extract the x y exponents respectively
    std::string sub_poly_expr = poly_expr.substr(starting_index);
    extract_poly_exponents(sub_poly_expr, term);
}

void AlgebraicCurveExpressionParser::extract_poly_terms(std::vector<struct term>& poly_terms)
{
    pre_hanlde_poly_expr(this->poly_expr);
    
    boost::char_separator<char> sep(",");
    boost::tokenizer< boost::char_separator<char> > tokens(this->poly_expr, sep);
    BOOST_FOREACH (const std::string& str, tokens) {
        std::string output = str;
        output.erase(remove_if(output.begin(), output.end(), isspace), output.end());
//            std::cout << output << std::endl<< std::endl;
        
        struct term term;
        extract_poly_components(output, term);
        print_poly_term(term);
        poly_terms.push_back(term);
    }
}