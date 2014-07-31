// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
// $URL$
// $Id$
//
//
// Author(s)     : Salahuddin Pasha <s9mdpash@stud.uni-saarland.de> 
//                 
// ============================================================================


#include<iostream>    //cout
#include<stdio.h> //printf
#include<string.h>    //strlen
#include<string>  //string
#include<sys/socket.h>    //socket
#include<arpa/inet.h> //inet_addr
#include<netdb.h> //hostent
#include <CGAL/Polynomial/Sage/Sage_Connection.h>
#include <CGAL/Polynomial/Polynomial_type.h>
#include <cstdio>
#include <sstream>

//#include <CGAL/Polynomial/Polynomial.h>
//#include <CGAL/Polynomial_type_generator.h>
//#include <CGAL/Polynomial_traits_d.h>

#ifndef POLYNOMIAL_SAGE_TYPE_H
#define POLYNOMIAL_SAGE_TYPE_H

#define HOST "127.0.1.1";
#define PORT 12345


namespace CGAL {

  template < class NT_, class Rep_ > class Polynomial;
  template < class NT_> class Polynomial_Sage;

  namespace internal {

    template <class NT_, class Rep_> class Polynomial_rep;
    
    template <class NT_> class Polynomial_sage_rep;
    
    template < class NT_> 
      class Polynomial_sage_rep: Polynomial_rep<NT_, internal::Polynomial_sage_rep<NT_> > {
    public:
      typedef NT_ NT;
      typedef std::vector<NT> Vector;
      typedef typename Vector::size_type      size_type;
      typedef typename Vector::iterator       iterator;
      typedef typename Vector::const_iterator const_iterator;
      Vector coeff;
      
    Polynomial_sage_rep() : coeff(), dirty(false) {   }
    Polynomial_sage_rep(Creation_tag, size_type s) : coeff(s,NT(0)), dirty(false) { }
    Polynomial_sage_rep(size_type n, ...) {}
      
      template <class Forward_iterator>
	Polynomial_sage_rep(Forward_iterator first, Forward_iterator last) 
	: coeff(first,last), dirty(false)
	{}
      
      void reduce() {
	while ( coeff.size()>1 && CGAL::is_zero(coeff.back())) coeff.pop_back();
      }
      
      void simplify_coefficients() {
	typename Algebraic_structure_traits<NT>::Simplify simplify;
	for (iterator i = coeff.begin(); i != coeff.end(); i++) {
	  simplify(*i);
	}
      }
      
      friend class Polynomial_Sage<NT>;
      
    private:
      size_t sage_id;
      bool dirty;
      std::string sage_internal;
      //Polynomial_rep<NT_> polynomial_rep;
      
    public:    
      //should return void 
        std::string push_to_sage() 
	{
	   dirty = true;

	  //convert Polynomial_rep to Sage_rep
	}
      
      std::string pull_from_sage(std::string &param)
	{
	  SageConnection connectionToSage;
	  std::string host;
	  
	  host = HOST
	    
	    //connect to host
	    connectionToSage.conn(host , PORT);
	  //send some data
	  connectionToSage.send_data(param);
	  
	  std::ostringstream oStringReceive;
	  oStringReceive << connectionToSage.receive(1024);
	  std::string dataoutput = oStringReceive.str();
	  
	  sage_id = dataoutput;
	  
	  return dataoutput;
	}
      
    };
    
  }


  template < class NT_ >
    class Polynomial_Sage : public Polynomial< NT_, internal::Polynomial_sage_rep<NT_> > {

    std::string sage_internal;
    //size_t sage_address;
    std::string sage_address;
    bool sage_update_status;
    bool cgal_update_status;

    typedef typename internal::Innermost_coefficient_type<NT_>::Type Innermost_coefficient_type; 
    
  public: 
    typedef NT_ NT; 
    typedef internal::Polynomial_sage_rep<NT> Rep;
    typedef Polynomial< NT, internal::Polynomial_sage_rep<NT> > Base;
    typedef typename Rep::Vector    Vector;
    typedef typename Rep::size_type size_type;
    typedef typename Rep::iterator  iterator;
    typedef typename Rep::const_iterator const_iterator;
    typedef Polynomial_Sage<NT> Self; 
    
  protected:
        Vector& coeffs() { return this->ptr()->coeff; }
    const Vector& coeffs() const { return this->ptr()->coeff; }
    
  Polynomial_Sage(internal::Creation_tag f, size_type s)
    : Base(internal::Polynomial_sage_rep<NT>(f,s) )
      {}

    NT& coeff(unsigned int i) {
      CGAL_precondition(!this->is_shared() && i<(this->ptr()->coeff.size()));
      return this->ptr()->coeff[i]; 
    }
    
    void reduce() { this->ptr()->reduce(); }
    
    void reduce_warn() {
      CGAL_precondition( this->ptr()->coeff.size() > 0 );
      if (this->ptr()->coeff.back() == NT(0)) {
        CGAL_warning_msg(false, "unexpected degree loss (zero divisor?)");
        this->ptr()->reduce();
      }
    }
    
  private:
    static Self& get_default_instance(){
#ifdef CGAL_HAS_THREADS  
      static boost::thread_specific_ptr< Self > safe_x_ptr;
      if (safe_x_ptr.get() == NULL) 
	safe_x_ptr.reset(new Self(0));
      return *safe_x_ptr.get();  
#else
      static Self x = Self(0);
      return x;
#endif        
    }

  public:
     Polynomial_Sage() : Base(static_cast<const Base&>(get_default_instance())) {}

     template <class T>
      explicit Polynomial_Sage(const T& a0)
      : Base(Rep(internal::Creation_tag(), 1))
      { coeff(0) = NT(a0); reduce(); simplify_coefficients(); } 
    

  //! construct the constant polynomial a0
    explicit Polynomial_Sage(const NT& a0)
      : Base(Rep (1, &a0))
      { }
      
    //! construct the polynomial a0 + a1*x
    Polynomial_Sage(const NT& a0, const NT& a1)
      : Base(Rep(2, &a0,&a1))
      { reduce(); simplify_coefficients(); }

    //! construct the polynomial a0 + a1*x + a2*x^2
    Polynomial_Sage(const NT& a0,const NT& a1,const NT& a2)
      : Base(Rep(3, &a0,&a1,&a2))
      { reduce(); simplify_coefficients(); }
      
    //! construct the polynomial a0 + a1*x + ... + a3*x^3
    Polynomial_Sage(const NT& a0,const NT& a1,const NT& a2, const NT& a3)
      : Base(Rep(4, &a0,&a1,&a2,&a3))
      { reduce(); simplify_coefficients(); }

    //! construct the polynomial a0 + a1*x + ... + a4*x^4
    Polynomial_Sage(const NT& a0,const NT& a1,const NT& a2, const NT& a3,
        const NT& a4)
      : Base(Rep(5, &a0,&a1,&a2,&a3,&a4))
      { reduce(); simplify_coefficients(); }
      
    //! construct the polynomial a0 + a1*x + ... + a5*x^5
    Polynomial_Sage(const NT& a0,const NT& a1,const NT& a2, const NT& a3,
        const NT& a4, const NT& a5)
      : Base(Rep(6, &a0,&a1,&a2,&a3,&a4,&a5))
      { reduce(); simplify_coefficients(); }

    //! construct the polynomial a0 + a1*x + ... + a6*x^6
    Polynomial_Sage(const NT& a0,const NT& a1,const NT& a2, const NT& a3,
        const NT& a4, const NT& a5, const NT& a6)
      : Base(Rep(7, &a0,&a1,&a2,&a3,&a4,&a5,&a6))
      { reduce(); simplify_coefficients(); }

    //! construct the polynomial a0 + a1*x + ... + a7*x^7
    Polynomial_Sage(const NT& a0,const NT& a1,const NT& a2, const NT& a3,
        const NT& a4, const NT& a5, const NT& a6, const NT& a7)
      : Base(Rep(8, &a0,&a1,&a2,&a3,&a4,&a5,&a6,&a7))
      { reduce(); simplify_coefficients(); }
      
    //! construct the polynomial a0 + a1*x + ... + a8*x^8
    Polynomial_Sage(const NT& a0,const NT& a1,const NT& a2, const NT& a3,
        const NT& a4, const NT& a5, const NT& a6, const NT& a7,
        const NT& a8)
      : Base(Rep(9, &a0,&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8))
      { reduce(); simplify_coefficients(); }



   template <class Forward_iterator>
    Polynomial_Sage(Forward_iterator first, Forward_iterator last)
      : Base(Rep(first,last)) 
      { reduce(); simplify_coefficients(); }


private:
    // NTX not decomposable
    template <class NTX, class TAG >
      CGAL::Sign sign_at_(const NTX& x, TAG) const{
      CGAL_precondition(degree()>=0);
      return CGAL::sign(evaluate(x));
    }
    // NTX decomposable
    
    template <class NTX>
      CGAL::Sign sign_at_(const NTX& x, CGAL::Tag_true) const{
      CGAL_precondition(degree()>=0);
      typedef Fraction_traits<NTX> FT;
      typedef typename FT::Numerator_type Numerator_type;
      typedef typename FT::Denominator_type Denominator_type;
      Numerator_type num;
      Denominator_type den;
      typename FT::Decompose decompose;
      decompose(x,num,den);
      CGAL_precondition(CGAL::sign(den) == CGAL::POSITIVE);

      typedef Coercion_traits< Numerator_type , Denominator_type > CT;
      typename CT::Cast cast;
      return CGAL::sign(evaluate_homogeneous(cast(num),cast(den)));
    }

  public:
    void simplify_coefficients() { this->ptr()->simplify_coefficients(); }
    int degree() const {
      std::cout << push_to_sage() << std::endl; 
      std::string sage_degree = push_to_sage();
      return atoi( sage_degree.c_str() );
      //return static_cast<int>(this->ptr()->coeff.size())-1;
     }

    void convert_to_sage_format()
    {
      //CGAL::set_pretty_mode(std::cout);
      std::ostringstream tmp_os_stream;
      CGAL::set_pretty_mode(tmp_os_stream);
      tmp_os_stream << this;
      sage_internal = tmp_os_stream.str();
      //std::cout << sage_internal << std::endl;
      sage_update_status = true;
    }
    
    std::string push_to_sage() 
    {
	convert_to_sage_format();
	
	std::ostringstream oStringForSage;
	//need to implement respect to others 
	//oStringForSage << "R.<y> = PolynomialRing(ZZ)\nR.<x> = PolynomialRing(ZZ)\nb=(x^2+x+1)\np=hex(id(b))";
	oStringForSage << "R.<y> = PolynomialRing(ZZ)\nR.<x> = PolynomialRing(ZZ)\nb=(" << sage_internal << ")\np=hex(id(b))";
	
	SageConnection p;
	std::string stringForSage = oStringForSage.str();
	std::string dataFromSage = getDataFromSage( stringForSage );
	
	//std::cout << dataFromSage << std::endl;
	
	return dataFromSage; 
    }
    
    std::string getDataFromSage(std::string &param)
    {
	SageConnection connectionToSage;
	std::string host;
	
	host = "127.0.1.1";
	
	//connect to host
	connectionToSage.conn(host , 12345);
	
	//send some data
	connectionToSage.send_data(param);
	
	std::ostringstream oStringReceive;
	oStringReceive << connectionToSage.receive(1024);
	std::string dataoutput = oStringReceive.str();
	
	sage_address = dataoutput;
	
	return dataoutput;
    }


  };
  
  
}

#endif

// branch build
// derivation
// modify polynomial test
// compile -Wall
// redefine member function
// degree, differentiate


//branch build
// modify template arg
//
