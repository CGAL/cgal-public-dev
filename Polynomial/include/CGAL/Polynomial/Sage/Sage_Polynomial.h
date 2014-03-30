#include<iostream>    //cout
#include<stdio.h> //printf
#include<string.h>    //strlen
#include<string>  //string
#include<sys/socket.h>    //socket
#include<arpa/inet.h> //inet_addr
#include<netdb.h> //hostent
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/Sage/Sage_Connection.h>

#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Polynomial_traits_d.h>

#ifndef SAGE_POLYNOMIAL_H
#define SAGE_POLYNOMIAL_H

template <class T, int d>
  class Sage_Polynomial : public CGAL::Polynomial_type_generator<T,d> {

 private:
  typename CGAL::Polynomial_type_generator<T,d>::Type cgal_internal;
  std::string sage_internal;
  //size_t sage_address;
  std::string sage_address;
  bool sage_update_status;
  bool cgal_update_status;

 public:
  Sage_Polynomial(typename CGAL::Polynomial_type_generator<T,d>::Type cgal_input) 
    {
      cgal_internal = cgal_input;
      cgal_update_status = true;
      sage_update_status = false;
    }


  ~Sage_Polynomial() //destructor
    {
      std::cout << "desctructor called" << std::endl;
    };


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

#endif
