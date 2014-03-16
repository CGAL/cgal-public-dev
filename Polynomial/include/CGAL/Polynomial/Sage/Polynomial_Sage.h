#include<iostream>    //cout
#include<stdio.h> //printf
#include<string.h>    //strlen
#include<string>  //string
#include<sys/socket.h>    //socket
#include<arpa/inet.h> //inet_addr
#include<netdb.h> //hostent
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/Sage/Sage_Connection.h>

#ifndef POLYNOMIAL_SAGE_H
#define POLYNOMIAL_SAGE_H

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
  
  return dataoutput;
}


std::string degree(const CGAL::Polynomial_Sage_type_generator<int,2>::Type &F) 
{
  CGAL::set_pretty_mode(std::cout);
  std::ostringstream oStringForSage;
  CGAL::set_pretty_mode(oStringForSage);
  
  //need to implement respect to others 
  oStringForSage << "R.<y> = PolynomialRing(ZZ)\nR.<x> = PolynomialRing(ZZ)\np=(" << F << ").degree()";

  SageConnection p;

  std::string stringForSage = oStringForSage.str();
  std::string dataFromSage = getDataFromSage( stringForSage );
  return dataFromSage;
}

#endif
