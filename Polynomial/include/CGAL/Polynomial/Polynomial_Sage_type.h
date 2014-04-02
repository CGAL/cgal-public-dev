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

//#include <CGAL/Polynomial/Polynomial.h>
//#include <CGAL/Polynomial_type_generator.h>
//#include <CGAL/Polynomial_traits_d.h>

#ifndef POLYNOMIAL_SAGE_TYPE_H
#define POLYNOMIAL_SAGE_TYPE_H

#define HOST "127.0.1.1";
#define PORT 12345


template <class NT_> 
class Polynomial_sage_rep : public Polynomial_rep< NT_ > {
 private:
  size_t sage_id;
  bool dirty;
  std::string sage_internal;

public:

  Polynomial_sage_rep() 
    {
      bool dirty = false;
      Polynomial_rep< NT_ >::Polynomial_rep();
    }

  Polynomial_sage_rep(size_type n, ...) 
    {
      bool dirty = false;
      Polynomial_rep< NT_ >::Polynomial_rep(size_type n, ...);
    }


  ~Polynomial_sage_rep() //destructor
    {

    };

  std::string push_to_sage() 
    {
      //convert Polynomial_rep to Sage_rep
    }
  
  std::string getDataFromSage(std::string &param)
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

      sage_address = dataoutput;
      
      return dataoutput;
    }

};
 

template < class NT_ >
class Sage_polynomial : public Polynomial< NT_, Polynomial_sage_rep< NT_ > > {
 
};


#endif
