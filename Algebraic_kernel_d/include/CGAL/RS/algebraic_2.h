// Copyright (c) 2011 National and Kapodistrian University of Athens (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_ALGEBRAIC_2_H
#define CGAL_RS_ALGEBRAIC_2_H

#include <CGAL/Gmpfi.h>

namespace CGAL{
namespace RS3{

class Algebraic_2_rep{
        public:
        // the intervals containing the x and y solutions of the system
        Gmpfi _x,_y;
        Algebraic_2_rep(){};
        ~Algebraic_2_rep(){};
        private:
        Algebraic_2_rep(const Algebraic_2_rep&);
        Algebraic_2_rep& operator=(const Algebraic_2_rep);
};

class Algebraic_2:Handle_for<Algebraic_2_rep>{
        private:
        typedef Handle_for<Algebraic_2_rep>     Base;
        public:

        // constructor from the output of the isolation
        Algebraic_2(const Gmpfi &x,const Gmpfi &y){
                ptr()->_x=x;
                ptr()->_y=y;
        };

        // constructor from two numbers from which Gmpfi can be constructed
        template <class T>
        Algebraic_2(const T &xcoord,const Y &ycoord){
                Gmpfi x(xcoord),y(ycoord);
                ptr()->_x=x;
                ptr()->_y=y;
        };

        // get the x-interval
        Gmpfi& get_x()const{
                return Ptr()->_x;
        }

        // get the y-interval
        Gmpfi& get_y()const{
                return Ptr()->_y;
        }

}; // class Algebraic_2

} // namespace RS3
} // namespace CGAL

#endif  // CGAL_RS_ALGEBRAIC_2_H
