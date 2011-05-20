// Copyright (c) 2011 INRIA Nancy-Grand Est (France).
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
// Authors: Yacine Bouzidi <yacine.bouzidi@inria.fr>
//          Luis Peñaranda <luis.penaranda@gmx.com>
//          Marc Pouget <marc.pouget@inria.fr>
//          Fabrice Rouillier <fabrice.rouillier@inria.fr>

#ifndef CGAL_RS_ALGEBRAIC_2_H
#define CGAL_RS_ALGEBRAIC_2_H

#include <CGAL/RS/rur_2.h>
#include <CGAL/Gmpfi.h>

namespace CGAL{
namespace RS3{

template <class Polynomial_>
class Algebraic_2{
        private:
        typedef Polynomial_                             Polynomial_2;
        typedef Rur_2<Polynomial_2>                     Rur;
        Rur _rur;
        Gmpfi _x,_y;
        public:

        // constructor from the output of the isolation
        // TODO: remove or rewrite in order to obtain a good RUR
        Algebraic_2(const Gmpfi &x,const Gmpfi &y):_x(x),_y(y){};

        // constructor from RUR and isolating intervals
        Algebraic_2(const Rur &r,const Gmpfi &x,const Gmpfi &y):
                _rur(r),_x(x),_y(y){};

        // constructor from two numbers from which Gmpfi can be constructed
        // TODO: remove or rewrite in order to obtain a good RUR
        template <class T>
        Algebraic_2(const T &xcoord,const T &ycoord):_x(xcoord),_y(ycoord){};

        // get the RUR
        const Rur& get_rur()const{
                return _rur;
        }

        // get the x-interval
        const Gmpfi& get_x()const{
                return _x;
        }

        // get the y-interval
        const Gmpfi& get_y()const{
                return _y;
        }

}; // class Algebraic_2

// write an algebraic number to a stream
// TODO: binary mode
template <class Polynomial_>
inline std::ostream& operator<<(std::ostream &o,
                                const Algebraic_2<Polynomial_> &x){
        if(is_pretty(o)){
                o<<"Algebraic_real_2("<<
                        x.get_x()<<','<<
                        x.get_y()<<','<<
                        x.get_rur()<<']';
        }else{
                o<<'['<<x.get_x()<<','<<x.get_y()<<']';
        }
        return o;
}

} // namespace RS3
} // namespace CGAL

#endif  // CGAL_RS_ALGEBRAIC_2_H
