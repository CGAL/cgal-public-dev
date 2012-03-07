// Copyright 2011-2012 National and Kapodistrian University of Athens,
// Greece.
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
// Authors:     Vissarion Fisikopoulos <vissarion@di.uoa.gr>
//              Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_KERNEL_D_LIFTED_POINT_D_H
#define CGAL_KERNEL_D_LIFTED_POINT_D_H

#include <CGAL/Origin.h>

namespace CGAL{

template <class _P,class _K>
class Lifted_point:public _P{
        private:
        typedef _P                                              Base_point;
        typedef _K                                              Base_kernel;
        typedef typename Base_kernel::RT                        RT;
        typedef typename Base_kernel::FT                        FT;
        public:
        typedef Lifted_point<Base_point,Base_kernel>            Self;
        typedef typename Base_point::R                          R;

        // Constructors. They are the same as in Indexed_point_d.

        Lifted_point(int d=0):Base_point(d){}

        Lifted_point(int d,const Origin &o):Base_point(d,o){}

        Lifted_point(int a,int b,int c=1):Base_point(RT(a),RT(b),RT(c)){}

        Lifted_point(const RT& a,const RT& b,const RT& c=1):
                Base_point(a,b,c){}

        Lifted_point(int a,int b,int c,int d):
                Base_point(RT(a),RT(b),RT(c),RT(d)){}

        Lifted_point(const RT &a,const RT &b,const RT &c,const RT &d):
                Base_point(a,b,c,d){}

        template <class InputIterator>
        Lifted_point(int d,InputIterator first,InputIterator last):
                Base_point(d,first,last){}

        template <class InputIterator>
        Lifted_point(int d,InputIterator first,InputIterator last,const RT &D):
                Base_point(d,first,last,D){}

        Lifted_point(const Self &p):Base_point(p){}

        explicit Lifted_point(const Base_point &p):Base_point(p){}

        void set_entry(int i,const FT &x){this->entry(i)=x;}
};

} // namespace CGAL

#endif // CGAL_KERNEL_D_LIFTED_POINT_D_H
