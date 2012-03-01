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
// $URL: svn+ssh://penarand@scm.gforge.inria.fr/svn/cgal/branches/features/Kernel_d-new_models-penarand_vfisikop/Kernel_d/include/CGAL/Kernel_d/Indexed_point_d.h $
// $Id: Indexed_point_d.h 67920 2012-03-01 16:00:51Z penarand $
//
//
// Authors:     Vissarion Fisikopoulos <vissarion@di.uoa.gr>
//              Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_KERNEL_D_INDEXED_POINT_D_H
#define CGAL_KERNEL_D_INDEXED_POINT_D_H

#include <CGAL/Origin.h>

namespace CGAL{

template <class _P,class _K>
class Indexed_point:public _P{
        private:
        typedef _P                                              Base_point;
        typedef _K                                              Base_kernel;
        typedef typename Base_kernel::RT                        RT;
        typedef typename Base_kernel::FT                        FT;
        public:
        typedef Indexed_point<Base_point,Base_kernel>           Self;
        typedef typename Base_point::R                          R;

        // Constructors. They are the same as original Point_d, but each
        // constructor creates a new index.

        Indexed_point(int d=0):
                Base_point(d),
                _index(Base_kernel::get_and_increment_index()){}

        Indexed_point(int d,const Origin &o):
                Base_point(d,o),
                _index(Base_kernel::get_and_increment_index()){}

        Indexed_point(int a,int b,int c=1):
                Base_point(RT(a),RT(b),RT(c)),
                _index(Base_kernel::get_and_increment_index()){}

        Indexed_point(const RT& a,const RT& b,const RT& c=1):
                Base_point(a,b,c),
                _index(Base_kernel::get_and_increment_index()){}

        Indexed_point(int a,int b,int c,int d):
                Base_point(RT(a),RT(b),RT(c),RT(d)),
                _index(Base_kernel::get_and_increment_index()){}

        Indexed_point(const RT &a,const RT &b,const RT &c,const RT &d):
                Base_point(a,b,c,d),
                _index(Base_kernel::get_and_increment_index()){}

        template <class InputIterator>
        Indexed_point(int d,InputIterator first,InputIterator last):
                Base_point(d,first,last),
                _index(Base_kernel::get_and_increment_index()){}

        template <class InputIterator>
        Indexed_point(int d,
                      InputIterator first,
                      InputIterator last,
                      const RT &D):
                Base_point(d,first,last,D),
                _index(Base_kernel::get_and_increment_index()){}

        explicit Indexed_point(const Self &p):
                Base_point(p),_index(p.index()){}

        explicit Indexed_point(const Base_point &p):
                Base_point(p),
                _index(Base_kernel::get_and_increment_index()){}

        Self operator+(const Vector_d<R> &v)const
                {return Self(Base_point::operator+(v));}
        Self operator-(const Vector_d<R> &v)const
                {return Self(Base_point::operator-(v));}
        Self& operator+=(const Vector_d<R> &v)
                {return static_cast<Self&>(Self(Base_point::operator+=(v)));}
        Self& operator-=(const Vector_d<R> &v)
                {return static_cast<Self&>(Self(Base_point::operator-=(v)));}

        size_t index()const{return _index;}

        private:
        size_t _index;
};

} // namespace CGAL

#endif // CGAL_KERNEL_D_INDEXED_POINT_D_H

