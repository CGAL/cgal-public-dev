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
// $URL: svn+ssh://penarand@scm.gforge.inria.fr/svn/cgal/branches/features/Kernel_d-new_models-penarand_vfisikop/Kernel_d/include/CGAL/Indexed_kernel_d.h $
// $Id: Indexed_kernel_d.h 67920 2012-03-01 16:00:51Z penarand $
//
//
// Authors:     Vissarion Fisikopoulos <vissarion@di.uoa.gr>
//              Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_KERNEL_D_INDEXED_KERNEL_D_H
#define CGAL_KERNEL_D_INDEXED_KERNEL_D_H

#include "Indexed_point_d.h"

namespace CGAL{

#ifdef CGAL_HAS_THREADS
#include <boost/thread/tss.hpp>
static boost::thread_specific_ptr<size_t> _max_index;
#else
static size_t _max_index=0;
#endif

template <class _K>
class Indexed_point_kernel_d:public _K{
        private:
        typedef _K                                              Base;
        typedef Indexed_point_kernel_d<Base>                    Self;
        typedef typename Base::Point_d                          Base_point;
        public:
        typedef Indexed_point<Base_point,Self>                  Point_d;

        public:
        static size_t get_index(){
#ifdef CGAL_HAS_THREADS
                if(_max_index.get()==NULL)
                        _max_index.reset(new size_t(0));
                return *_max_index;
#else
                return _max_index;
#endif
        }

        public:
        static size_t get_and_increment_index(){
#ifdef CGAL_HAS_THREADS
                if(_max_index.get()==NULL)
                        _max_index.reset(new size_t(0));
                _max_index.reset(new size_t(*_max_index+1));
                return *_max_index-1;
#else
                return _max_index++;
#endif
        }
};

} // namespace CGAL

#endif // CGAL_KERNEL_D_INDEXED_KERNEL_D_H
