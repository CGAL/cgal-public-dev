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

#ifndef CGAL_KERNEL_D_HASHED_VOLUME_D_D_H
#define CGAL_KERNEL_D_HASHED_VOLUME_D_D_H

#include "Hashed_determinant_d.h"
#include <CGAL/assertions.h>

namespace CGAL{

template <class _K>
struct HashedVolume{
        typedef _K                                              K;
        typedef typename K::FT                                  FT;
        typedef HashedDeterminant<K>                            HD;

        template <class ForwardIterator>
        FT operator()(ForwardIterator first,ForwardIterator last)const{
                size_t d=std::distance(first,last);
                CGAL_assertion_msg(
                        first->dimension()+1==d||first->dimension()==d,
                        "Hashed_volume_d: needs d or d+1 points");
                FT det;
                bool correct_sign=HD()(first,last,det);
                return correct_sign?det:-det;
        }
};

} // namespace CGAL

#endif // CGAL_KERNEL_D_HASHED_VOLUME_D_D_H
