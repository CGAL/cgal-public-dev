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

#ifndef CGAL_KERNEL_D_HASHED_ORIENTATION_D_H
#define CGAL_KERNEL_D_HASHED_ORIENTATION_D_H

#include <CGAL/assertions.h>

namespace CGAL{

template <class _K>
struct HashedOrientation{
        typedef _K                                              K;
        typedef typename K::Point_d                             Point_d;
        typedef typename K::FT                                  FT;
        typedef typename K::Index                               Index;
        typedef typename K::LA                                  LA;

        template <class ForwardIterator>
        void determ(ForwardIterator first,
                  ForwardIterator last,
                  const Index &idx,
                  const Index &pidx,
                  FT &ret)const{
                if(K::get_table().count(pidx)>0){
                        ret=K::get_table()[pidx];
                        return;
                }
                // TODO: implement this correctly; now it calls CGAL
                // determinant.
                size_t d=idx.size();
                typename LA::Matrix M(d);
                for(size_t col=0;col<d;++col){
                        for(size_t row=0;row<d-1;++row)
                                M(col,row)=(*(first+idx[col]))[row];
                        M(col,d-1)=1;
                }
                ret=LA::determinant(M);
                K::get_table()[pidx]=ret;
                return;
        };

        template <class ForwardIterator>
        Orientation operator()(ForwardIterator first,ForwardIterator last){
                int d=static_cast<int>(std::distance(first,last));
                CGAL_assertion_msg(first->dimension()+1==d,
                                   "Hashed_orientation_d: needs d+1 points");
                Index indices,p_indices;
                indices.reserve(d-1);
                p_indices.reserve(d-1);
                for(size_t i=1;i<d;++i)
                        indices.push_back(i);
                for(ForwardIterator k=first+1;k!=last;++k)
                        p_indices.push_back(k->index());
                FT det(0),minor;
                size_t i=0;
                int lift_row=d-2; // The index of points lift coordinate.
                for(ForwardIterator j=first;j!=last;++j){
                        determ(first,last,indices,p_indices,minor);
                        if((lift_row+i)%2)
                                det-=j->cartesian(lift_row)*minor;
                        else
                                det+=j->cartesian(lift_row)*minor;
                        indices[i]=i;
                        p_indices[i++]=j->index();
                }
                return(det==0?ZERO:det<0?POSITIVE:NEGATIVE);
        }
};

} // namespace CGAL

#endif // CGAL_KERNEL_D_HASHED_ORIENTATION_D_H
