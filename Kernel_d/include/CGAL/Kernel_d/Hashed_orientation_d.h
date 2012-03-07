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
#include "sort_swap.h"

namespace CGAL{

template <class _K>
class HashedOrientation{
        private:
        typedef _K                                              K;
        typedef typename K::Point_d                             Point_d;
        typedef typename K::FT                                  FT;
        typedef typename K::Table                               Table;
        typedef typename Table::iterator                        It;
        typedef typename K::Index                               Index;
        typedef typename K::LA                                  LA;

        private:
        template <class ForwardIterator>
        void determ(ForwardIterator first,
                  ForwardIterator last,
                  const Index &idx,
                  FT &ret)const{
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
                return;
        };

        public:
        template <class ForwardIterator>
        Orientation operator()(ForwardIterator first,ForwardIterator last){
                int d=static_cast<int>(std::distance(first,last));
                CGAL_assertion_msg(first->dimension()+1==d,
                                   "Hashed_orientation_d: needs d+1 points");
                // The vector all_p_ind contains all the indices of the
                // points whose orientation must be computed. all_ind[i] is
                // defined to i. all_p_ind will be sorted, and all_ind's
                // elements will be swapped in the same way.
                Index all_p_ind,all_ind;
                all_p_ind.reserve(d);
                all_ind.reserve(d);
                for(ForwardIterator k=first;k!=last;++k)
                        all_p_ind.push_back(k->index());
                for(size_t i=0;i<d;++i)
                        all_ind.push_back(i);
                bool swap=inplace_sort_permute_count(all_p_ind,all_ind);
                // The vector index p_ind contains the indices (point
                // indices) that correspond to the minor. The vector index
                // ind contains the indices of those points in the input
                // vector.
                Index p_ind,ind;
                p_ind.reserve(d-1);
                ind.reserve(d-1);
                for(size_t i=1;i<d;++i){
                        p_ind.push_back(all_p_ind[i]);
                        ind.push_back(all_ind[i]);
                }
                FT det(0);//,minor;
                size_t i=0;
                int lift_row=d-2; // The index of points lift coordinate.
                for(ForwardIterator j=first;j!=last;++j){
                        // I am not sure if this works correctly. I don't
                        // understand the semantics of the insert function.
                        // It returns true when:
                        // (i) the value was not hashed, or
                        // (ii) there were no value with the same hash
                        // value?
                        // I assume (i), but I'm not 100% sure. I have to
                        // check this.
                        std::pair<It,bool> ib=K::get_table().insert(
                                std::make_pair(p_ind,FT()));
                        if(ib.second) // The minor is not hashed.
                                determ(first,last,ind,(ib.first)->second);
                        // These commented lines are the original version:
                        //if(K::get_table().count(p_ind)>0){
                        //        minor=K::get_table()[p_ind];
                        //}else{
                        //        determ(first,last,ind,minor);
                        //        K::get_table()[p_ind]=minor;
                        //}
                        if((lift_row+i)%2)
                                det-=j->cartesian(lift_row)*(ib.first)->second;
                        else
                                det+=j->cartesian(lift_row)*(ib.first)->second;
                        ind[i]=all_ind[i];
                        p_ind[i]=all_p_ind[i];
                        ++i;
                }
                std::cout<<"det = "<<det<<std::endl;
                if(det==0)
                        return ZERO;
                if(swap)
                        return (det>0?POSITIVE:NEGATIVE);
                else
                        return (det<0?POSITIVE:NEGATIVE);
        }
};

} // namespace CGAL

#endif // CGAL_KERNEL_D_HASHED_ORIENTATION_D_H
