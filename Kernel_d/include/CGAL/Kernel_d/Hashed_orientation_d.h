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

#include "sort_swap.h"
#include <CGAL/determinant.h>
#include <CGAL/assertions.h>

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
                // This constructs a matrix by subtracting the last column
                // from all the others. Its determinant is the same as the
                // original orientation matrix, but it's one dimension
                // smaller.
                size_t d=idx.size()-1;
                typename LA::Matrix M(d);
                for(size_t col=0;col<d;++col){
                        for(size_t row=0;row<d;++row)
                                M(col,row)=(*(first+idx[col]))[row]-
                                           (*(first+idx[d]))[row];
                }
                switch(d){
                        case 1:
                                ret=M(0,0);
                                break;
                        case 2:
                                ret=CGAL::determinant(M(0,0),M(1,0),
                                                      M(0,1),M(1,1));
                                break;
                        case 3:
                                ret=CGAL::determinant(M(0,0),M(1,0),M(2,0),
                                                      M(0,1),M(1,1),M(2,1),
                                                      M(0,2),M(1,2),M(2,2));
                                break;
                        case 4:
                                ret=CGAL::determinant(
                                        M(0,0),M(1,0),M(2,0),M(3,0),
                                        M(0,1),M(1,1),M(2,1),M(3,1),
                                        M(0,2),M(1,2),M(2,2),M(3,2),
                                        M(0,3),M(1,3),M(2,3),M(3,3));
                                break;
                        case 5:
                                ret=CGAL::determinant(
                                        M(0,0),M(1,0),M(2,0),M(3,0),M(4,0),
                                        M(0,1),M(1,1),M(2,1),M(3,1),M(4,1),
                                        M(0,2),M(1,2),M(2,2),M(3,2),M(4,2),
                                        M(0,3),M(1,3),M(2,3),M(3,3),M(4,3),
                                        M(0,4),M(1,4),M(2,4),M(3,4),M(4,4));
                                break;
                        // In theory, Laplace is faster than Gaussian
                        // elimination for d<6. In practice, we observed
                        // that here it is also faster for d=6.
                        case 6:
                                ret=CGAL::determinant(
                                M(0,0),M(1,0),M(2,0),M(3,0),M(4,0),M(5,0),
                                M(0,1),M(1,1),M(2,1),M(3,1),M(4,1),M(5,1),
                                M(0,2),M(1,2),M(2,2),M(3,2),M(4,2),M(5,2),
                                M(0,3),M(1,3),M(2,3),M(3,3),M(4,3),M(5,3),
                                M(0,4),M(1,4),M(2,4),M(3,4),M(4,4),M(5,4),
                                M(0,5),M(1,5),M(2,5),M(3,5),M(4,5),M(5,5));
                                break;
                        default:
                                // TODO: use something faster than CGAL, Eigen?
                                ret=LA::determinant(M);
                }
                return;
        }

        public:
        template <class ForwardIterator>
        Orientation operator()(ForwardIterator first,
                               ForwardIterator last)const{
                size_t d=std::distance(first,last);
                CGAL_assertion_msg(first->dimension()+1==d,
                                   "Hashed_orientation_d: needs d+1 points");
                // The vector all_p_ind contains all the indices of the
                // points whose orientation must be computed. all_ind[i] is
                // defined to i. all_p_ind will be sorted, and all_ind's
                // elements will be swapped in the same way.
                Index all_p_ind,all_ind;
                all_p_ind.reserve(d);
                all_ind.reserve(d);
                for(size_t i=0;i<d;++i){
                        all_p_ind.push_back((first+i)->index());
                        all_ind.push_back(i);
                }
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
                FT det(0);
                int lift_row=d-2; // The index of points lift coordinate.
                for(size_t k=0;k<d;++k){
                        // I am not sure if this works correctly. I don't
                        // understand the semantics of the insert function.
                        // It returns true when:
                        // (i) the value was not hashed, or
                        // (ii) there were no value with the same hash
                        // value?
                        // I assume (i), but I'm not 100% sure. I have to
                        // check this.
                        if((*(first+all_ind[k]))[lift_row]!=0){
                                std::pair<It,bool> ib=
                                        K::get_table().insert(
                                                std::make_pair(p_ind,FT()));
                                if(ib.second) // The minor is not hashed.
                                        determ(first,last,ind,ib.first->second);
                                if((lift_row+k)%2)
                                        det-=(*(first+all_ind[k]))[lift_row]*
                                                (ib.first)->second;
                                else
                                        det+=(*(first+all_ind[k]))[lift_row]*
                                                (ib.first)->second;
                        }
                        if(k+1!=d){
                                ind[k]=all_ind[k];
                                p_ind[k]=all_p_ind[k];
                        }
                }
                if(det==0)
                        return COPLANAR;
                if(swap)
                        return (det>0?COUNTERCLOCKWISE:CLOCKWISE);
                else
                        return (det<0?COUNTERCLOCKWISE:CLOCKWISE);
        }
};

} // namespace CGAL

#endif // CGAL_KERNEL_D_HASHED_ORIENTATION_D_H
