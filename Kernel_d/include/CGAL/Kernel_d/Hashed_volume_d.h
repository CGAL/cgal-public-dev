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

#include "sort_swap.h"
#include <CGAL/determinant.h>
#include <CGAL/assertions.h>

#ifndef CGAL_HASH_TABLE_SIZE_LIMIT
#error At this point, CGAL_HASH_TABLE_SIZE_LIMIT must be defined!
#endif

namespace CGAL{

template <class _K>
class HashedVolume{
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
        void volume(ForwardIterator first,
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
        FT operator()(ForwardIterator first,ForwardIterator last)const{
                size_t d=std::distance(first,last);
                CGAL_assertion_msg(first->dimension()==d,
                                   "Hashed_volume_d: needs d points");
#ifndef CGAL_HASH_TABLE_DONT_CLEAR
                // Clear the table when it consumes much memory.
                if(K::get_table().size()>CGAL_HASH_TABLE_SIZE_LIMIT)
                        K::get_table().clear();
#endif
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
                std::pair<It,bool> ib=
                        K::get_table().insert(std::make_pair(all_p_ind,FT()));
                if(ib.second) // The volume is not hashed.
                        volume(first,last,all_ind,ib.first->second);
                return swap?(ib.first->second):-(ib.first->second);
        }
};

} // namespace CGAL

#endif // CGAL_KERNEL_D_HASHED_VOLUME_D_D_H
