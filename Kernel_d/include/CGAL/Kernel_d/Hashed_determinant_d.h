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

#ifndef CGAL_KERNEL_D_HASHED_DETERMINANT_D_H
#define CGAL_KERNEL_D_HASHED_DETERMINANT_D_H

#include "sort_swap.h"
#include <CGAL/determinant.h>
#include <CGAL/assertions.h>

#ifndef CGAL_HASH_TABLE_SIZE_LIMIT
// This value comes from empirical observations. When working with rational
// points, coming from doubles, this limits the process size to around 2Gb.
#define CGAL_HASH_TABLE_SIZE_LIMIT 7999999
#endif

namespace CGAL{

template <class _K>
struct HashedDeterminant{
        typedef _K                                              K;
        typedef typename K::Point_d                             Point_d;
        typedef typename K::FT                                  FT;
        typedef typename K::Table                               Table;
        typedef typename Table::iterator                        It;
        typedef typename K::Index                               Index;
        typedef typename K::LA                                  LA;

        // This functor stores in &det the determinant of the matrix formed
        // by the points in the range [first,last). It returns true iff the
        // sign of &det is correct.
        template <class ForwardIterator>
        bool operator()(ForwardIterator first,
                        ForwardIterator last,
                        FT &det)const{
                size_t d=std::distance(first,last);
                CGAL_assertion_msg(
                        first->dimension()+1==d||first->dimension()==d,
                        "Hashed_determinant_d: needs d or d+1 points");
#ifndef CGAL_HASH_TABLE_DONT_CLEAR
                // Clear the table when it consumes much memory.
                if(K::get_table().size()>CGAL_HASH_TABLE_SIZE_LIMIT)
                        K::get_table().clear();
#endif
                det=0;

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

                // The function has two flavors, depending on the input.
                // When given d points, the determinant is looked up in the
                // hash table, computed and inserted if needed, and
                // returned. When given d+1 points, it is necessary to
                // expand along the lifting row, computing and storing
                // minors when needed.
                if(first->dimension()==d){
                        std::pair<It,bool> ib=K::get_table().insert(
                                std::make_pair(all_p_ind,FT()));
                        if(ib.second) // The determinant is not hashed.
                                K::compute_determinant(
                                        first,last,all_ind,ib.first->second);
                        det=ib.first->second;
                }else{
                        // The vector index p_ind contains the indices
                        // (point indices) that correspond to the minor.
                        // The vector index ind contains the indices of
                        // those points in the input vector.
                        Index p_ind,ind;
                        p_ind.reserve(d-1);
                        ind.reserve(d-1);
                        for(size_t i=1;i<d;++i){
                                p_ind.push_back(all_p_ind[i]);
                                ind.push_back(all_ind[i]);
                        }
                        int lift_row=d-2; // Index of points lift coordinate.
                        for(size_t k=0;k<d;++k){
                                // I am not sure if this works correctly. I
                                // don't understand the semantics of the
                                // insert function. It returns true when:
                                // (i) the value was not hashed, or (ii)
                                // there were no value with the same hash
                                // value? I assume (i), but I'm not 100%
                                // sure. I have to check this.
                                if((*(first+all_ind[k]))[lift_row]!=0){
                                        std::pair<It,bool> ib=
                                                K::get_table().insert(
                                                std::make_pair(p_ind,FT()));
                                        if(ib.second) // Minor is not hashed.
                                                K::compute_determinant(
                                                        first,last,ind,
                                                        ib.first->second);
                                        if((lift_row+k)%2)
                                                det-=(*(first+all_ind[k]))
                                                        [lift_row]*
                                                        (ib.first)->second;
                                        else
                                                det+=(*(first+all_ind[k]))
                                                        [lift_row]*
                                                        (ib.first)->second;
                                }
                                if(k+1!=d){
                                        ind[k]=all_ind[k];
                                        p_ind[k]=all_p_ind[k];
                                }
                        }
                }

                return swap;
        }
};

} // namespace CGAL

#endif // CGAL_KERNEL_D_HASHED_DETERMINANT_D_H
