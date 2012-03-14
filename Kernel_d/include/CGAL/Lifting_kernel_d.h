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

#ifndef CGAL_LIFTING_KERNEL_D
#define CGAL_LIFTING_KERNEL_D

#include "Kernel_d/Indexed_kernel_d.h"
#include "Kernel_d/Lifted_point_d.h"
#include "Kernel_d/Hashed_orientation_d.h"
#include "Kernel_d/Hashed_volume_d.h"
#include "Kernel_d/sort_swap.h"
#include <CGAL/determinant.h>

// The boost implementation of hash tables appeared in version 1.36. If the
// installed version is older, we abort compilation.
#include <boost/version.hpp>
#if BOOST_VERSION < 103600
#error This needs Boost 1.36 or newer
#endif

#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <vector>

#ifndef CGAL_HASH_TABLE_SIZE_LIMIT
// This value comes from empirical observations. When working with rational
// points, coming from doubles, this limits the process size to around 2Gb.
#define CGAL_HASH_TABLE_SIZE_LIMIT 7999999
#endif

namespace CGAL{

// The determinants hash table is a global (or thread-local) variable.
// Since at this point it is not possible to detect the number type used
// for determinants, it will be created as a void*, and later casted to the
// correct type.
#ifdef CGAL_HAS_THREADS
#include <boost/thread/tss.hpp>
static boost::thread_specific_ptr<void*> _det_table;
#else
static void* _det_table=NULL;
#endif

template <class _IK>
class Lifting_kernel_d:public Indexed_point_kernel_d<_IK>{
        template <class _SomeKernel> friend class HashedOrientation;
        template <class _SomeKernel> friend class HashedVolume;
        private:
        typedef Indexed_point_kernel_d<_IK>                     Base_kernel;
        public:
        typedef Lifting_kernel_d<Base_kernel>                   Self;
        typedef typename Base_kernel::Point_d                   Base_point;
        typedef Lifted_point<Base_point,Self>                   Point_d;
        typedef typename Base_kernel::Vector_d                  Vector_d;
        typedef std::vector<size_t>                             Index;
        typedef typename Base_kernel::FT                        FT;
        typedef typename Base_kernel::LA                        LA;
        typedef boost::unordered_map<Index,FT>                  Table;
        typedef typename Table::iterator                        It;

        // To override some functors from the base kernel:
        public:

        typedef HashedOrientation<Self>                         Orientation_d;
        Orientation_d orientation_d_object()const{return Orientation_d();}

        typedef HashedVolume<Self>                              Volume_d;
        Volume_d volume_d_object()const{return Volume_d();}

        struct Point_to_vector_d{
                Vector_d operator()(const Point_d& p)const{
                        return p-CGAL::ORIGIN;
                }
        };
        Point_to_vector_d point_to_vector_d_object()const{
                return Point_to_vector_d();
        }
        struct Vector_to_point_d{
                Point_d operator()(const Vector_d &v)const{
                        return CGAL::ORIGIN+v;
                }
        };
        Vector_to_point_d vector_to_point_d_object()const{
                return Vector_to_point_d();
        }

        public:
        template <class T>
        static void set_lifting(Point_d &point,const T &l){
                point.set_lifting(l);
        }

        private:
        template <class ForwardIterator>
        static void compute_determinant(ForwardIterator first,
                                        ForwardIterator last,
                                        const Index &idx,
                                        FT &ret){
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

        protected:
        template <class ForwardIterator>
        static bool determinant(ForwardIterator first,
                                ForwardIterator last,
                                FT &det){
                size_t d=std::distance(first,last);
                CGAL_assertion_msg(
                        first->dimension()+1==d||first->dimension()==d,
                        "Hashed_determinant_d: needs d or d+1 points");
#ifndef CGAL_HASH_TABLE_DONT_CLEAR
                // Clear the table when it consumes much memory.
                if(get_table().size()>CGAL_HASH_TABLE_SIZE_LIMIT)
                        get_table().clear();
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
                        std::pair<It,bool> ib=get_table().insert(
                                std::make_pair(all_p_ind,FT()));
                        if(ib.second) // The determinant is not hashed.
                                compute_determinant(
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
                                                get_table().insert(
                                                std::make_pair(p_ind,FT()));
                                        if(ib.second) // Minor is not hashed.
                                                compute_determinant(
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

        protected:
        static Table& get_table(){
#ifdef CGAL_HAS_THREADS
                if(_det_table.get()==NULL)
                        _det_table.reset((void**)(new Table()));
                return (Table&)(*_det_table);
#else
                if(_det_table==NULL)
                        _det_table=(void*)(new Table());
                return (Table&)*(Table*)_det_table;
#endif
        };
};

} // namespace CGAL

#endif // CGAL_LIFTING_KERNEL_D
