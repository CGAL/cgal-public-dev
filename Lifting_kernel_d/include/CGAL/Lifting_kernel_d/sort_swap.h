#ifndef CGAL_LIFTING_KERNEL_D_SORT_SWAP_H
#define CGAL_LIFTING_KERNEL_D_SORT_SWAP_H

#include <CGAL/assertions.h>
#include <boost/version.hpp>
#if BOOST_VERSION < 104300
#include <boost/detail/algorithm.hpp>
#else
#include <boost/range/algorithm_ext/is_sorted.hpp>
#endif

// This function sorts the index vector v and permutes the input vector
// perm exactly in the same way v was permuted. It returns true iff the
// number of swaps in the permutation is even.
template <class Index,class PermutationVector>
bool inplace_sort_permute_count(Index &v,PermutationVector &perm){
        typedef typename Index::value_type                      elt_t;
        size_t n=v.size();
        CGAL_assertion(v.size()==perm.size());
        bool swaps=true; // The number of swaps used is even.
        for(size_t i=0;i+1<n;++i){
                size_t min=i;
                for(size_t j=i+1;j<n;++j)
                        if(v[j]<v[min])
                                min=j;
                if(min!=i){
                        elt_t tmp=v[min];
                        v[min]=v[i];
                        v[i]=tmp;
                        swaps=!swaps;
                        tmp=perm[min];
                        perm[min]=perm[i];
                        perm[i]=tmp;
                }
        }
        CGAL_assertion(boost::is_sorted(v));
        return swaps;
}

#endif // CGAL_LIFTING_KERNEL_D_SORT_SWAP_H
