#include <CGAL/Lifting_kernel_d/sort_swap.h>
//#include <CGAL/Lifting_kernel_d/determinant.h>

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

template <class _K,
          typename _P=typename _K::Point_d,
          typename _FT=typename _K::FT>
class Lifting_kernel_d:
public _K{
        private:
        typedef _K                                              Kernel_d;
        typedef _FT                                             FT;
        public:
        typedef _P                                              Point_d;
        private:
        typedef std::vector<const Point_d&>                     Index;
        typedef boost::unordered_map<Index,FT>                  Table;
        typedef typename Table::iterator                        It;

        public:
        struct Orientation_d;
};

template <class _K,typename _P,typename _FT>
struct Lifting_kernel_d<_K,_P,_FT>::Orientation_d{

        typedef _K                                              Kernel_d;

        template <class PointInputIterator>
        CGAL::Orientation
        operator()(PointInputIterator first,PointInputIterator last)const{
                return Kernel_d().orientation_d_object()(first,last);
        }

        template <class PointInputIterator,class LiftingInputIterator>
        CGAL::Orientation
        operator()(PointInputIterator,PointInputIterator,
                   LiftingInputIterator,LiftingInputIterator)const{
                // TODO
                return COLLINEAR;
        }
};

} // namespace CGAL
