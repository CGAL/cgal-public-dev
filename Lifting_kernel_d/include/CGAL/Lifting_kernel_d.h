#include <CGAL/Lifting_kernel_d/determinant.h>

template <class _K,class _P=_K::Point_d,class _NT=P::NT>
class Lifting_kernel_d:
public K{
        private:
        typedef typename _K                     Kernel_d;
        typedef typename _P                     Point_d;
        typedef typename _NT                    NT;

        private:
        struct Orientation_d;
}

template <class _K,class _P,class _NT>
struct Lifting_kernel_d::Orientation_d{
        template <class PointInputIterator,LiftingInputIterator>
        CGAL::Orientation
        operator()(PointInputIterator,PointInputIterator,
                   LiftingInputIterator,LiftingInputIterator){
                // TODO
        }
}
