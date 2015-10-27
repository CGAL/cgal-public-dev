// The lifting kernel is not conceived to work with inexact construction
// kernels. However, we are interested in testing its interface with
// Epick_d, so as to make sure that Lifting_kernel_d will work with
// eventual constructions kernels in the future.

#include <CGAL/Epick_d.h>
#include <CGAL/Lifting_kernel_d.h>
#include <CGAL/assertions.h>

int main(){
        typedef double                                          NT;
        typedef CGAL::Epick_d<CGAL::Dimension_tag<2> >          Base;
        typedef CGAL::Lifting_kernel_d<Base>                    K;
        typedef K::Point_d                                      Point;
        typedef K::Orientation_d                                Ori;

        Point p(5,-1),q(3,-4),r(1,2),s(4,2),t(3,0);

        std::vector<Point> points;
        points.push_back(p);
        points.push_back(q);
        points.push_back(r);

        // Test the original Kernel_d Orientation predicate.
        CGAL_assertion(Ori()(points.begin(),points.end())==CGAL::CLOCKWISE);

        // Test another Kernel_d functor, called from the lifting kernel.
        CGAL_assertion(K::Midpoint_d()(q,t)==Point(3,-2));

        // Test the lifted Orientation_d predicate.
        CGAL_assertion(Ori()(points.begin(),points.end(),
                             points.begin(),points.end())==
                       CGAL::COLLINEAR);

        return 0;
}
