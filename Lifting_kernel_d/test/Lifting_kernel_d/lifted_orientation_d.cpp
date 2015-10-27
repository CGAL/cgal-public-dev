#include <CGAL/Cartesian_d.h>
#include <CGAL/Lifting_kernel_d.h>
#include <CGAL/assertions.h>

int main(){
        typedef double                                          NT;
        typedef CGAL::Cartesian_d<NT>                           Base;
        typedef Base::Point_d                                   Point;
        typedef CGAL::Lifting_kernel_d<Base>                    K;
        typedef K::Orientation_d                                Ori;

        Point p(5,-1),q(3,-4),r(1,2),s(4,2),t(3,0);

        K mykernel;
        std::vector<Point> points;
        points.push_back(p);
        points.push_back(q);
        points.push_back(r);

        // Test Orientation_d predicates with d+1 points.
        CGAL_assertion(Ori()(points.begin(),points.end(),
                             points.begin(),points.end())==
                       CGAL::COUNTERCLOCKWISE);

        // Test Orientation_d predicates with d points.
        /*points.pop_back();
        CGAL_assertion(Ori()(points.begin(),points.end())==
                       CGAL::COLLINEAR);
         */

        return 0;
}
