#include <CGAL/Cartesian_d.h>
#include <CGAL/Lifting_kernel_d.h>
#include <CGAL/assertions.h>

int main(){
        typedef double                                          NT;
        typedef CGAL::Cartesian_d<NT>                           Base;
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

				std::vector<Point> lpoints;
        lpoints.push_back(p);
        lpoints.push_back(q);
        lpoints.push_back(r);
				lpoints.push_back(s);
        
        NT liftarray[] = {16,2,77,29};
        std::vector<NT> lifting (liftarray, liftarray + sizeof(liftarray) / sizeof(NT) );
        
        // Test the lifted Orientation_d predicate with d+1 points.
        CGAL_assertion(Ori()(lpoints.begin(),lpoints.end(),
                             lifting.begin(),lifting.end())==
                       CGAL::COLLINEAR);

        return 0;
}
