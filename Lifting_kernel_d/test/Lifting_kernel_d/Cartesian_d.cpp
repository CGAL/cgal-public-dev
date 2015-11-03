#include <CGAL/Cartesian_d.h>
#include <CGAL/Lifting_kernel_d.h>
#include <CGAL/assertions.h>

int main(){
        typedef double                                          NT;
        typedef CGAL::Cartesian_d<NT>                           Base;
        typedef CGAL::Lifting_kernel_d<Base>                    K;
        typedef K::Point_d                                      Point;
        typedef K::Orientation_d                                Ori;

        // Standard Kernel

        NT a0[] = {5,-1,16};
        std::vector<NT> b0 (a0, a0 + sizeof(a0) / sizeof(NT) );
        NT a1[] = {3,-4,2};
        std::vector<NT> b1 (a1, a1 + sizeof(a1) / sizeof(NT) );
        NT a2[] = {1,2,77};
        std::vector<NT> b2 (a2, a2 + sizeof(a2) / sizeof(NT) );
        NT a3[] = {4,2,29};
        std::vector<NT> b3 (a3, a3 + sizeof(a3) / sizeof(NT) );
        
        Point p0(3,b0.begin(),b0.end());
        Point p1(3,b1.begin(),b1.end());
        Point p2(3,b2.begin(),b2.end());
        Point p3(3,b3.begin(),b3.end());
        
        std::vector<Point> points;
        points.push_back(p0);
        points.push_back(p1);
        points.push_back(p2);
        points.push_back(p3);

        // Test the original Kernel_d Orientation predicate.
        //std::cout << Ori()(points.begin(),points.end()) << "\n";
        
        // Test another Kernel_d functor, called from the lifting kernel.
        //CGAL_assertion(K::Midpoint_d()(q,t)==Point(3,-2));


        // Lifted Kernel
        // Test the lifted Orientation_d predicate with d+2 points.
       
        Point q0(2,b0.begin(),b0.end()-1);
        Point q1(2,b1.begin(),b1.end()-1);
        Point q2(2,b2.begin(),b2.end()-1);
        Point q3(2,b3.begin(),b3.end()-1);
        
				std::vector<Point> lpoints;
        lpoints.push_back(q0);
        lpoints.push_back(q1);
        lpoints.push_back(q2);
				lpoints.push_back(q3);
        
        NT liftarray[] = {16,2,77,29};
        std::vector<NT> lifting (liftarray, liftarray + sizeof(liftarray) / sizeof(NT) );
        
        //std::cout << Ori()(lpoints.begin(),lpoints.end(),lifting.begin(),lifting.end()) << "\n";
        
        //IMPORTANT NOTE
        //our code gives the opposite sign than the Kernel_d if dim is odd and the same otherwise, this is because we compute the homogenous det while the original kernel computes the  det of the matrix created after subdividing from all columns (points) the first one 
         
        CGAL_assertion(Ori()(points.begin(),points.end())==Ori()(lpoints.begin(),lpoints.end(),lifting.begin(),lifting.end()));

        
        return 0;
}
