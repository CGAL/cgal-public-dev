#include <cassert>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Vector_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef Kernel::FT         FT;
typedef Kernel::Point_2    Point;
typedef Kernel::Vector_2   Vector;

using std::cout; using std::endl;

int main() {
    FT x, y;
    FT a(sqrt(2.432543553246234));
    FT b(exp(3.432452344532652));

    /// Test 1
    x = FT(exp(CGAL::to_double(a)));
    y = FT(sqrt(CGAL::to_double(b)));
    Point p(x, y);
    std::cout<<p.x()<<" "<<p.y()<<std::endl;

    Point q(y,x);
    Vector v(p,q);

    /// Test 2
    FT length = static_cast<FT >(sqrt(CGAL::to_double(v.squared_length())));
    std::cout<<length<<std::endl;

    /// Test 3
    FT dot_product(a*b);
    FT exponent = static_cast<FT >(exp(CGAL::to_double(-dot_product)) );
    std::cout<<exponent<<std::endl;



    return EXIT_SUCCESS;
}
