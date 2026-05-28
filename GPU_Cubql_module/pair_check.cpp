#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Triangle_3 = K::Triangle_3;

int main(int argc, char** argv)
{
    // Define Triangle A with full precision
    Triangle_3 tri_A = Triangle_3(
        Point_3(0.26442988429999997, 0.56879958470000003, 0.71239627849999998),
        Point_3(0.24675713530000001, 0.48177735389999998, 0.75283349430000002),
        Point_3(0.35325841520000001, 0.53954155599999998, 0.78478380130000003)
    );

    // Mesh B triangles from all 6 pairs
    Triangle_3 pairs[] = {
        // Pair 1 (Triangle 12466)
        Triangle_3(
            Point_3(0.28172562919999999, 0.57147973100000005, 0.77439370819999997),
            Point_3(0.32367477150000001, 0.53441684369999998, 0.76945816600000005),
            Point_3(0.31350119389999997, 0.55834523759999999, 0.82075131040000004)
        ),
        // Pair 2 (Triangle 12729)
        Triangle_3(
            Point_3(0.32367477150000001, 0.53441684369999998, 0.76945816600000005),
            Point_3(0.32506398180000001, 0.57775415460000001, 0.78167728790000002),
            Point_3(0.31350119389999997, 0.55834523759999999, 0.82075131040000004)
        ),
        // Pair 3 (Triangle 12465)
        Triangle_3(
            Point_3(0.25024791400000002, 0.60007261690000002, 0.7689953383),
            Point_3(0.32367477150000001, 0.53441684369999998, 0.76945816600000005),
            Point_3(0.28172562919999999, 0.57147973100000005, 0.77439370819999997)
        ),
        // Pair 4 (Triangle 12728)
        Triangle_3(
            Point_3(0.32367477150000001, 0.53441684369999998, 0.76945816600000005),
            Point_3(0.29182013979999999, 0.59704641700000005, 0.74041401610000002),
            Point_3(0.32506398180000001, 0.57775415460000001, 0.78167728790000002)
        ),
        // Pair 5 (Triangle 12727)
        Triangle_3(
            Point_3(0.31132459089999998, 0.54386790289999998, 0.77399981799999995),
            Point_3(0.29182013979999999, 0.59704641700000005, 0.74041401610000002),
            Point_3(0.32367477150000001, 0.53441684369999998, 0.76945816600000005)
        ),
        // Pair 6 (Triangle 12464)
        Triangle_3(
            Point_3(0.25024791400000002, 0.60007261690000002, 0.7689953383),
            Point_3(0.31132459089999998, 0.54386790289999998, 0.77399981799999995),
            Point_3(0.32367477150000001, 0.53441684369999998, 0.76945816600000005)
        )
    };

    // Check intersections for all pairs
    std::cout << "Checking intersections:\n" << std::endl;
    for (int i = 0; i < 6; ++i) {
        bool intersects = CGAL::do_intersect(tri_A, pairs[i]);
        std::cout << "Pair " << (i + 1) << ": " << (intersects ? "INTERSECT" : "NO INTERSECT") << std::endl;
    }

    return EXIT_SUCCESS;
}