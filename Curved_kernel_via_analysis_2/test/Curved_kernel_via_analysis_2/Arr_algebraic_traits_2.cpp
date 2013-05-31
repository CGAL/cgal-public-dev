
// Needed temporarily because of local changes
#include <CGAL/Timer.h>
CGAL::Timer overall_timer,symb_timer, refine_timer, solve_timer, test_timer,
    lift_timer,inter_timer;

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

typedef CGAL::CORE_arithmetic_kernel AK;
typedef AK::Integer Integer;
typedef CGAL::Arr_algebraic_segment_traits_2<Integer> Arr_traits_2;
typedef CGAL::Arrangement_2<Arr_traits_2> Arrangement_2;
typedef Arrangement_2::X_monotone_curve_2 Segment_2;
typedef Arrangement_2::Point_2 Point_2;
typedef Arr_traits_2::Curve_2 Curve_2;
typedef Arr_traits_2::Polynomial_2 Polynomial_2;
typedef Arr_traits_2::Algebraic_real_1 Algebraic_real_1;

int main(int argc, char** argv) {

    Arr_traits_2 arr_traits;

    
    Arr_traits_2::Construct_curve_2 construct_curve_2
        = arr_traits.construct_curve_2_object();

    Arr_traits_2::Construct_point_2 construct_point_2
        = arr_traits.construct_point_2_object();

    Arr_traits_2::Construct_segments_2 construct_segments_2
        = arr_traits.construct_segments_2_object();

    CGAL::Polynomial_parser_d<Polynomial_2> parse;

    {
        Arrangement_2 arr(&arr_traits);
        Polynomial_2 f;
        parse("(x^2+y^2-1)*(x^2+(y-1)^2-2)*(2*x-2*y-1)*(x-y-20)",f);
        Curve_2 cv = construct_curve_2(f);
        Point_2 a = construct_point_2(Algebraic_real_1(-5),cv,0);
        Point_2 b = construct_point_2(Algebraic_real_1(5),cv,0);
        std::vector<Segment_2> segments;
        construct_segments_2(cv,a,b,std::back_inserter(segments));
        CGAL::insert(arr,segments.begin(),segments.end());
        parse("(x-5)",f);
        CGAL::insert(arr,construct_curve_2(f));
        parse("(x+5)",f);
        CGAL::insert(arr,construct_curve_2(f));
        std::cout << arr.number_of_vertices() << ","
                  << arr.number_of_edges() << ","
                  << arr.number_of_faces() << std::endl;
    }
    {
        Arrangement_2 arr(&arr_traits);
        Polynomial_2 f;
        parse("x^4+y^4-1",f);
        Curve_2 cv = construct_curve_2(f);
        Point_2 a = construct_point_2(Algebraic_real_1(-1),cv,0);
        Point_2 b = construct_point_2(Algebraic_real_1(1),cv,0);
        Point_2 c = construct_point_2(Algebraic_real_1(0),cv,1);
        std::vector<Segment_2> segments;
        construct_segments_2(cv,c,CGAL::ZERO,std::back_inserter(segments));
        CGAL::insert(arr,segments.begin(),segments.end());
        parse("((y+1)^3-x^2)",f);
        cv = construct_curve_2(f);
        segments.clear();
        construct_segments_2(cv,a,b,std::back_inserter(segments));
        CGAL::insert(arr,segments.begin(),segments.end());
        std::cout << arr.number_of_vertices() << ","
                  << arr.number_of_edges() << ","
                  << arr.number_of_faces() << std::endl;
    }

    {
        Arrangement_2 arr(&arr_traits);
        Polynomial_2 f;
        parse("(x-y)*y*(x-2-y^2)",f);
        Curve_2 cv = construct_curve_2(f);
        Point_2 a = construct_point_2(Algebraic_real_1(1),cv,0);
        std::vector<Segment_2> segments;
        construct_segments_2(cv,a,CGAL::ZERO,std::back_inserter(segments));
        CGAL::insert(arr,segments.begin(),segments.end());
        parse("(x+y)*y*(x+1+y^2)",f);
        cv = construct_curve_2(f);
        segments.clear();
        Point_2 b = construct_point_2(Algebraic_real_1(-1),cv,1);
        Point_2 c = construct_point_2(Algebraic_real_1(0),cv,0);
        construct_segments_2(cv,b,c,std::back_inserter(segments));
        CGAL::insert(arr,segments.begin(),segments.end());
        std::cout << arr.number_of_vertices() << ","
                  << arr.number_of_edges() << ","
                  << arr.number_of_faces() << std::endl;
    }

    {
        Arrangement_2 arr(&arr_traits);
        Polynomial_2 f;
        parse("(y-x^3)*(y^3-x^2)*((x+3)^2+y^2-1)",f);
        Curve_2 cv = construct_curve_2(f);
        Point_2 a = construct_point_2(Algebraic_real_1(-1),cv,0);
        Point_2 b = construct_point_2(Algebraic_real_1(-1),cv,1);
        Point_2 c = construct_point_2(Algebraic_real_1(3),cv,0);
        Point_2 d = construct_point_2(Algebraic_real_1(2),cv,1);
        std::vector<Segment_2> segments;
        construct_segments_2(cv,a,CGAL::ZERO,std::back_inserter(segments));
        construct_segments_2(cv,b,CGAL::ZERO,std::back_inserter(segments));
        construct_segments_2(cv,c,CGAL::NEGATIVE,std::back_inserter(segments));
        construct_segments_2(cv,d,CGAL::POSITIVE,std::back_inserter(segments));

        CGAL::insert(arr,segments.begin(),segments.end());
        std::cout << arr.number_of_vertices() << ","
                  << arr.number_of_edges() << ","
                  << arr.number_of_faces() << std::endl;
    }



/*

    for(int i = 1; i < argc; i++) {
        Polynomial_2 f;
        if(! parse(std::string(argv[i]), f)) {
            std::cerr << "Bad polynomial format" << std::endl;
            std::exit(1);
        } else {
            std::cout << "Adding curve " << argv[i] << std::endl;
        }
        CGAL::insert(arr,construct_curve_2(f));
    }

    std::cout << arr.number_of_vertices() << ","
              << arr.number_of_edges() << ","
              << arr.number_of_faces() << std::endl;
*/
    return 0;
}
