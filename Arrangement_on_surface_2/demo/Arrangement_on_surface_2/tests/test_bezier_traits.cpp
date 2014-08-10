#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;

typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;

typedef Rat_kernel::Point_2                             Rat_point_2;
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
  Bezier_traits_2;

typedef Bezier_traits_2::Point_2 Point_2;
typedef Bezier_traits_2::X_monotone_curve_2 X_monotone_curve_2;
typedef Bezier_traits_2::Construct_x_monotone_curve_2 Construct_x_monotone_curve_2;
typedef Bezier_traits_2::Approximate_2 Approximate_2;
typedef Bezier_traits_2::Approximate_number_type Approximate_number_type;

int main( int argc, char *argv[] )
{
  Bezier_traits_2 traits;
  Rat_point_2 p1(1, 2);
  Rat_point_2 p2(1, 3);
  Construct_x_monotone_curve_2 construct_curve =
    traits.construct_x_monotone_curve_2_object( );
  X_monotone_curve_2 xcv = construct_curve( p1, p2 );
  Point_2 bp1 = xcv.source( );
  Point_2 bp2 = xcv.target( );
  Approximate_2 approx = traits.approximate_2_object( );
  Approximate_number_type x_approx = approx( bp1, 0 );
  Approximate_number_type y_approx = approx( bp1, 1 );
  return 0;
}
