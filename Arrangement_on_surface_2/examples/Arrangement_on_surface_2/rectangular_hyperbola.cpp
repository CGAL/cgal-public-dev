#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef Kernel::FT                                          NT;
typedef CGAL::Arr_rectangular_hyperbola_with_iso_asymptotes_traits_2<Kernel>
                                                            Traits;
typedef Traits::X_monotone_curve_2                          X_monotone_curve_2;
typedef Traits::Curve_2                                     Curve_2;
typedef Traits::Point_2                                     Point_2;
typedef CGAL::Arrangement_2<Traits>                         Arrangement;

int main()
{
  Traits traits;

  Traits::Construct_x_monotone_curve_2 ctr_x_curve =
    traits.construct_x_monotone_curve_2_object();
  Traits::Construct_curve_2 ctr_curve = traits.construct_curve_2_object();

  Point_2 left(0,0), right(1,1);

  std::cout << "1" << std::endl;

  X_monotone_curve_2 xc1 =
    ctr_x_curve(false, 0, 0, -1, left, right, true, true, true, true, true);

  std::cout << "2" << std::endl;
  
  Point_2 p1(0,0), p2(0,1);
  Curve_2 c1 = ctr_curve(true, 0, 0, -1, p1, p2, true, true, true, true, true);

  std::cout << "3" << std::endl;
  
  Arrangement arr(&traits);
  CGAL::insert(arr, xc1);

  std::cout << "4" << std::endl;

  CGAL::insert(arr, c1);

  std::cout << "5" << std::endl;

  return 0;
}
