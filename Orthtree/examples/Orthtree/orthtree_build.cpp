#include <iostream>

#include <CGAL/Epick_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Orthtree.h>
#include <CGAL/Orthtree_traits_point.h>
#include <CGAL/Random.h>

// Type Declarations
using Dimension = CGAL::Dimension_tag<4>;
using Kernel = CGAL::Epick_d<Dimension>;
using Point_d = Kernel::Point_d;
using Point_vector = std::vector<Point_d>;
using Traits = CGAL::Orthtree_traits_point<Kernel, Point_vector>;
using Traits2 = CGAL::Orthtree_traits_point<Kernel, Point_vector, CGAL::Identity_property_map<typename std::iterator_traits<typename Point_vector::iterator>::value_type>, CGAL::Dimension_tag<4>>;
using Orthtree = CGAL::Orthtree<Traits2>;

int main()
{
  CGAL::Random r;

  Point_vector points_dd;
  for (std::size_t i = 0; i < 500; ++ i)
  {
    std::array<double, Dimension::value> init{};
    for (double& v : init)
      v = r.get_double(-1., 1.);
    points_dd.emplace_back (init.begin(), init.end());
  }

  Orthtree orthtree(points_dd);
  orthtree.refine(10, 5);

  std::cout << orthtree.bbox(orthtree.root()).min()[0] << std::endl;
  std::cout << orthtree << std::endl;

  return EXIT_SUCCESS;
}
