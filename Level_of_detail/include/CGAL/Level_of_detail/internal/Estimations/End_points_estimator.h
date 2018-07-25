#ifndef CGAL_LEVEL_OF_DETAIL_END_POINTS_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_END_POINTS_ESTIMATOR_H

// LOD includes.
#include <CGAL/Level_of_detail/internal/utils.h>

#include <CGAL/barycenter.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits>
class End_points_estimator {

public:
  using Kernel = GeomTraits;
            
  using FT      = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;

  typename Kernel::Compute_squared_distance_2 squared_distance_2;

  struct Pair_from_point
  {
    typedef Point_2 argument_type;
    typedef std::pair<Point_2, FT> result_type;
    
    result_type operator() (const argument_type& p) const
    {
      return std::make_pair (p, FT(1));
    }
    
  };

  void find_end_points_wrt_barycentre_2(const std::vector<Point_2>& elements, Point_2 &a, Point_2 &b) const {
    CGAL_precondition(elements.size() > 1);

    // Find barycentre.
    Point_2 barycentre
      = CGAL::barycenter 
        (boost::make_transform_iterator
         (elements.begin(), Pair_from_point()),
         boost::make_transform_iterator
         (elements.end(), Pair_from_point()));

    // Find the first end point.
    find_furthest_point_2(elements, barycentre, a);

    // Find the second end point.
    find_furthest_point_2(elements, a, b);
  }

  void find_furthest_point_2(const std::vector<Point_2>& elements, const Point_2 &reference_point, Point_2 &resulting_point) const {
    
    CGAL_precondition(elements.size() > 0);

    FT max_distance = FT(0);
    std::size_t idx = 0;

    for (std::size_t i = 0; i < elements.size(); ++ i) {
      const Point_2 &point = elements[i];

      const FT distance = squared_distance_2(point, reference_point);
      if (distance > max_distance) {
                        
        max_distance = distance;
        idx = i;
      }
    }
    resulting_point = elements[idx];
  }
};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_END_POINTS_ESTIMATOR_H
