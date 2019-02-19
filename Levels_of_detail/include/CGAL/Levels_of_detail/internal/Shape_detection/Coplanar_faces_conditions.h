#ifndef CGAL_LEVELS_OF_DETAIL_COPLANAR_FACES_CONDITIONS_H
#define CGAL_LEVELS_OF_DETAIL_COPLANAR_FACES_CONDITIONS_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Coplanar_faces_conditions {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;

    typename Traits::Construct_cross_product_vector_3 cross_product_3;
    typename Traits::Compute_scalar_product_3 scalar_product_3;

    Coplanar_faces_conditions(
      const std::vector< std::vector<Point_3> >& walls) : 
    m_walls(walls),
    m_tolerance(FT(1) / FT(100000)) 
    { }

    bool belongs_to_region(
      const std::size_t query_index, 
      const std::vector<std::size_t>& region) const {

      CGAL_precondition(region.size() > 0);
      
      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_walls.size());

      return are_coplanar(query_index, region[0]);
    }

    bool is_valid_region(const std::vector<std::size_t>& region) const {
      return region.size() >= 1;
    }

    void update(const std::vector<std::size_t>&) {
      // skipped!
    }

  private:
    const std::vector< std::vector<Point_3> >& m_walls;
    const FT m_tolerance;

    bool are_coplanar(const std::size_t id1, const std::size_t id2) const {

      const std::vector<Point_3>& w1 = m_walls[id1];
			const std::vector<Point_3>& w2 = m_walls[id2];

			std::size_t count = 0;
			for (std::size_t i = 0; i < w1.size(); ++i) {

				const std::size_t ip = (i + 1) % w1.size();
				const std::size_t ipp = (i + 2) % w1.size();

				for (std::size_t j = 0; j < w2.size(); ++j)
					if (is_coplanar(w1[i], w1[ip], w1[ipp], w2[j])) 
            ++count;
			}
			return count == w1.size() * w2.size();
    }

    bool is_coplanar(
      const Point_3& p1, const Point_3& p2, 
      const Point_3& p3, const Point_3& p4) const {

			const Vector_3 v1 = Vector_3(p1, p2);
			const Vector_3 v2 = Vector_3(p1, p3);
			const Vector_3 v3 = Vector_3(p1, p4);

			const Vector_3 v4 = cross_product_3(v2, v3);
			const FT result = scalar_product_3(v1, v4);

			return CGAL::abs(result) < m_tolerance;
		}

  }; // Coplanar_faces_conditions

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_COPLANAR_FACES_CONDITIONS_H
