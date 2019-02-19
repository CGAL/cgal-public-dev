#ifndef CGAL_LEVELS_OF_DETAIL_COPLANAR_FACES_CONNECTIVITY_H
#define CGAL_LEVELS_OF_DETAIL_COPLANAR_FACES_CONNECTIVITY_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Coplanar_faces_connectivity {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;

    Coplanar_faces_connectivity(
      const std::vector< std::vector<Point_3> >& walls) :
    m_walls(walls),
    m_tolerance(FT(1) / FT(100000)) { 

			m_neighbors.clear();
      m_neighbors.reserve(m_walls.size());

			std::vector<std::size_t> neighbors;
			for (std::size_t i = 0; i < m_walls.size(); ++i) {
					
				find_neighbors(m_walls[i], i, neighbors);
				m_neighbors.push_back(neighbors);
			}
    }

    void get_neighbors(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      CGAL_precondition(query_index >= 0);
      CGAL_precondition(query_index < m_neighbors.size());

      neighbors = m_neighbors[query_index];
    }

  private:
    const std::vector< std::vector<Point_3> >& m_walls;
    const FT m_tolerance;

    std::vector< std::vector<std::size_t> > m_neighbors;

    void find_neighbors(
      const std::vector<Point_3>& wall,
      const std::size_t curr_index,
      std::vector<std::size_t>& neighbors) {

      neighbors.clear();
			for (std::size_t i = 0; i < m_walls.size(); ++i)
				if (i != curr_index && share_edge(wall, m_walls[i]))
					neighbors.push_back(i);
    }

		bool share_edge(
      const std::vector<Point_3>& w1, 
      const std::vector<Point_3>& w2) const {

			for (std::size_t i = 0; i < w1.size(); ++i) {
				const std::size_t ip = (i + 1) % w1.size();

				for (std::size_t j = 0; j < w2.size(); ++j) {
					const std::size_t jp = (j + 1) % w2.size();

					if (are_equal_edges(w1[i], w1[ip], w2[j], w2[jp])) 
						return true;
				}
			}
			return false;
		}

		bool are_equal_edges(
      const Point_3& p1, const Point_3& p2, 
      const Point_3& q1, const Point_3& q2) const {
      
      return 
      (are_equal_points(p1, q1) && are_equal_points(p2, q2)) || 
      (are_equal_points(p1, q2) && are_equal_points(p2, q1));
    }

    bool are_equal_points(const Point_3& p, const Point_3& q) const {

      const FT eps = m_tolerance;
      return 
      (CGAL::abs(p.x() - q.x()) < eps) && 
      (CGAL::abs(p.y() - q.y()) < eps) && 
      (CGAL::abs(p.z() - q.z()) < eps);
    }

  }; // Coplanar_faces_connectivity

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_COPLANAR_FACES_CONNECTIVITY_H
