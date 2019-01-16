#ifndef CGAL_LEVELS_OF_DETAIL_GROUND_H
#define CGAL_LEVELS_OF_DETAIL_GROUND_H

// STL includes.
#include <vector>

namespace CGAL {

namespace Levels_of_detail {

namespace internal {

  template<class DataStructure>
  class Ground {

  public:
    using Data_structure = DataStructure;

    using Traits = typename DataStructure::Traits;
    using Point_3 = typename Traits::Point_3;

    Ground(Data_structure &data_structure) :
    m_data(data_structure)
    { }

    void make_planar() {

      m_data.planar_ground.clear();

      std::vector<Point_3> bbox(4);
      bbox[0] = Point_3(0,0,0);
      bbox[1] = Point_3(1,0,0);
      bbox[2] = Point_3(1,1,0);
      bbox[3] = Point_3(0,1,0);

      m_data.planar_ground = bbox;
    }

    template<typename VerticesOutputIterator>
    void return_as_polygon(VerticesOutputIterator vertices) const {

      for (size_t i = 0; i < m_data.planar_ground.size(); ++i)
        *(vertices++) = m_data.planar_ground[i];
    }

  private:
    Data_structure &m_data;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_GROUND_H