#ifndef CGAL_LEVEL_OF_DETAIL_TREE_H
#define CGAL_LEVEL_OF_DETAIL_TREE_H

namespace CGAL {
namespace Level_of_detail {

template<typename GeomTraits, typename InputRange, typename PointMap>
class Tree
{
private:

  typedef typename InputRange::const_iterator const_iterator;
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;

  std::vector<const_iterator> m_inliers;

  Point_3 m_center;
  FT m_radius;

public:

  Tree() { }

  std::vector<const_iterator>& inliers() { return m_inliers; }
  const std::vector<const_iterator>& inliers() const { return m_inliers; }

  void set_center_and_radius (const Point_3& center, FT radius)
  {
    m_center = center;
    m_radius = radius;
  }

  const Point_3& center() const { return m_center; }
  FT radius() const { return m_radius; }

};


}

}


#endif // CGAL_LEVEL_OF_DETAIL_TREE_H
