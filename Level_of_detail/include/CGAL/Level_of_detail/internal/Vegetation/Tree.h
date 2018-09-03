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
  typedef typename GeomTraits::Point_2 Point_2;
  typedef typename GeomTraits::Vector_2 Vector_2;

  std::vector<const_iterator> m_inliers;

  Point_3 m_center;
  FT m_radius;
  FT m_height;

public:

  Tree() { }

  std::vector<const_iterator>& inliers() { return m_inliers; }
  const std::vector<const_iterator>& inliers() const { return m_inliers; }

  void set_center_radius_and_height (const Point_3& center, FT radius, FT height)
  {
    m_center = center;
    m_radius = radius;
    m_height = height;
  }

  const Point_3& center() const { return m_center; }
  
  Point_2 center_2() const { return Point_2(m_center.x(), m_center.y()); }
  
  FT radius() const { return m_radius; }
  FT height() const { return m_height; }

  template <typename PartitionElement>
  void to_partition_face (PartitionElement& el, FT edge_length)
  {
    Point_2 center = center_2();

    std::size_t nb_pts = std::size_t (m_radius / edge_length);
    if (nb_pts < 3) nb_pts = 3;

    for (std::size_t j = 0; j < nb_pts; ++ j)
    {
      FT angle = 2. * CGAL_PI * (j / FT(nb_pts));
      Point_2 p = center + Vector_2 (m_radius * std::cos(angle), m_radius * std::sin(angle));
      el.push_back(p);
    }

    el.visibility_label() = Visibility_label::VEGETATION;
  }


};


}

}


#endif // CGAL_LEVEL_OF_DETAIL_TREE_H
