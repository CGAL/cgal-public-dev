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

public:
  struct Model_3
  {
    cpp11::array<FT, 5> height;
    cpp11::array<FT, 3> width;
  };

private:

  std::vector<const_iterator> m_inliers;

  Point_3 m_center;
  FT m_radius;
  FT m_height;

  Model_3 m_model;

public:

  Tree() { }

  std::vector<const_iterator>& inliers() { return m_inliers; }
  const std::vector<const_iterator>& inliers() const { return m_inliers; }

  Model_3& model() { return m_model; }

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
    if (nb_pts < 12) nb_pts = 12;

    for (std::size_t j = 0; j < nb_pts; ++ j)
    {
      FT angle = 2. * CGAL_PI * (j / FT(nb_pts));
      Point_2 p = center + Vector_2 (m_radius * std::cos(angle), m_radius * std::sin(angle));
      el.push_back(p);
    }

    el.visibility_label() = Visibility_label::VEGETATION;
  }

  template <typename OutputIterator>
  void trunk_to_2d_outline (FT edge_length, OutputIterator output)
  {
    Point_2 center = center_2();
    FT radius = m_radius / FT(6);

    std::size_t nb_pts = std::size_t (radius / edge_length);
    if (nb_pts < 12) nb_pts = 12;

    for (std::size_t j = 0; j < nb_pts; ++ j)
    {
      FT angle = 2. * CGAL_PI * (j / FT(nb_pts));
      Point_2 p = center + Vector_2 (radius * std::cos(angle), radius * std::sin(angle));
      *(output ++) = p;
    }
  }

  template <typename FaceHandle, typename PointOutputIterator, typename TriangleOutputIterator>
  void to_3d_model(const std::vector<FaceHandle>& ground_faces,
                   PointOutputIterator points,
                   TriangleOutputIterator faces) const
  {
    std::map<Point_2, FT> ground_points;
    for (std::size_t i = 0; i < ground_faces.size(); ++ i)
      for (std::size_t j = 0; j < 3; ++ j)
      {
        const Point_2& p = ground_faces[i]->vertex(j)->point();
        ground_points.insert (std::make_pair (p, ground_faces[i]->info().height(j)));
      }

    std::size_t nb_pts = ground_points.size();

    Point_2 center = center_2();
    FT radius = m_radius / FT(6);
    
    Point_2 p2_max = center;
    FT h_max = m_model.height[4];
    Point_3 p_max (p2_max.x(), p2_max.y(), h_max);
    *(points ++) = p_max;
      
    for (std::size_t j = 0; j < nb_pts; ++ j)
    {
      FT angle = 2. * CGAL_PI * (j / FT(nb_pts));
      
      Point_2 p2_trunk = center + Vector_2 (radius * std::cos(angle), radius * std::sin(angle));
      FT h_trunk = ground_points[p2_trunk];
      Point_3 p_trunk (p2_trunk.x(), p2_trunk.y(), h_trunk);
      *(points ++) = p_trunk;
      
      Point_2 p2_min = center + Vector_2 (radius * std::cos(angle), radius * std::sin(angle));
      FT h_min = m_model.height[0];
      Point_3 p_min (p2_min.x(), p2_min.y(), h_min);
      *(points ++) = p_min;

      for (std::size_t i = 1; i < 4; ++ i)
      {
        Point_2 p2 = center + Vector_2 (m_model.width[i-1] * std::cos(angle), m_model.width[i-1] * std::sin(angle));
        FT h = m_model.height[i];
        Point_3 p (p2.x(), p2.y(), h);
        *(points ++) = p;
      }
    }

    
    for (std::size_t j = 0; j < nb_pts; ++ j)
    {
      std::size_t l = 1 + 5 * j;
      std::size_t r = 1 + 5 * ((j+1) % nb_pts);

      for (std::size_t i = 0; i < 4; ++ i)
      {
        *(faces ++) =  (make_array(l+i, r+i, l+1+i));
        *(faces ++) = (make_array(l+1+i, r+i, r+1+i));
      }
      *(faces ++) = (make_array(l+4, r+4, 0));
    }
  }


};


}

}


#endif // CGAL_LEVEL_OF_DETAIL_TREE_H
