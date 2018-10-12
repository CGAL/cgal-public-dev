#ifndef CGAL_LEVEL_OF_DETAIL_TREE_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_TREE_ESTIMATOR_H

#include <CGAL/Level_of_detail/internal/Vegetation/Tree.h>

#include <CGAL/Classification/Image.h>
#include <CGAL/barycenter.h>

#include <boost/make_shared.hpp>

// TO REMOVE
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>


namespace CGAL {
namespace Level_of_detail {

template<typename GeomTraits, typename InputRange, typename PointMap>
class Tree_estimator
{
public:

  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Circle_2 Circle_2;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename InputRange::const_iterator const_iterator;

  typedef Tree<GeomTraits, InputRange, PointMap> Tree_item;
  
private:

  std::vector<Tree_item>& m_trees;
  PointMap m_point_map;

  struct Point_from_iterator_and_pmap
  {
    typedef const_iterator argument_type;
    typedef std::pair<Point_3, FT> return_type;

    PointMap m_point_map;

    Point_from_iterator_and_pmap (PointMap point_map)
      : m_point_map (point_map) { }

    return_type operator() (const argument_type& arg) const
    {
      return std::make_pair (get (m_point_map, *arg), FT(1));
    }
  };

public:

  Tree_estimator (std::vector<Tree_item>& input, PointMap point_map)
    : m_trees (input), m_point_map (point_map)
  {
  }

  void estimate (const FT& minimum_radius)
  {
    for (std::size_t i = 0; i < m_trees.size(); ++ i)
    {
      Tree_item& tree = m_trees[i];

      Point_3 center
        = CGAL::barycenter (boost::make_transform_iterator
                            (tree.inliers().begin(), Point_from_iterator_and_pmap(m_point_map)),
                             boost::make_transform_iterator
                            (tree.inliers().end(), Point_from_iterator_and_pmap(m_point_map)));

      FT radius = FT(0);

      FT height = -std::numeric_limits<FT>::max();
      
      for (std::size_t j = 0; j < tree.inliers().size(); ++ j)
      {
        radius += CGAL::squared_distance (center, get (m_point_map, *(tree.inliers()[j])));
        height = (std::max)(height, get (m_point_map, *(tree.inliers()[j])).z());
      }

      radius = CGAL::approximate_sqrt (radius / tree.inliers().size());

      tree.set_center_radius_and_height (center, radius, height);

    }


    std::sort (m_trees.begin(), m_trees.end(),
               [](const Tree_item& a, const Tree_item& b) -> bool
               {
                 return a.radius() > b.radius();
               });
    std::vector<Tree_item> filtered_trees;

    for (std::size_t i = 0; i < m_trees.size(); ++ i)
    {
      Tree_item& tree_a = m_trees[i];

      if (tree_a.radius() < minimum_radius)
        continue;
      
      Circle_2 circle_a (Point_2 (tree_a.center().x(), tree_a.center().y()),
                         tree_a.radius() * tree_a.radius());

      bool okay = true;
      for (std::size_t j = 0; j < filtered_trees.size(); ++ j)
      {
        Tree_item& tree_b = filtered_trees[j];

        Circle_2 circle_b (Point_2 (tree_b.center().x(), tree_b.center().y()),
                           tree_b.radius() * tree_b.radius());
        if (CGAL::do_intersect (circle_a, circle_b) ||
            circle_b.has_on_bounded_side (circle_a.center()))
        {
          okay = false;
          break;
        }
      }

      if (okay)
        filtered_trees.push_back (tree_a);
    }

    m_trees.swap (filtered_trees);

    std::ofstream f ("trees.polylines.txt");
    f.precision(18);
    for (std::size_t i = 0; i < m_trees.size(); ++ i)
    {
      Tree_item& tree = m_trees[i];
    
      f << "13 ";
      for (std::size_t j = 0; j <= 12; ++ j)
      {
        FT angle = 2. * CGAL_PI * (j / FT(12));
        Point_3 p = tree.center() + Vector_3 (tree.radius() * std::cos(angle),
                                              tree.radius() * std::sin(angle),
                                              0.);
        f << p << " ";
      }
      f << std::endl;
    }
  }

  void fit_3d_models()
  {
    for (std::size_t i = 0; i < m_trees.size(); ++ i)
    {
      Tree_item& tree = m_trees[i];

      const Point_3& center = tree.center();
      Point_2 center_2 = tree.center_2();
      FT radius = tree.radius();
      FT height = tree.height();

      std::sort (tree.inliers().begin(), tree.inliers().end(),
                 [&](const const_iterator& a,
                     const const_iterator& b) -> bool
                 {
                   return get (m_point_map, *a).z() < get (m_point_map, *b).z();
                 });


      typename Tree_item::Model_3& model = tree.model();

      model.height[0] = get(m_point_map, *(tree.inliers().front())).z();
      model.height[1] = get(m_point_map, *(tree.inliers()[tree.inliers().size() / 10])).z();
      model.height[2] = get(m_point_map, *(tree.inliers()[tree.inliers().size() / 2])).z();
      model.height[3] = get(m_point_map, *(tree.inliers()[std::size_t(9. * tree.inliers().size() / 10.)])).z();
      model.height[4] = get(m_point_map, *(tree.inliers().back())).z();

      cpp11::array<FT, 4> width = make_array(FT(0),FT(0),FT(0),FT(0));
      cpp11::array<std::size_t, 4> nb = make_array(std::size_t(0),std::size_t(0),std::size_t(0),std::size_t(0));
      
      for (std::size_t j = 0; j < tree.inliers().size(); ++ j)
      {
        const Point_3& p = get(m_point_map, *(tree.inliers()[j]));
        Point_2 p2 (p.x(), p.y());

        std::size_t idx = 0;
        if (p.z() < model.height[1])
          idx = 0;
        else if (p.z() < model.height[2])
          idx = 1;
        else if (p.z() < model.height[3])
          idx = 2;
        else
          idx = 3;

        width[idx] += CGAL::squared_distance (p2, center_2);
        nb[idx] ++;
      }

      for (std::size_t j = 0; j < width.size(); ++ j)
        if (nb[j] != 0)
          width[j] = CGAL::approximate_sqrt(width[j] / nb[j]);

      model.width[0] = (width[0] + width[1]) / 2.;
      model.width[1] = (width[1] + width[2]) / 2.;
      model.width[2] = (width[2] + width[3]) / 2.;
    }
  }


};
  
} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TREE_ESTIMATOR_H
