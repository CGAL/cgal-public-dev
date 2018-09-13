#ifndef CGAL_LEVEL_OF_DETAIL_VEGETATION_SEGMENTOR_H
#define CGAL_LEVEL_OF_DETAIL_VEGETATION_SEGMENTOR_H

#include <CGAL/Level_of_detail/internal/Vegetation/Tree.h>

#include <CGAL/Classification/Image.h>

#include <boost/make_shared.hpp>

// TO REMOVE
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>


namespace CGAL {
namespace Level_of_detail {

template<typename GeomTraits, typename InputRange, typename PointMap>
class Vegetation_segmentor
{
public:

  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename InputRange::const_iterator const_iterator;

  typedef Tree<GeomTraits, InputRange, PointMap> Tree_item;
  
  typedef std::vector<const_iterator> Iterators;
  typedef std::pair<std::size_t, std::size_t> Coord;

  typedef CGAL::Classification::Image<Iterators> Image;
private:

  class Persistent_component
  {
    std::vector<Coord> m_inliers;
    FT m_starting_height;
    FT m_ending_height;
  public:

    Persistent_component (const Coord& coord, FT height)
      : m_inliers (1, coord), m_starting_height (height)
    {
    }

    std::vector<Coord>& inliers() { return m_inliers; }

    void add (const Coord& coord)
    {
      m_inliers.push_back (coord);
    }

    const FT& starting_height() const { return m_starting_height; }
    FT& ending_height() { return m_ending_height; }

    FT size() const { return m_starting_height - m_ending_height; }
    
  };

  typedef boost::shared_ptr<Persistent_component> Persistent_component_ptr;
  typedef CGAL::Classification::Image<Persistent_component_ptr> Tree_image;
    
  InputRange m_input;
  PointMap m_point_map;

public:

  Vegetation_segmentor (const InputRange& input, PointMap point_map)
    : m_input (input), m_point_map (point_map)
  {
  }

  void segment (const FT& grid_cell_size, const FT& minimum_height,
                std::vector<Tree_item>& output)
  {
    CGAL::Bbox_3 bbox
      = CGAL::bbox_3 (boost::make_transform_iterator
                      (m_input.begin(), Property_map_to_unary_function<PointMap>(m_point_map)),
                      boost::make_transform_iterator
                      (m_input.end(), Property_map_to_unary_function<PointMap>(m_point_map)));

    std::size_t width = (std::size_t)((bbox.xmax() - bbox.xmin()) / grid_cell_size) + 1;
    std::size_t height = (std::size_t)((bbox.ymax() - bbox.ymin()) / grid_cell_size) + 1;

    Image image (width, height);

    for (const_iterator ce_it = m_input.begin();
         ce_it != m_input.end(); ++ce_it)
    {
      const Point_3 &point = get(m_point_map, *ce_it);
      std::size_t x = (std::size_t)((point.x() - bbox.xmin()) / grid_cell_size);
      std::size_t y = (std::size_t)((point.y() - bbox.ymin()) / grid_cell_size);
      image(x,y).push_back (ce_it);
    }

    // Sort cells by height
    std::vector<Coord> sorted;
    for (std::size_t x = 0; x < image.width(); ++ x)
      for (std::size_t y = 0; y < image.height(); ++ y)
        if (!image(x,y).empty())
        {
          std::sort (image(x,y).begin(), image(x,y).end(),
                     [&](const const_iterator& a, const const_iterator& b) -> bool
                     {
                       return (get(m_point_map, *a).z()
                               > get(m_point_map, *b).z());
                     });
          sorted.push_back (std::make_pair(x,y));
        }

    std::sort (sorted.begin(), sorted.end(),
               [&](const Coord& a, const Coord& b) -> bool
               {
                 return (get(m_point_map, *(image(a.first, a.second).front())).z()
                         > get(m_point_map, *(image(b.first, b.second).front())).z());
               });

    // Init tree map
    Tree_image tree_map (width, height);
    for (std::size_t x = 0; x < image.width(); ++ x)
      for (std::size_t y = 0; y < image.height(); ++ y)
        tree_map(x,y) = Persistent_component_ptr();

    std::vector<Persistent_component_ptr> components;
        
    // pick cells one by one and track life time of each component
    for (std::size_t i = 0; i < sorted.size(); ++ i)
    {
      std::size_t x = sorted[i].first;
      std::size_t y = sorted[i].second;

      if (tree_map(x,y) != Persistent_component_ptr()) // Already handled
        continue;
          
      Iterators& cell = image(x,y);

      std::set<Persistent_component_ptr> local_trees;

      int xmin = int(x) - 1; if (xmin < 0) xmin = 0;
      int xmax = int(x) + 1; if (xmax == image.width())  xmax = image.width() - 1;
      int ymin = int(y) - 1; if (ymin < 0) ymin = 0;
      int ymax = int(y) + 1; if (ymax == image.height()) ymax = image.height() - 1;
          
      for (int xx = xmin; xx <= xmax; ++ xx)
        for (int yy = ymin; yy <= ymax; ++ yy)
        {
          if (xx == x && yy == y) // same cell
            continue;

          Iterators& neighbor = image(xx,yy);
          if (neighbor.empty())
            continue;

          if (tree_map(xx, yy) != Persistent_component_ptr()) // not unused cell
            local_trees.insert (tree_map(xx, yy));
        }

      if (local_trees.empty())
      {
        components.push_back
          (boost::make_shared<Persistent_component>
           (sorted[i], get(m_point_map, *(image(x,y).front())).z()));

        tree_map(x,y) = components.back();
      }
      else if (local_trees.size() == 1)
      {
        Persistent_component_ptr local_tree = *(local_trees.begin());
        tree_map(x,y) = local_tree;
        local_tree->add (sorted[i]);
      }
      else // merge happens
      {

        FT height = get(m_point_map, *(image(x,y).front())).z();
        Persistent_component_ptr chosen;
        FT size_max = -std::numeric_limits<FT>::max();

        // Keep highest component
        for (typename std::set<Persistent_component_ptr>::iterator it = local_trees.begin();
             it != local_trees.end(); ++ it)
        {
          Persistent_component_ptr neighbor = *it;

          neighbor->ending_height() = height;
          
          if (neighbor->size() > size_max)
          {
            size_max = neighbor->size();
            chosen = neighbor;
          }
        }

        // Add current cell
        chosen->add (sorted[i]);
        tree_map(x,y) = chosen;
        
        // Merge other components
        for (typename std::set<Persistent_component_ptr>::iterator it = local_trees.begin();
             it != local_trees.end(); ++ it)
        {
          Persistent_component_ptr neighbor = *it;
          if (neighbor == chosen)
            continue;

          // component with size above threshold are trees
          if (neighbor->size() < minimum_height)
          {
            for (std::size_t n = 0; n < neighbor->inliers().size(); ++ n)
            {
              const Coord& inlier = neighbor->inliers()[n];
              tree_map(inlier.first, inlier.second) = chosen;
              chosen->add (inlier);
            }
            neighbor->inliers().clear();
          }
        }
      }
    }
    
    for (std::size_t i = 0; i < components.size(); ++ i)
    {
      if (components[i]->inliers().empty())
        continue;

      output.push_back (Tree_item());

      for (std::size_t j = 0; j < components[i]->inliers().size(); ++ j)
      {
        const Coord& coord = components[i]->inliers()[j];
        const Iterators& cell = image(coord.first, coord.second);
        std::copy (cell.begin(), cell.end(), std::back_inserter (output.back().inliers()));
      }
    }
  }

private:

};
  
} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_VEGETATION_SEGMENTOR_H
