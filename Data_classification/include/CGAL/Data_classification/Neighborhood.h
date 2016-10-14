#ifndef CGAL_DATA_CLASSIFICATION_NEIGHBORHOOD_H
#define CGAL_DATA_CLASSIFICATION_NEIGHBORHOOD_H

#include <vector>

#include <boost/iterator/counting_iterator.hpp>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/centroid.h>
#include <CGAL/grid_simplify_point_set.h>

#include <CGAL/Data_classification/Image.h>

namespace CGAL {

namespace Data_classification {

  /*!
    \ingroup PkgDataClassification

    \brief Class that precomputes spatial searching structures and
    gives access to local neighborhoods of points.

    \tparam Kernel The geometric kernel used.
    \tparam RandomAccessIterator Iterator over the input.
    \tparam PointPMap Property map to access the input points.
  */
template <typename Kernel, typename RandomAccessIterator, typename PointPMap>
class Neighborhood
{
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;
  
  class My_point_property_map{
    RandomAccessIterator begin;
    PointPMap point_pmap;
    
  public:
    typedef Point value_type;
    typedef const value_type& reference;
    typedef std::size_t key_type;
    typedef boost::lvalue_property_map_tag category;
    My_point_property_map () { }
    My_point_property_map (RandomAccessIterator begin, PointPMap point_pmap)
      : begin (begin), point_pmap (point_pmap) { }
    reference operator[] (key_type k) const { return get(point_pmap, begin[k]); }
    friend inline reference get (const My_point_property_map& ppmap, key_type i) 
    { return ppmap[i]; }
  };

  typedef Search_traits_3<Kernel> SearchTraits_3;
  typedef Search_traits_adapter <std::size_t, My_point_property_map, SearchTraits_3> Search_traits;
  typedef Sliding_midpoint<Search_traits> Splitter;
  typedef Distance_adapter<std::size_t, My_point_property_map, Euclidean_distance<SearchTraits_3> > Distance;
  typedef Kd_tree<Search_traits, Splitter, Tag_true> Tree;
  typedef Fuzzy_sphere<Search_traits> Sphere;
  typedef Orthogonal_k_neighbor_search<Search_traits, Distance, Splitter, Tree> Knn;


  Tree* m_tree;
  Distance m_distance;

  std::vector<std::vector<std::size_t> > m_precomputed_neighbors;
  
public:
  
  Neighborhood () : m_tree (NULL) { }

  /*!
    \brief Constructs a neighborhood object based on the input range.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_pmap Property map to access the input points
  */
  Neighborhood (const RandomAccessIterator& begin,
                const RandomAccessIterator& end,
                PointPMap point_pmap)
    : m_tree (NULL)
  {
    std::size_t size = end - begin;
    
    My_point_property_map pmap (begin, point_pmap);
    m_tree = new Tree (boost::counting_iterator<std::size_t> (0),
                       boost::counting_iterator<std::size_t> (size),
                       Splitter(),
                       Search_traits (pmap));
    m_distance = Distance (pmap);
    m_tree->build();
  }

  /*!
    \brief Constructs a simplified neighborhood object based on the input range.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_pmap Property map to access the input points
  */
  Neighborhood (const RandomAccessIterator& begin,
                const RandomAccessIterator& end,
                PointPMap point_pmap,
                double voxel_size)
    : m_tree (NULL)
  {
    // First, simplify
    std::size_t size = end - begin;
    std::vector<std::size_t> indices (size);
    for (std::size_t i = 0; i < indices.size(); ++ i)
      indices[i] = i;
    My_point_property_map pmap (begin, point_pmap);

    voxelize_point_set(indices, pmap, voxel_size);
    
    m_tree = new Tree (indices.begin(), indices.end(),
                       Splitter(),
                       Search_traits (pmap));
    m_distance = Distance (pmap);
    m_tree->build();
  }

  /// \cond SKIP_IN_MANUAL
  ~Neighborhood ()
  {
    if (m_tree != NULL)
      delete m_tree;
  }
  /// \endcond

  /*!
    \brief Gets the nearest neighbors computed in a local sphere of user defined radius.

    \tparam OutputIterator Where the indices of found neighbor points are stored.
    \param query The query point.
    \param radius_neighbors Radius of the query sphere.
  */
  template <typename OutputIterator>
  void range_neighbors (const Point& query, const FT radius_neighbors, OutputIterator output) const
  {
    CGAL_assertion (m_tree != NULL);
    Sphere fs (query, radius_neighbors, 0, m_tree->traits());
    m_tree->search (output, fs);
  }

  /*!
    \brief Gets the K nearest neighbors.

    \tparam OutputIterator Where the indices of found neighbor points are stored.
    \param index Index of the query point.
    \param k Number of nearest neighbors.
  */
  template <typename OutputIterator>
  void k_neighbors (const Point& query, const std::size_t k, OutputIterator output) const
  {
    CGAL_assertion (m_tree != NULL);
    Knn search (*m_tree, query, k, 0, true, m_distance);
    for (typename Knn::iterator it = search.begin(); it != search.end(); ++ it)
      *(output ++) = it->first;
  }

private:

  template <typename PointMap>
  void voxelize_point_set (std::vector<std::size_t>& indices, PointMap point_map,
                           double voxel_size)
  {
    std::map<Point, std::vector<std::size_t> > grid;

    for (std::size_t i = 0; i < indices.size(); ++ i)
      {
        const Point& p = get(point_map, indices[i]);
        Point ref (std::floor(p.x() / voxel_size),
                   std::floor(p.y() / voxel_size),
                   std::floor(p.z() / voxel_size));
        typename std::map<Point, std::vector<std::size_t> >::iterator it;
        boost::tie (it, boost::tuples::ignore)
          = grid.insert (std::make_pair (ref, std::vector<std::size_t>()));
        it->second.push_back (indices[i]);
      }
    indices.clear();
    for (typename std::map<Point, std::vector<std::size_t> >::iterator
           it = grid.begin(); it != grid.end(); ++ it)
      {
        const std::vector<std::size_t>& pts = it->second;
        Point centroid = CGAL::centroid (boost::make_transform_iterator
                                         (pts.begin(),
                                          CGAL::Property_map_to_unary_function<PointMap>(point_map)),
                                         boost::make_transform_iterator
                                         (pts.end(),
                                          CGAL::Property_map_to_unary_function<PointMap>(point_map)));
        std::size_t chosen = 0;
        double min_dist = std::numeric_limits<double>::max();
        for (std::size_t i = 0; i < pts.size(); ++ i)
          {
            double dist = CGAL::squared_distance(get(point_map, pts[i]), centroid);
            if (dist < min_dist)
              {
                min_dist = dist;
                chosen = pts[i];
              }
          }
        indices.push_back (chosen);
      }
  }
};
  

}
  
}


#endif // CGAL_DATA_CLASSIFICATION_NEIGHBORHOOD_H
