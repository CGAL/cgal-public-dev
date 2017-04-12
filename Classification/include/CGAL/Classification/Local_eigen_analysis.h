// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
#define CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H

#include <vector>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/centroid.h>
#include <CGAL/PCA_util.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#include <tbb/mutex.h>
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassificationDataStructures

    \brief Class that precomputes and stored the eigenvectors and
    eigenvalues of the covariance matrices of all points of a point
    set using a local neighborhood.

  */
class Local_eigen_analysis
{
public:
  typedef CGAL::cpp11::array<float, 3> Eigenvalues; ///< Eigenvalues (sorted in ascending order)
  
private:

#ifdef CGAL_LINKED_WITH_TBB
  template <typename PointRange, typename PointMap, typename NeighborQuery, typename DiagonalizeTraits>
  class Compute_eigen_values
  {
    Local_eigen_analysis& m_eigen;
    const PointRange& m_input;
    PointMap m_point_map;
    const NeighborQuery& m_neighbor_query;
    float& m_mean_range;
    tbb::mutex& m_mutex;
    
  public:
    
    Compute_eigen_values (Local_eigen_analysis& eigen,
                          const PointRange& input,
                          PointMap point_map,
                          const NeighborQuery& neighbor_query,
                          float& mean_range,
                          tbb::mutex& mutex)
      : m_eigen (eigen), m_input (input), m_point_map (point_map),
        m_neighbor_query (neighbor_query), m_mean_range (mean_range), m_mutex (mutex)
    { }
    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t i = r.begin(); i != r.end(); ++ i)
      {
        std::vector<std::size_t> neighbors;
        m_neighbor_query (get(m_point_map, *(m_input.begin()+i)), std::back_inserter (neighbors));

        std::vector<typename PointMap::value_type> neighbor_points;
        for (std::size_t j = 0; j < neighbors.size(); ++ j)
          neighbor_points.push_back (get(m_point_map, *(m_input.begin()+neighbors[j])));

        m_mutex.lock();
        m_mean_range += CGAL::sqrt (CGAL::squared_distance (get(m_point_map, *(m_input.begin() + i)),
                                                            get(m_point_map, *(m_input.begin() + neighbors.back()))));
        m_mutex.unlock();
          
        m_eigen.compute<typename PointMap::value_type,
                        DiagonalizeTraits> (i, get(m_point_map, *(m_input.begin()+i)), neighbor_points);
      }
    }

  };

  template <typename FaceListGraph, typename NeighborQuery, typename DiagonalizeTraits>
  class Compute_eigen_values_graph
  {
    typedef typename boost::graph_traits<FaceListGraph>::face_descriptor face_descriptor;
    
    Local_eigen_analysis& m_eigen;
    const FaceListGraph& m_input;
    const NeighborQuery& m_neighbor_query;
    float& m_mean_range;
    tbb::mutex& m_mutex;
    
  public:
    
    Compute_eigen_values_graph (Local_eigen_analysis& eigen,
                                const FaceListGraph& input,
                                const NeighborQuery& neighbor_query,
                                float& mean_range,
                                tbb::mutex& mutex)
      : m_eigen (eigen), m_input (input),
        m_neighbor_query (neighbor_query), m_mean_range (mean_range), m_mutex (mutex)
    { }
    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for (std::size_t i = r.begin(); i != r.end(); ++ i)
      {
        face_descriptor fd(i);
        std::vector<face_descriptor> neighbors;
        m_neighbor_query (fd, std::back_inserter (neighbors));

        m_mutex.lock();
        m_mean_range += m_eigen.face_radius(fd, m_input);
        m_mutex.unlock();
        
        m_eigen.compute_triangles<FaceListGraph, DiagonalizeTraits>
          (m_input, fd, neighbors);
      }
    }

  };
#endif

  typedef CGAL::cpp11::array<float, 3> float3;
  std::vector<float3> m_eigenvalues;
  std::vector<float> m_sum_eigenvalues;
  std::vector<float3> m_centroids;
  std::vector<float3> m_smallest_eigenvectors;
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
  std::vector<float3> m_middle_eigenvectors;
  std::vector<float3> m_largest_eigenvectors;
#endif
  float m_mean_range;

  
public:

  /// \cond SKIP_IN_MANUAL
  Local_eigen_analysis () { }
  /// \endcond

  /*! 
    \brief Computes the local eigen analysis of an input range based
    on a local neighborhood.

    \tparam PointRange model of `ConstRange`. Its iterator type is
    `RandomAccessIterator`.
    \tparam PointMap model of `ReadablePropertyMap` whose key
    type is the value type of the iterator of `PointRange` and value type
    is `Geom_traits::Point_3`.
    \tparam NeighborQuery model of `NeighborQuery`
    \tparam ConcurrencyTag enables sequential versus parallel
    algorithm. Possible values are `Parallel_tag` (default value is %CGAL
    is linked with TBB) or `Sequential_tag` (default value otherwise).
    \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
    for matrix diagonalization.

    \param input input range.
    \param point_map property map to access the input points
    \param neighbor_query object used to access neighborhoods of points.
  */
  template <typename PointRange,
            typename PointMap,
            typename NeighborQuery,
#if defined(DOXYGEN_RUNNING)
            typename ConcurrencyTag,
#elif defined(CGAL_LINKED_WITH_TBB)
            typename ConcurrencyTag = CGAL::Parallel_tag,
#else
            typename ConcurrencyTag = CGAL::Sequential_tag,
#endif
            typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<float, 3> >
    Local_eigen_analysis (const PointRange& input,
                          PointMap point_map,
                          const NeighborQuery& neighbor_query,
                          const ConcurrencyTag& = ConcurrencyTag(),
                          const DiagonalizeTraits& = DiagonalizeTraits())
  {
    m_eigenvalues.resize (input.size());
    m_sum_eigenvalues.resize (input.size());
    m_centroids.resize (input.size());
    m_smallest_eigenvectors.resize (input.size());
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
    m_middle_eigenvectors.resize (input.size());
    m_largest_eigenvectors.resize (input.size());
#endif
    
    m_mean_range = 0.;
      
#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
    {
      tbb::mutex mutex;
      Compute_eigen_values<PointRange, PointMap, NeighborQuery, DiagonalizeTraits>
        f(*this, input, point_map, neighbor_query, m_mean_range, mutex);
      tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size ()), f);
    }
    else
#endif
    {
      for (std::size_t i = 0; i < input.size(); i++)
      {
        std::vector<std::size_t> neighbors;
        neighbor_query (get(point_map, *(input.begin()+i)), std::back_inserter (neighbors));

        std::vector<typename PointMap::value_type> neighbor_points;
        for (std::size_t j = 0; j < neighbors.size(); ++ j)
          neighbor_points.push_back (get(point_map, *(input.begin()+neighbors[j])));

        m_mean_range += CGAL::sqrt (CGAL::squared_distance (get(point_map, *(input.begin() + i)),
                                                            get(point_map, *(input.begin() + neighbors.back()))));
        
        compute<typename PointMap::value_type, DiagonalizeTraits>
          (i, get(point_map, *(input.begin()+i)), neighbor_points);
      }
    }
    m_mean_range /= input.size();
  }


#if defined(DOXYGEN_RUNNING)
  template <typename FaceListGraph,
            typename NeighborQuery,
            typename ConcurrencyTag,
            typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<float, 3> >
  Local_eigen_analysis (const FaceListGraph& input,
                        const NeighborQuery& neighbor_query,
                        const ConcurrencyTag& = ConcurrencyTag(),
                        const DiagonalizeTraits& = DiagonalizeTraits())
  { }
#endif

  /// \cond SKIP_IN_MANUAL
  
  // To remove the ambiguity between this constructor and the point
  // set based one, we explicitly make one constructor with
  // CGAL::Sequential_tag and one with CGAL::Parallel_tag (the
  // parallel version is default if TBB is enable, sequential default
  // if it is not, exactly what is documented with the fake templated
  // version).
  
  template <typename FaceListGraph,
            typename NeighborQuery,
            typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<float, 3> >
  Local_eigen_analysis (const FaceListGraph& input,
                        const NeighborQuery& neighbor_query,
#ifdef CGAL_LINKED_WITH_TBB
                        const CGAL::Sequential_tag&,
#else
                        const CGAL::Sequential_tag& = CGAL::Sequential_tag&,
#endif
                        const DiagonalizeTraits& = DiagonalizeTraits())
  {
    typedef typename boost::graph_traits<FaceListGraph>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<FaceListGraph>::face_iterator face_iterator;
    typedef typename CGAL::Iterator_range<face_iterator> Face_range;

    Face_range range (faces(input));

    m_eigenvalues.resize (range.size());
    m_sum_eigenvalues.resize (range.size());
    m_centroids.resize (range.size());
    m_smallest_eigenvectors.resize (range.size());
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
    m_middle_eigenvectors.resize (range.size());
    m_largest_eigenvectors.resize (range.size());
#endif
    
    m_mean_range = 0.;
      
    BOOST_FOREACH(face_descriptor fd, range)
    {
      std::vector<face_descriptor> neighbors;
      neighbor_query (fd, std::back_inserter (neighbors));

      m_mean_range += face_radius(fd);
        
      compute_triangles<FaceListGraph, DiagonalizeTraits>
        (input, fd, neighbors);
    }

    m_mean_range /= range.size();
  }

#ifdef CGAL_LINKED_WITH_TBB
  template <typename FaceListGraph,
            typename NeighborQuery,
            typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<float, 3> >
  Local_eigen_analysis (const FaceListGraph& input,
                        const NeighborQuery& neighbor_query,
                        const CGAL::Parallel_tag& = CGAL::Parallel_tag(),
                        const DiagonalizeTraits& = DiagonalizeTraits())
  {
    typedef typename boost::graph_traits<FaceListGraph>::face_iterator face_iterator;
    typedef typename CGAL::Iterator_range<face_iterator> Face_range;
    
    Face_range range (faces(input));

    m_eigenvalues.resize (range.size());
    m_sum_eigenvalues.resize (range.size());
    m_centroids.resize (range.size());
    m_smallest_eigenvectors.resize (range.size());
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
    m_middle_eigenvectors.resize (range.size());
    m_largest_eigenvectors.resize (range.size());
#endif
    
    m_mean_range = 0.;

    tbb::mutex mutex;
    Compute_eigen_values_graph<FaceListGraph, NeighborQuery, DiagonalizeTraits>
        f(*this, input, neighbor_query, m_mean_range, mutex);

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, range.size()), f);
    m_mean_range /= range.size();
  }
#endif

  /// \endcond


  /*!
    \brief Returns the estimated unoriented normal vector of the point at position `index`.
    \tparam Geom_traits model of \cgal Kernel.
  */
  template <typename Geom_traits>
  typename Geom_traits::Vector_3 normal_vector (std::size_t index) const
  {
    return typename Geom_traits::Vector_3(double(m_smallest_eigenvectors[index][0]),
                                          double(m_smallest_eigenvectors[index][1]),
                                          double(m_smallest_eigenvectors[index][2]));
  }

  /*!
    \brief Returns the estimated local tangent plane of the point at position `index`.
    \tparam Geom_traits model of \cgal Kernel.
  */
  template <typename Geom_traits>
  typename Geom_traits::Plane_3 plane (std::size_t index) const
  {
    return typename Geom_traits::Plane_3
      (typename Geom_traits::Point_3 (double(m_centroids[index][0]),
                                      double(m_centroids[index][1]),
                                      double(m_centroids[index][2])),
       typename Geom_traits::Vector_3 (double(m_smallest_eigenvectors[index][0]),
                                       double(m_smallest_eigenvectors[index][1]),
                                       double(m_smallest_eigenvectors[index][2])));
  }

  /*!
    \brief Returns the normalized eigenvalues of the point at position `index`.
  */
  const Eigenvalues& eigenvalue (std::size_t index) const { return m_eigenvalues[index]; }

  /*!
    \brief Returns the sum of eigenvalues of the point at position `index`.
  */
  const float& sum_of_eigenvalues (std::size_t index) const { return m_sum_eigenvalues[index]; }

  /// \cond SKIP_IN_MANUAL
  float mean_range() const { return m_mean_range; }
  /// \endcond

private:

  template <typename FaceListGraph>
  float face_radius (typename boost::graph_traits<FaceListGraph>::face_descriptor& fd,
                     const FaceListGraph& g)
  {
    typedef typename boost::graph_traits<FaceListGraph>::halfedge_descriptor halfedge_descriptor;
    
    float out = 0.f;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, g), g))
    {
      out = (std::max)(out,
                       float(CGAL::squared_distance (get(get (CGAL::vertex_point, g), source(hd,g)),
                                                     get(get (CGAL::vertex_point, g), target(hd,g)))));
    }
    return out;
  }

  template <typename Point, typename DiagonalizeTraits>
  void compute (std::size_t index, const Point& query, std::vector<Point>& neighbor_points)
  {
    typedef typename Kernel_traits<Point>::Kernel::Vector_3 Vector;
      
    if (neighbor_points.size() == 0)
    {
      Eigenvalues v = {{ 0.f, 0.f, 0.f }};
      m_eigenvalues[index] = v;
      m_centroids[index] = {{ float(query.x()), float(query.y()), float(query.z()) }};
      m_smallest_eigenvectors[index] = {{ 0.f, 0.f, 1.f }};
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
      m_middle_eigenvectors[index] = {{ 0.f, 1.f, 0.f }}; 
      m_largest_eigenvectors[index] = {{ 1.f, 0.f, 0.f }};
#endif
      return;
    }

    Point centroid = CGAL::centroid (neighbor_points.begin(), neighbor_points.end());
    m_centroids[index] = {{ float(centroid.x()), float(centroid.y()), float(centroid.z()) }};
    
    CGAL::cpp11::array<float, 6> covariance = {{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }};
      
    for (std::size_t i = 0; i < neighbor_points.size(); ++ i)
    {
      Vector d = neighbor_points[i] - centroid;
      covariance[0] += d.x () * d.x ();
      covariance[1] += d.x () * d.y ();
      covariance[2] += d.x () * d.z ();
      covariance[3] += d.y () * d.y ();
      covariance[4] += d.y () * d.z ();
      covariance[5] += d.z () * d.z ();
    }

    CGAL::cpp11::array<float, 3> evalues = {{ 0.f, 0.f, 0.f }};
    CGAL::cpp11::array<float, 9> evectors = {{ 0.f, 0.f, 0.f,
                                               0.f, 0.f, 0.f,
                                               0.f, 0.f, 0.f }};

    DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
      (covariance, evalues, evectors);

    // Normalize
    float sum = evalues[0] + evalues[1] + evalues[2];
    if (sum > 0.f)
      for (std::size_t i = 0; i < 3; ++ i)
        evalues[i] = evalues[i] / sum;
    m_sum_eigenvalues[index] = float(sum);
    m_eigenvalues[index] = {{ float(evalues[0]), float(evalues[1]), float(evalues[2]) }};
    m_smallest_eigenvectors[index] = {{ float(evectors[0]), float(evectors[1]), float(evectors[2]) }};
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
    m_middle_eigenvectors[index] = {{ float(evectors[3]), float(evectors[4]), float(evectors[5]) }};
    m_largest_eigenvectors[index] = {{ float(evectors[6]), float(evectors[7]), float(evectors[8]) }};
#endif
  }

  template <typename FaceListGraph, typename DiagonalizeTraits>
  void compute_triangles (const FaceListGraph& g,
                          typename boost::graph_traits<FaceListGraph>::face_descriptor& query,
                          std::vector<typename boost::graph_traits<FaceListGraph>::face_descriptor>& neighbor_faces)
  {
    typedef typename boost::property_map<FaceListGraph, boost::vertex_point_t>::type::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    typedef typename Kernel::Triangle_3 Triangle;
        
    typedef typename boost::graph_traits<FaceListGraph>::face_descriptor face_descriptor;

    if (neighbor_faces.size() == 0)
    {
      Eigenvalues v = {{ 0.f, 0.f, 0.f }};
      m_eigenvalues[query] = v;

      CGAL::cpp11::array<Triangle,1> tr
        = {{ Triangle (get(get (CGAL::vertex_point, g), target(halfedge(query, g), g)),
                       get(get (CGAL::vertex_point, g), target(next(halfedge(query, g), g), g)),
                       get(get (CGAL::vertex_point, g), target(next(next(halfedge(query, g), g), g), g))) }};
      Point c = CGAL::centroid(tr.begin(),
                               tr.end(), Kernel(), CGAL::Dimension_tag<2>());

      m_centroids[query] = {{ float(c.x()), float(c.y()), float(c.z()) }};
      
      m_smallest_eigenvectors[query] = {{ 0.f, 0.f, 1.f }};
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
      m_middle_eigenvectors[query] = {{ 0.f, 1.f, 0.f }}; 
      m_largest_eigenvectors[query] = {{ 1.f, 0.f, 0.f }};
#endif
      return;
    }

    std::vector<Triangle> triangles;
    triangles.reserve(neighbor_faces.size());
    for (std::size_t i = 0; i < neighbor_faces.size(); ++ i)
    {
      const face_descriptor& fd = neighbor_faces[i];
      triangles.push_back
        (Triangle (get(get (CGAL::vertex_point, g), target(halfedge(fd, g), g)),
                   get(get (CGAL::vertex_point, g), target(next(halfedge(fd, g), g), g)),
                   get(get (CGAL::vertex_point, g), target(next(next(halfedge(fd, g), g), g), g))));
    }

    CGAL::cpp11::array<float, 6> covariance = {{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }};
    Point c = CGAL::centroid(triangles.begin(),
                             triangles.end(), Kernel(), CGAL::Dimension_tag<2>());
    
    CGAL::internal::assemble_covariance_matrix_3 (triangles.begin(), triangles.end(), covariance,
                                                  c, Kernel(), (Triangle*)NULL, CGAL::Dimension_tag<2>(),
                                                  DiagonalizeTraits());
      
    m_centroids[query] = {{ float(c.x()), float(c.y()), float(c.z()) }};
    
    CGAL::cpp11::array<float, 3> evalues = {{ 0.f, 0.f, 0.f }};
    CGAL::cpp11::array<float, 9> evectors = {{ 0.f, 0.f, 0.f,
                                               0.f, 0.f, 0.f,
                                               0.f, 0.f, 0.f }};

    DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
      (covariance, evalues, evectors);

    // Normalize
    float sum = evalues[0] + evalues[1] + evalues[2];
    if (sum > 0.f)
      for (std::size_t i = 0; i < 3; ++ i)
        evalues[i] = evalues[i] / sum;
    m_sum_eigenvalues[query] = float(sum);
    m_eigenvalues[query] = {{ float(evalues[0]), float(evalues[1]), float(evalues[2]) }};
    m_smallest_eigenvectors[query] = {{ float(evectors[0]), float(evectors[1]), float(evectors[2]) }};
#ifdef CGAL_CLASSIFICATION_EIGEN_FULL_STORAGE
    m_middle_eigenvectors[query] = {{ float(evectors[3]), float(evectors[4]), float(evectors[5]) }};
    m_largest_eigenvectors[query] = {{ float(evectors[6]), float(evectors[7]), float(evectors[8]) }};
#endif
  }

};
  

}
  
}


#endif // CGAL_CLASSIFICATION_LOCAL_EIGEN_ANALYSIS_H
