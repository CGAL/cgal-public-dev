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

#ifndef CGAL_CLASSIFICATION_MESH_NEIGHBORHOOD_H
#define CGAL_CLASSIFICATION_MESH_NEIGHBORHOOD_H

#include <vector>

#include <boost/iterator/counting_iterator.hpp>

/// \cond SKIP_IN_MANUAL

namespace CGAL {

namespace Classification {

  /*!
    \ingroup PkgClassificationDataStructures

  */
template <typename FaceGraph>
class Mesh_neighborhood
{
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  const FaceGraph& m_mesh;
  
public:

  class One_ring_neighbor_query
  {
  public:
    typedef Mesh_neighborhood::face_descriptor value_type; ///<
  private:
    const Mesh_neighborhood& neighborhood;

  public:

    One_ring_neighbor_query (const Mesh_neighborhood& neighborhood)
      : neighborhood (neighborhood) { }

    /// \cond SKIP_IN_MANUAL
    template <typename OutputIterator>
    OutputIterator operator() (const value_type& query, OutputIterator output) const
    {
      neighborhood.one_ring_neighbors (query, output);
      return output;
    }
    /// \endcond
  };

  friend class One_ring_neighbor_query;


  Mesh_neighborhood (const FaceGraph& mesh) : m_mesh (mesh)
  {
  }


  /// \cond SKIP_IN_MANUAL
  ~Mesh_neighborhood ()
  {
  }
  /// \endcond

  One_ring_neighbor_query one_ring_neighbor_query () const
  {
    return One_ring_neighbor_query (*this);
  }


private:

  template <typename OutputIterator>
  void direct_neighbors (const face_descriptor& query, OutputIterator output) const
  {
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(query, m_mesh), m_mesh))
      {
        *(output ++ ) = face(opposite(hd, m_mesh), m_mesh);
      }
  }
  
  template <typename OutputIterator>
  void one_ring_neighbors (const face_descriptor& query, OutputIterator output) const
  {
    std::set<face_descriptor> done;
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(query, m_mesh), m_mesh))
      {
        BOOST_FOREACH(face_descriptor fd, faces_around_target(hd, m_mesh))
          {
            if (fd == boost::graph_traits<FaceGraph>::null_face())
              continue;
            if (done.insert(fd).second)
              *(output ++) = fd;
          }
      }
  }

};
  

}
  
}

/// \endcond

#endif // CGAL_CLASSIFICATION_MESH_NEIGHBORHOOD_H
