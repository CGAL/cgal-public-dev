// Copyright (c) 2023
// INRIA Sophia-Antipolis (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef COLLISION_MESH_3_H
#define COLLISION_MESH_3_H

#include <CGAL/Surface_mesh.h>
#include <AABB_triangle_trajectory_primitive.h>

namespace CGAL {

template <class K_>
class Collision_mesh {

    public:
        typedef          K_                                 K;
        typedef typename K::Point_3                         Point;
        typedef typename K::Vector_3                        Vector;
        typedef          Surface_mesh<Point>                Base;
        typedef typename Base::Vertex_index                 Vertex_index;
        typedef typename Base::Edge_index                   Edge_index;
        typedef typename Base::Halfedge_index               Halfedge_index;
        typedef typename Base::Face_index                   Face_index;
        typedef typename Base::Face_range                   Face_range;
        typedef typename Base::Vertex_around_face_range     Vertex_around_face_range;
        typedef typename Base::size_type                    size_type;
        typedef          Triangle_trajectory<K, Face_index> Trajectory;

    public:
        Base & mesh_;
        decltype(mesh_.template add_property_map<Vertex_index, Vector>("v:velocity").first) vvelocity_;
        decltype(mesh_.template add_property_map<Vertex_index, Point>("v:next_point").first) vnext_point_;
        Collision_mesh(Base & mesh) : mesh_{mesh}
        {
            vvelocity_ = mesh_.template add_property_map<Vertex_index, Vector>("v:velocity").first;
            vnext_point_ = mesh_.template add_property_map<Vertex_index, Point>("v:next_point").first;
            for(Vertex_index vd : mesh_.vertices()){
                put(vvelocity_, vd, Vector(0, 0, 0));
                put(vnext_point_, vd, point(vd));
            }
        }

        const Vector& velocity(Vertex_index v) const;
        Vector& velocity(Vertex_index v);

        const Point& point(Vertex_index v) const;
        Point& point(Vertex_index v);

        const Point& next_point(Vertex_index v) const;
        Point& next_point(Vertex_index v);

        Face_range faces() const;
        Vertex_around_face_range vertices_around_face(Halfedge_index h) const;
        Halfedge_index halfedge(Face_index f) const;

        size_type num_faces() const;

};

    /// returns the velocity associated to vertex `v`.
    template <class K>
    const typename K::Vector_3& Collision_mesh<K>::velocity(Vertex_index v) const { return vvelocity_[v]; }

    template <class K>
    typename K::Vector_3& Collision_mesh<K>::velocity(Vertex_index v) { return vvelocity_[v]; }

    template <class K>
    const typename K::Point_3& Collision_mesh<K>::point(Vertex_index v) const { return mesh_.point(v); }

    template <class K>
    typename K::Point_3& Collision_mesh<K>::point(Vertex_index v) { return mesh_.point(v); }

    template <class K>
    const typename K::Point_3& Collision_mesh<K>::next_point(Vertex_index v) const { return vnext_point_[v]; }

    template <class K>
    typename K::Point_3& Collision_mesh<K>::next_point(Vertex_index v) { return vnext_point_[v]; }

    template <class K>
    typename Collision_mesh<K>::Face_range Collision_mesh<K>::faces() const { return mesh_.faces(); }

    template <class K>
    typename Collision_mesh<K>::Vertex_around_face_range Collision_mesh<K>::vertices_around_face(Halfedge_index h) const
    {
      return mesh_.vertices_around_face(h);
    }

    template <class K>
    typename Collision_mesh<K>::Halfedge_index Collision_mesh<K>::halfedge(Face_index f) const
    {
        return mesh_.halfedge(f);
    }

    template <class K>
    typename Collision_mesh<K>::size_type Collision_mesh<K>::num_faces() const
    {
        return mesh_.num_faces();
    }

} // end CGAL

#endif