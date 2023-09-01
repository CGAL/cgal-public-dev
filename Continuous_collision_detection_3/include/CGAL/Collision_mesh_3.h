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
#include <CGAL/Trajectories.h>

namespace CGAL {



  template <typename K_>
  class Collision_mesh: public Surface_mesh<typename K_::Point_3> {

    public:
        struct Point_3_trajectory;
        struct Segment_3_trajectory;
        struct Triangle_3_trajectory;

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

        typedef            ::CGAL::Point_3_trajectory<K>    Point_trajectory;
        typedef            ::CGAL::Segment_3_trajectory<K>  Segment_trajectory;
        typedef            ::CGAL::Triangle_3_trajectory<K> Triangle_trajectory;

        // typedef decltype(Base().template add_property_map<Vertex_index, Vector>("v:velocity").first)        Vector_map;
        typedef decltype(Base().template add_property_map<Vertex_index, Point>("v:next_point").first)       Point_map;

        typedef decltype(Base().template add_property_map<Vertex_index, CGAL::IO::Color>("v:color").first)  Vertex_color_map;
        typedef decltype(Base().template add_property_map<Edge_index, CGAL::IO::Color>("e:color").first)    Edge_color_map;
        typedef decltype(Base().template add_property_map<Face_index, CGAL::IO::Color>("f:color").first)    Face_color_map;

        typedef typename boost::graph_traits<Collision_mesh<K>>::vertex_descriptor      vertex_descriptor;
        typedef typename boost::graph_traits<Collision_mesh<K>>::edge_descriptor        edge_descriptor;
        typedef typename boost::graph_traits<Collision_mesh<K>>::halfedge_descriptor    halfedge_descriptor;
        typedef typename boost::graph_traits<Collision_mesh<K>>::face_descriptor        face_descriptor;

        using Base::point;
        using Base::source;
        using Base::target;
        using Base::halfedge;
        using Base::next;

    private:
        Vector_map vvelocity_;
        Point_map vnext_point_;

        Vertex_color_map vcolor_;
        Edge_color_map ecolor_;
        Face_color_map fcolor_;

    public:
        Collision_mesh();

        Collision_mesh(const Base& surface_mesh);
        Collision_mesh(Base&& surface_mesh);

        Collision_mesh(const Collision_mesh& collision_mesh);
        Collision_mesh(Collision_mesh&& collision_mesh);

        Collision_mesh<K>& operator=(const Collision_mesh<K>& rhs); // assigns `rhs` to `*this`. Performs a deep copy of all properties.
        Collision_mesh<K>& operator=(Collision_mesh<K>&& collision_mesh); // move assignment

        const Vector& velocity(Vertex_index v) const;
        Vector& velocity(Vertex_index v);

        const Point& next_point(Vertex_index v) const;
        Point& next_point(Vertex_index v);

        Point_trajectory point_trajectory(Vertex_index v) const;
        Segment_trajectory edge_trajectory(Halfedge_index h) const;
        Triangle_trajectory face_trajectory(Face_index f) const;

        void color(const Face_index& fi, CGAL::IO::Color c);
  };

    // CONSTRUCTORS
    // ============
    template <class K>
    Collision_mesh<K>::Collision_mesh() : Base() {}

    template <class K>
    Collision_mesh<K>::Collision_mesh(const Base& surface_mesh) : Base{surface_mesh} {
        vvelocity_ = this->template add_property_map<Vertex_index, Vector>("v:velocity").first;
        vnext_point_ = this->template add_property_map<Vertex_index, Point>("v:next_point").first;
        vcolor_ = this->template add_property_map<Vertex_index, CGAL::IO::Color>("v:color").first;

        //
        // Initialize vertex maps
        for(const vertex_descriptor& vd : this->vertices()){
            put(vcolor_, vd, CGAL::IO::black());
            put(vvelocity_, vd, Vector(0, 0, 0));
            put(vnext_point_, vd, this->point(vd));
        }

        fcolor_ = this-> template add_property_map<Face_index, CGAL::IO::Color>("f:color").first;
        for( const face_descriptor& fd : this->faces()){
            put(fcolor_, fd, CGAL::IO::white());
        }
    }

    template <class K>
    Collision_mesh<K>::Collision_mesh(Base&& surface_mesh) : Base{std::move(surface_mesh)} {
        vvelocity_ = this-> template add_property_map<Vertex_index, Vector>("v:velocity").first;
        vnext_point_ = this-> template add_property_map<Vertex_index, Point>("v:next_point").first;
        vcolor_ = this-> template add_property_map<Vertex_index, CGAL::IO::Color>("v:color").first;

        for(const vertex_descriptor& vd : this->vertices()){
            put(vvelocity_, vd, Vector(0, 0, 0));
            put(vnext_point_, vd, point(vd));
            put(vcolor_, vd, CGAL::IO::black());
        }

        fcolor_ = this-> template add_property_map<Face_index, CGAL::IO::Color>("f:color").first;
        for( const face_descriptor& fd : this->faces()){
            put(fcolor_, fd, CGAL::IO::white());
        }
    }

    template <class K>
    Collision_mesh<K>::Collision_mesh(const Collision_mesh<K>& collision_mesh){ *this = collision_mesh; }

    template <class K>
    Collision_mesh<K>::Collision_mesh(Collision_mesh<K>&& collision_mesh)
        : Base{std::move(collision_mesh)}
        , vnext_point_(std::move(collision_mesh.vnext_point_))
        , vvelocity_(std::move(collision_mesh.vvelocity_))
        , vcolor_(std::move(collision_mesh.vcolor_))
        , fcolor_(std::move(collision_mesh.fcolor_))
    {}


    // OPERATORS
    // =========
    template <class K>
    Collision_mesh<K>& Collision_mesh<K>::operator=(const Collision_mesh<K>& rhs)
    {
        Base::operator=(rhs);
        if (this != &rhs)
        {
            vnext_point_ = this-> template property_map<Vertex_index, Point>("v:next_point").first;
            vvelocity_   = this-> template property_map<Vertex_index, Vector>("v:velocity").first;
            vcolor_      = this-> template property_map<Vertex_index, CGAL::IO::Color>("v:color").first;
            fcolor_      = this-> template property_map<Face_index, CGAL::IO::Color>("f:color").first;
        }
        return *this;
    }

    template <class K>
    Collision_mesh<K>& Collision_mesh<K>::operator=(Collision_mesh<K>&& collision_mesh)
    {
        Base::operator=(std::move(collision_mesh));

        vnext_point_    = std::move(collision_mesh.vnext_point_);
        vvelocity_      = std::move(collision_mesh.vvelocity_);
        vcolor_         = std::move(collision_mesh.vcolor_);
        fcolor_         = std::move(collision_mesh.fcolor_);

        return *this;
    }


    // MEMBER ROUTINES
    // ===============
    template <class K>
    const typename K::Vector_3& Collision_mesh<K>::velocity(Vertex_index v) const { return vvelocity_[v]; }

    template <class K>
    typename K::Vector_3& Collision_mesh<K>::velocity(Vertex_index v) { return vvelocity_[v]; }

    template <class K>
    const typename K::Point_3& Collision_mesh<K>::next_point(Vertex_index v) const { return vnext_point_[v]; }

    template <class K>
    typename K::Point_3& Collision_mesh<K>::next_point(Vertex_index v) { return vnext_point_[v]; }

    template <class K>
    auto Collision_mesh<K>::point_trajectory(Vertex_index v) const -> Point_trajectory
    {
        return Point_trajectory(point(v), next_point(v));
    }

    template <class K>
    auto Collision_mesh<K>::edge_trajectory(Halfedge_index h) const -> Segment_trajectory
    {
        return Segment_trajectory(
            point_trajectory(source(h)),
            point_trajectory(target(h))
        );
    }

    template <class K>
    auto Collision_mesh<K>::face_trajectory(Face_index f) const -> Triangle_trajectory
    {
        Halfedge_index h = halfedge(f);
        return Triangle_trajectory(
            point_trajectory(source(h)),
            point_trajectory(target(h)),
            point_trajectory(target(next(h)))
        );
    }

    template <class K>
    void Collision_mesh<K>::color(const Face_index& fi, CGAL::IO::Color c) {
        put(fcolor_, fi, c);
        return;
    }



} // namespace CGAL

#define CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS typename K
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME CGAL::Collision_mesh<K>
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<typename K::Point_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

#endif