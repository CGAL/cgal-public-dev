// Copyright (c) 2023 GeometryFactory (France).
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef COLLISION_MESH_3_H
#define COLLISION_MESH_3_H

#include <CGAL/Surface_mesh.h>
#include <CGAL/Trajectories.h>

namespace CGAL {


/// @{
/// \ingroup PkgCollisions3Classes

/// @brief The class `Collision_mesh` derives from `Surface_mesh`. It's primary purpose 
/// is to formalize the requirement that both the current and next positions of 
/// each point are known. This requirement enables the use of collision detection
/// algorithms.
template <typename Kernel>
class Collision_mesh : public Surface_mesh<typename Kernel::Point_3> {

    // TODO: Change template to use Point_3 and align with Surface_mesh.
    public:

        /// \name Types
        /// @{

        /// @brief The underlying Kernel class.
        typedef          Kernel                             K;

        /// @brief The type of `Point_3`.
        typedef typename K::Point_3                         Point;

        /// @brief The type of `Vector_3`.
        typedef typename K::Vector_3                        Vector;

        /// @brief The type of `Surface_mesh` from which this class derives.
        typedef          Surface_mesh<Point>                Base;

        /// @brief An alias for `Surface_mesh::Vertex_index`.
        typedef typename Base::Vertex_index                 Vertex_index;

        /// @brief An alias for `Surface_mesh::Edge_index`.
        typedef typename Base::Edge_index                   Edge_index;

        /// @brief An alias for `Surface_mesh::Halfedge_index`.
        typedef typename Base::Halfedge_index               Halfedge_index;

        /// @brief An alias for `Surface_mesh::Face_index`.
        typedef typename Base::Face_index                   Face_index;

        /// @brief An alias for `Surface_mesh::size_type`.
        typedef typename Base::size_type                    size_type;

        /// @brief The type of `Point_3_trajectory`.
        /// @details Contains current and next positions of a `Point_3` object.
        typedef          ::CGAL::Point_3_trajectory<K>      Point_trajectory;

        /// @brief The type of `Triangle_3_trajectory`.
        /// @details Contains current and next positions of a `Triangle_3` object.
        typedef          ::CGAL::Segment_3_trajectory<K>    Segment_trajectory;

        /// @brief The type of `Segment_3_trajectory`.
        /// @details Contains current and next positions of a `Segment_3` object.
        typedef          ::CGAL::Triangle_3_trajectory<K>   Triangle_trajectory;

        /// @}

    private:
        typedef decltype(Base().template add_property_map<Vertex_index, Point>("v:next_point").first)       Point_map;
        typedef decltype(Base().template add_property_map<Vertex_index, CGAL::IO::Color>("v:color").first)  Vertex_color_map;
        typedef decltype(Base().template add_property_map<Edge_index, CGAL::IO::Color>("e:color").first)    Edge_color_map;
        typedef decltype(Base().template add_property_map<Face_index, CGAL::IO::Color>("f:color").first)    Face_color_map;

        typedef typename boost::graph_traits< Collision_mesh<K> >::vertex_descriptor      vertex_descriptor;
        typedef typename boost::graph_traits< Collision_mesh<K> >::edge_descriptor        edge_descriptor;
        typedef typename boost::graph_traits< Collision_mesh<K> >::halfedge_descriptor    halfedge_descriptor;
        typedef typename boost::graph_traits< Collision_mesh<K> >::face_descriptor        face_descriptor;

        Point_map vnext_point_;
        Vertex_color_map vcolor_;
        Edge_color_map ecolor_;
        Face_color_map fcolor_;

    public:
        /// \name Creation
        /// @{

        /// @brief Creates an empty collision mesh
        Collision_mesh();

        /// @brief Specializes a surface mesh to a collision mesh.
        /// @details The next position of each point defaults to its current position.
        Collision_mesh(const Base& surface_mesh);

        /// @brief Specializes a surface mesh to a collision mesh.
        /// @details The next position of each point defaults to its current position.
        Collision_mesh(Base&& surface_mesh);

        /// @brief Copy constructor.
        Collision_mesh(const Collision_mesh& collision_mesh);

        /// @brief Move constructor.
        Collision_mesh(Collision_mesh&& collision_mesh);

        /// @}


        /// \name Operators
        /// @{

        /// @brief Copy assignment.
        Collision_mesh<K>& operator=(const Collision_mesh<K>& rhs); // assigns `rhs` to `*this`. Performs a deep copy of all properties.
        
        /// @brief Move assignment.
        Collision_mesh<K>& operator=(Collision_mesh<K>&& collision_mesh); // move assignment

        /// @}


        /// \name Methods
        /// @{

        /// @brief Returns the velocity vector computed as the difference between
        /// the current and next position of the vertex associated with that vertex index.
        const Vector& velocity(Vertex_index v) const;

        /// @brief Returns the velocity vector computed as the difference between
        /// the current and next position of the vertex associated with that vertex index.
        Vector& velocity(Vertex_index v);

        /// @brief Returns the next position associated with the vertex index.
        const Point& next_point(Vertex_index v) const;

        /// @brief Returns the next position associated with the vertex index.
        Point& next_point(Vertex_index v);

        /// @brief Returns the `Point_trajectory_3` associated with the vertex index.
        Point_trajectory point_trajectory(Vertex_index v) const;

        /// @brief Returns `Segment_trajectory_3` associated with the halfedge index.
        Segment_trajectory edge_trajectory(Halfedge_index h) const;

        /// @brief Returns the `Triangle_trajectory_3` associated with the face index.
        Triangle_trajectory face_trajectory(Face_index f) const;

        /// @brief Assigns the specified color to the face associated with the face index.
        void color(const Face_index& fi, CGAL::IO::Color c);

        /// @}

        using Base::point;
        using Base::source;
        using Base::target;
        using Base::halfedge;
        using Base::next;
};

/// @}

// ============
// CONSTRUCTORS
// ============
template <class K>
Collision_mesh<K>::Collision_mesh() : Base() {}

template <class K>
Collision_mesh<K>::Collision_mesh(const Base& surface_mesh) : Base{surface_mesh} {
    vnext_point_ = this->template add_property_map<Vertex_index, Point>("v:next_point").first;
    vcolor_ = this->template add_property_map<Vertex_index, CGAL::IO::Color>("v:color").first;

    //
    // Initialize vertex maps
    for(const vertex_descriptor& vd : this->vertices()){
        put(vcolor_, vd, CGAL::IO::black());
        put(vnext_point_, vd, this->point(vd));
    }

    fcolor_ = this-> template add_property_map<Face_index, CGAL::IO::Color>("f:color").first;
    for( const face_descriptor& fd : this->faces()){
        put(fcolor_, fd, CGAL::IO::white());
    }
}

template <class K>
Collision_mesh<K>::Collision_mesh(Base&& surface_mesh) : Base{std::move(surface_mesh)} {
    vnext_point_ = this-> template add_property_map<Vertex_index, Point>("v:next_point").first;
    vcolor_ = this-> template add_property_map<Vertex_index, CGAL::IO::Color>("v:color").first;

    for(const vertex_descriptor& vd : this->vertices()){
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
    , vcolor_(std::move(collision_mesh.vcolor_))
    , fcolor_(std::move(collision_mesh.fcolor_))
{}


// =========
// OPERATORS
// =========
template <class K>
Collision_mesh<K>& Collision_mesh<K>::operator=(const Collision_mesh<K>& rhs)
{
    Base::operator=(rhs);
    if (this != &rhs)
    {
        vnext_point_ = this-> template property_map<Vertex_index, Point>("v:next_point").first;
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
    vcolor_         = std::move(collision_mesh.vcolor_);
    fcolor_         = std::move(collision_mesh.fcolor_);

    return *this;
}


// ===============
// MEMBER ROUTINES
// ===============
template <class K>
const typename K::Vector_3& Collision_mesh<K>::velocity(Vertex_index v) const { 
    return vnext_point_[v] - vpoint_[v]; 
}

template <class K>
typename K::Vector_3& Collision_mesh<K>::velocity(Vertex_index v)  { 
    return vnext_point_[v] - vpoint_[v]; 
}

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