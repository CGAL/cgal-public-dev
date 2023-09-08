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

#ifndef COLLISION_SCENE_3_H
#define COLLISION_SCENE_3_H

#include <CGAL/Continuous_collision_detection_3/internal/AABB_triangle_trajectory_primitive.h>

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Trajectories.h>
#include <CGAL/Collision_mesh_3.h>

#include <vector>
#include <utility>


namespace CGAL {


/// \ingroup PkgCollisions3Classes
/// @{

/// @brief The class `Collision_scene` serves as an aggregator of `Collision_mesh` objects
/// to be evaluated for mutual collisions. 
/// @details The class provides standardized indexing across the primitives of all contained
/// collision meshes, as well as an AABB Tree that tracks the trajectories of these primitives. 
/// `Collision_scene` serves as the standard data structure for querying and 
/// visualizing collisions between an arbitrary number of collision meshes.
template <class Kernel>
class Collision_scene {

    public:


        /// \name Member Classes
        /// @{

        /// @brief a simple index wrapper provided to enable pretty-printing of mesh indices.
        struct Mesh_index;
        
        /// @brief A template that combines `Mesh_index` and `Local_index`
        /// to create a unique index for primitives within the scene.
        /// \tparam Local_index is any of the indices provided through `Surface_mesh`, e.g., `Surface_mesh::Vertex_index`
        template<class Local_index> struct Scene_index;

        /// @}

        /// \name Types
        /// @{

        /// @brief The underlying Kernel used by the contained `Collision_mesh` objects.
        typedef          Kernel                         K;

        /// @brief The type of collision mesh contained in the scene.
        typedef          Collision_mesh<K>              Mesh;

        /// @brief An alias for `Collision_mesh::Point`.
        typedef typename K::Point_3                     Point;

        /// @brief An alias for `Collision_mesh::Vector`.
        typedef typename K::Vector_3                    Vector;

        /// @brief An alias for `Collision_mesh::Vertex_index`.
        typedef typename Mesh::Vertex_index             Vertex_index;

        /// @brief An alias for `Collision_mesh::Edge_index`.
        typedef typename Mesh::Edge_index               Edge_index;

        /// @brief An alias for `Collision_mesh::Halfedge_index`.
        typedef typename Mesh::Halfedge_index           Halfedge_index;

        /// @brief An alias for `Collision_mesh::Face_index`.
        typedef typename Mesh::Face_index               Face_index;

        /// @brief An alias for `Collision_mesh::Face_range`.
        typedef typename Mesh::Face_range               Face_range;

        /// @brief An alias for `Collision_mesh::Vertex_range`.
        typedef typename Mesh::Vertex_range             Vertex_range;

        /// @brief An index type that uniquely identifies a vertex by the combination of 
        /// its local `Surfaced_mesh::Vertex_index` and the corresponding `Mesh_index`.
        typedef          Scene_index<Vertex_index>          Scene_vertex_index;

        /// @brief An index type that uniquely identifies a vertex by the combination of 
        /// its local `Surfaced_mesh::Face_index` and the corresponding `Mesh_index`.
        typedef          Scene_index<Face_index>            Scene_face_index;

        /// @brief A vector type containing `Scene_vertex_index` objects associated
        /// with vertices in the scene.
        typedef          std::vector<Scene_vertex_index>    Scene_vertex_range;

        /// @brief A vector type containing `Scene_face_index` objects associated
        /// with faces in the scene.
        typedef          std::vector<Scene_face_index>      Scene_face_range;

        /// @brief A type that contains six pointers to the `Point_3` objects that comprise a triangle trajectory
        typedef          ::CGAL::Collisions::internal::Triangle_trajectory_observer<K, Scene_face_index>        Trajectory;

        /// @brief A type that wraps `Collision_scene::Trajectory` for use in an `AABB_tree`.
        typedef          ::CGAL::Collisions::internal::AABB_Triangle_trajectory_primitive<K, Scene_face_index>  Trajectory_primitive;

        /// @brief A vector type containing `Collision_scene::Trajectory` objects associated
        /// with the scene.
        typedef          std::vector<Trajectory>                                    Trajectory_range;

        /// @brief The traits type used in `Collision_scene::Tree`
        typedef          ::CGAL::AABB_traits<K, Trajectory_primitive>               AABB_traits;

        /// @brief The type of `AABB_tree`.
        typedef          AABB_tree<AABB_traits>                                     Tree;

        /// @brief An alias for `AABB_tree::Primitive_id`.
        typedef typename Tree::Primitive_id                                         Primitive_id;

        /// @}

    private:
        Trajectory make_trajectory(const Mesh & mesh, const Scene_face_index& scene_face_index);
        Scene_vertex_range  vertices_;
        Scene_face_range    faces_;
        Trajectory_range    trajectories_;
        std::vector<Mesh>&  meshes_;
        Tree                trajectory_tree_;

    public:


        /// \name Creation
        /// @{

        /// @brief Instantiates a scene defined by the vector of collision meshes it contains, each
        /// of which are assigned a `Mesh_index` corresponding to their index in the vector.
        explicit Collision_scene(std::vector<Mesh> & meshes);

        /// @}

        /// \name Methods
        /// @{

        /// @brief Returns a vector of all the vertex indices contained in the scene.
        const Scene_vertex_range& vertices() const;

        /// @brief Returns a vector of all the face indices contained in the scene.
        const Scene_face_range& faces() const;

        /// @brief Returns a reference to the `AABB_tree` that covers all the triangle 
        /// trajectories in the scene.
        const Tree& tree() const;

        /// @brief returns a single collision mesh, the result of joining all collision 
        /// meshes contained in the scene.
        /// @details The join operation boostraps `Surface_mesh::join()`
        Collision_mesh<K> joined_meshes();

        /// @brief Applies the specified `CGAL::Color` to the face corresponding to 
        /// the given `Scene_face_index`.
        void color(const Scene_face_index& ti, CGAL::IO::Color c);

        /// @brief An interface by which a user can update any data associated with the vertex indices of the 
        /// scene using a functor that models the `CollisionSceneUpdateFunctor` concept.
        /// \tparam UpdateFunctor is a model of the concept `CollisionSceneUpdateFunctor`.
        template <class UpdateFunctor>
        void update_state(UpdateFunctor& update_functor, bool update_tree=false);

        /// @brief Recomputes the `AABB_tree` that covers all the triangle trajectories 
        /// contained in the scene.
        void update_tree();

        /// @}
};


// ============
// Constructors
// ============
template <class K>
Collision_scene<K>::Collision_scene(std::vector<Mesh> & meshes) : meshes_{meshes}
{
    size_t i{0};
    for( const auto & mesh_ : meshes_ )
    {
        Mesh_index current_mesh_index{i};

        vertices_.reserve(vertices_.size() + mesh_.num_vertices());
        for( const auto & current_vertex_index : mesh_.vertices()) {
            vertices_.push_back(
                Scene_vertex_index{
                    current_mesh_index,
                    current_vertex_index
                }
            );
        }

        size_t k{faces_.size()};
        faces_.reserve(faces_.size() + mesh_.num_faces());
        trajectories_.reserve(trajectories_.size() + mesh_.num_faces());
        for( const auto & current_face_index : mesh_.faces()) {
            faces_.push_back(
                Scene_face_index{
                    current_mesh_index,
                    current_face_index
                }
            );
            trajectories_.push_back(
                make_trajectory(mesh_, faces_[k])
            );
            ++k;
        }

        ++i;
    }

    trajectory_tree_ = Tree(trajectories_.begin(),trajectories_.end());
};

/// @}


// ========================
// Member Class definitions
// ========================
template <class K>
template <class Local_index>
struct Collision_scene<K>::Scene_index : std::pair<Mesh_index, Local_index> {

    typedef std::pair<Mesh_index, Local_index> Base;
    typedef std::size_t size_type;

    Scene_index() {}

    Scene_index(Mesh_index mi, Local_index li)
    : Base(mi,li)
    {}

    Scene_index(const Base& index_pair)
    : Base(index_pair)
    {}

    friend std::ostream& operator<<(std::ostream& os, Scene_index const& Scene_index)
    {
        return (os << "(" << Scene_index.first << ", " << Scene_index.second << ")" );
    }

    bool operator<(Scene_index const& other_Scene_index) const {
        bool less_than = this->first == other_Scene_index.first
            ? this->second < other_Scene_index.second
            : this->first < other_Scene_index.first;
        return less_than;
    }

    Mesh_index mesh_index() const {
        return this->first;
    }

    Local_index local_index() const {
        return this->second;
    }

};

template <class K>
struct Collision_scene<K>::Mesh_index {

    typedef std::size_t size_type;

    size_type idx_;

    Mesh_index() {}
    explicit Mesh_index(size_type idx) : idx_{ idx } {}

    friend std::ostream& operator<<(std::ostream& os, Mesh_index const& mesh_index)
    {
        return (os << "m" << mesh_index.idx_ );
    }

    bool operator==(Mesh_index const& other_mesh_index) const {
        return this->idx_ == other_mesh_index.idx_;
    }

    bool operator<(Mesh_index const& other_mesh_index) const {
        return this->idx_ < other_mesh_index.idx_;
    }

    operator size_type() const { return idx_; }

    size_type idx() const {
        return idx_;
    }

    size_type id() const {
        return idx_;
    }
};

// ===============
// Member routines
// ===============
template <class K>
const typename Collision_scene<K>::Tree& Collision_scene<K>::tree() const
{
    return this->trajectory_tree_;
}

template <class K>
const typename Collision_scene<K>::Scene_vertex_range& Collision_scene<K>::vertices() const
{
    return this->vertices_;
}

template <class K>
const typename Collision_scene<K>::Scene_face_range& Collision_scene<K>::faces() const
{
    return this->faces_;
}

template <class K>
auto Collision_scene<K>::make_trajectory(
    const Mesh& mesh,
    const Scene_face_index& scene_face_index
) -> Trajectory
{

    size_t j{0};
    std::vector<const Point*> trajectory_points(6);

    for( const auto & vert_index : mesh.vertices_around_face(mesh.halfedge(scene_face_index.local_index())) )
    {
        trajectory_points.at(j) = & mesh.point(vert_index);
        trajectory_points.at(j+1) = & mesh.next_point(vert_index);

        j += 2;
    }

    return Trajectory(
        trajectory_points[0], // v0_current
        trajectory_points[2], // v0_next
        trajectory_points[4], // v1_current
        trajectory_points[1], // v1_next
        trajectory_points[3], // v2_current
        trajectory_points[5], // v2_next
        scene_face_index
    );
};

template <class K>
Collision_mesh<K> Collision_scene<K>::joined_meshes()
{
    Mesh joined_mesh{meshes_[0]};

    auto begin = meshes_.begin();
    auto end = meshes_.end();
    begin++;

    std::for_each(
        begin,
        end,
        [&joined_mesh](const auto& m){
            joined_mesh += m;
        }
    );

    return joined_mesh;
};

template <class K>
void Collision_scene<K>::color(const Scene_face_index& ti, CGAL::IO::Color c) {
    meshes_[ti.mesh_index().idx()].color(ti.local_index(), c);
    return;
};

template <class K>
template <class UpdateFunctor>
void Collision_scene<K>::update_state(UpdateFunctor& update_functor, bool update_tree_)
{
    Mesh* mesh_ptr = nullptr;
    for( auto& scene_vertex_index : vertices_)
    {
        mesh_ptr = &(meshes_[scene_vertex_index.mesh_index()]);
        update_functor(mesh_ptr, scene_vertex_index);
    }

    if( update_tree_ ) {
        update_tree();
    }
}

template <class K>
void Collision_scene<K>::update_tree()
{
    for( auto& t : trajectories_ )
    {
        t.update();
    }
    // TODO: use the update functor that sloriot provided
    trajectory_tree_ = Tree(trajectories_.begin(),trajectories_.end());
}


}

#endif