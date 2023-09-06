
namespace CGAL {

/*!
\ingroup PkgCollisions3Classes

The class `Collision_scene` serves as an aggregator of `Collision_mesh` objects
to be evaluated for mutual collisions. The class provides standardized indexing
across the primitives of all contained collision meshes, as well as an AABB Tree
that tracks the trajectories of these primitives.

`Collision_scene` serves as the standard data structure for querying and 
visualizing collisions between an arbitrary number of collision meshes.
*/
template <class Kernel>
class Collision_scene {

    public:

        /// \name Member Classes
        /// @{
        
        /*!
        a template that combines `Mesh_index` and `Local_index`
        to create a unique index for primitives within the scene.

        \tparam Local_index is any of the indices provided through `Surface_mesh`, e.g., `Surface_mesh::Vertex_index`
        */
        template<class Local_index> struct Scene_index;
        
        /*!
        a simple index wrapper provided to enable pretty-printing of mesh indices.
        */
        struct Mesh_index;

        ///@}


        /// \name Types
        /// @{
        
        /*!
        the underlying Kernel used by the contained `Collision_mesh` objects.
        */
        typedef          Kernel                         K;
        
        /*!
        the type of collision mesh contained in the scene.
        */
        typedef          Collision_mesh<K>              Mesh;
        
        /*!
        an alias for `Collision_mesh::Point`.
        */
        typedef typename Mesh::Point                    Point;
        
        /*!
        an alias for `Collision_mesh::Vector`.
        */
        typedef typename Mesh::Vector                   Vector;
        
        /*!
        an alias for `Collision_mesh::Vertex_index`.
        */
        typedef typename Mesh::Vertex_index             Vertex_index;
        
        /*!
        an alias for `Collision_mesh::Edge_index`.
        */
        typedef typename Mesh::Edge_index               Edge_index;
        
        /*!
        an alias for `Collision_mesh::Halfedge_index`.
        */
        typedef typename Mesh::Halfedge_index           Halfedge_index;
        
        /*!
        an alias for `Collision_mesh::Face_index`.
        */
        typedef typename Mesh::Face_index               Face_index;
        
        /*!
        an alias for `Collision_mesh::Face_range`.
        */
        typedef typename Mesh::Face_range               Face_range;
        
        /*!
        an alias for `Collision_mesh::Vertex_range`.
        */
        typedef typename Mesh::Vertex_range             Vertex_range;
        
        /*!
        an index that uniquely identifies a vertex by the combination of its local `Surfaced_mesh::Vertex_index` and the corresponding `Mesh_index`.
        */
        typedef          Scene_index<Vertex_index>                                  Scene_vertex_index;
        
        /*!
        an index that uniquely identifies a vertex by the combination of its local `Surfaced_mesh::Face_index` and the corresponding `Mesh_index`.
        */
        typedef          Scene_index<Face_index>                                    Scene_face_index;
        
        /*!
        an iterator over all contained `Scene_vertex_index` objects.
        */
        typedef          std::vector<Scene_vertex_index>                            Scene_vertex_range;
        
        /*!
        an iterator over all contained `Scene_face_index` objects.
        */
        typedef          std::vector<Scene_face_index>                              Scene_face_range;
        
        /*!
        a container of pointers to the `Point_3` objects that comprise a triangle trajectory
        */
        typedef          ::CGAL::Collisions::internal::Triangle_trajectory_observer<K, Scene_face_index>        Trajectory;
        
        /*!
        a primitive that wraps `Collision_scene::Trajectory` for use in an `AABB_tree`.
        */
        typedef          ::CGAL::Collisions::internal::AABB_Triangle_trajectory_primitive<K, Scene_face_index>  Trajectory_primitive;
        
        /*!
        an iterator over all contained `Collision_scene::Trajectory` objects.
        */
        typedef          std::vector<Trajectory>                                    Trajectory_range;
        
        /*!
        an alias for `CGAL::AABB_traits`.
        */
        typedef          ::CGAL::AABB_traits<K, Trajectory_primitive>               AABB_traits;
        
        /*!
        an alias for `CGAL::AABB_tree<AABB_traits>`.
        */
        typedef          ::CGAL::AABB_tree<AABB_traits>                                     Tree;
        
        /*!
        an alias for `CGAL::AABB_tree::Primitive_id`.
        */
        typedef typename Tree::Primitive_id                                         Primitive_id;

        ///@}


        /// \name Creation
        /// @{

        /*!
        instantiates a scene defined by its containment of the specified vector of meshes, which
        are assigned a `Mesh_index` corresponding to their index in the vector.
        */
        explicit Collision_scene(std::vector<Mesh> & meshes);

        /// @}


        /// \name Methods
        /// @{
        
        /*!
        returns a range over all the vertex indices contained in the scene.
        */
        const Scene_vertex_range& vertices() const;
        
        /*!
        returns a range over all the face indices contained in the scene.
        */
        const Scene_face_range& faces() const;
        
        /*!
        returns a reference to the `AABB_tree` that covers all the triangle trajectories in the scene.
        */
        const Tree& tree() const;

        /*!
        returns a single collision mesh, the result of joining all collision meshes contained in the scene.
        */
        Collision_mesh<K> joined_meshes();

        /*!
        applies the specified `CGAL::Color` to the face corresponding to the given `Scene_face_index`.
        */
        void color(const Scene_face_index& ti, CGAL::IO::Color c);

        /*!
        provides an interface by which a user can update any data associated with the vertex indices of the 
        scene using a functor that models the `UpdateFunctor` concept.

        \tparam UpdateFunctor is a model of the concept `UpdateFunctor`.
        */
        template <class UpdateFunctor>
        void update_state(UpdateFunctor& update_functor, bool update_tree=false);

        /*!
        recomputes the `AABB_tree` that covers all the triangle trajectories contained in the scene.
        */
        void update_tree();

        ///@}
};

}