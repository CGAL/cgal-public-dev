
namespace CGAL {

/*!
\ingroup PkgCollisions3Classes

The class `Collision_mesh` derives from `Surface_mesh`. It's primary purpose 
is to formalize the requirement that both the current and next positions of 
each point are known. This requirement enables the use of collision detection
algorithms.
*/
template <typename Kernel>
class Collision_mesh : public Surface_mesh<typename Kernel::Point_3> {

    public:
        /// \name Types
        /// @{

        /*!
        the underlying Kernel class
        */
        typedef          Kernel                             K;

        /*!
        the type of Point_3
        */
        typedef typename K::Point_3                         Point;

        /*!
        the type of Vector_3
        */
        typedef typename K::Vector_3                        Vector;

        /*!
        the type of `Surface_mesh` from which this class derives
        */
        typedef          Surface_mesh<Point>                Base;

        /*!
        alias for `Surface_mesh::Vertex_index`
        */
        typedef typename Base::Vertex_index                 Vertex_index;

        /*!
        alias for `Surface_mesh::Edge_index`
        */
        typedef typename Base::Edge_index                   Edge_index;

        /*!
        alias for `Surface_mesh::Halfedge_index`
        */
        typedef typename Base::Halfedge_index               Halfedge_index;

        /*!
        alias for `Surface_mesh::Face_index`
        */
        typedef typename Base::Face_index                   Face_index;

        /*!
        alias for `Surface_mesh::size_type`
        */
        typedef typename Base::size_type                    size_type;

        /*!
        the type of `Point_3_trajectory`, which contains current and next positions of a `Point_3` object
        */
        typedef          ::CGAL::Point_3_trajectory<K>      Point_trajectory;
        
        /*!
        the type of `Segment_3_trajectory`, which contains current and next positions of a `Segment_3` object
        */
        typedef          ::CGAL::Segment_3_trajectory<K>    Segment_trajectory;
        
        /*!
        the type of `Triangle_3_trajectory`, which contains current and next positions of a `Triangle_3` object
        */
        typedef          ::CGAL::Triangle_3_trajectory<K>   Triangle_trajectory;
        
        /*!
        property map correlating a `Vertex_index` to a `Point_3` object
        */
        typedef unspecified_type Point_map;
        
        /*!
        property map correlating a `Vertex_index` to a color
        */
        typedef unspecified_type Vertex_color_map;
        
        /*!
        property map correlating a `Edge_index` to a color
        */
        typedef unspecified_type Edge_color_map;
        
        /*!
        property map correlating a `Face_index` to a color
        */
        typedef unspecified_type Face_color_map;

        /// @}
    
        // ===================================================================

        /// \name Creation
        /// @{

        /*!
        creates an empty collision mesh
        */
        Collision_mesh();

        /*!
        creates a collision mesh from a surface mesh, defaulting the next
        position of each point to its current position.
        */
        Collision_mesh(const Base& surface_mesh);

        /*!
        move constructor that creates a collision mesh from a surface mesh,
        defaulting the next position of each point to its current position.
        */
        Collision_mesh(Base&& surface_mesh);

        /*!
        copy constructor
        */
        Collision_mesh(const Collision_mesh& collision_mesh);    

        /*!
        move constructor
        */
        Collision_mesh(Collision_mesh&& collision_mesh);
    
        /// @}
    
        // ===================================================================

        /// \name Operators
        /// @{

        /*!
        copy assignment
        */
        Collision_mesh<K>& operator=(const Collision_mesh<K>& rhs); // assigns `rhs` to `*this`. Performs a deep copy of all properties.       

        /*!
        move assignment
        */
        Collision_mesh<K>& operator=(Collision_mesh<K>&& collision_mesh); // move assignment

        /// @}
    
        // ===================================================================

        /// \name Methods
        /// @{

        /*!
        returns the velocity vector computed as the difference between
        the current and nex position of the vertex.
        */
        const Vector& velocity(Vertex_index v) const;

        /*!
        returns the velocity vector computed as the difference between
        the current and nex position of the vertex.
        */
        Vector& velocity(Vertex_index v);

        /*!
        returns the next position associated with the vertex index.
        */
        const Point& next_point(Vertex_index v) const;

        /*!
        returns the next position associated with the vertex index.
        */
        Point& next_point(Vertex_index v);

        /*!
        returns the `Point_trajectory_3` associated with the vertex index.
        */
        Point_trajectory point_trajectory(Vertex_index v) const;

        /*!
        returns the `Segment_trajectory_3` associated with the halfedge index.
        */
        Segment_trajectory edge_trajectory(Halfedge_index h) const;

        /*!
        returns the `Triangle_trajectory_3` associated with the face index.
        */
        Triangle_trajectory face_trajectory(Face_index f) const;

        /*!
        assigns the specified color to the face associated with the face index.
        */
        void color(const Face_index& fi, CGAL::IO::Color c);

        /// @}
};

}