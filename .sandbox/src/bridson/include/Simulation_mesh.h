#include <CGAL/Surface_mesh.h>

namespace CGAL {

template <class K>
class Simulation_mesh : public Surface_mesh<typename K::Point_3>
{
    public:
        typedef typename K::Point_3     Point;
        typedef typename K::Vector_3    Vector;
        typedef Surface_mesh<Point>     Base;
        typedef Base::Vertex_index      vertex_descriptor;
        typedef Base::Edge_index        edge_descriptor;
        typedef Base::Halfedge_index    halfedge_descriptor;
        typedef Base::Face_index        face_descriptor;

    private:  
        Property_map<vertex_descriptor, Vector> vvelocity_;

    public:
        Simulation_mesh() : Base() {}
        Simulation_mesh(const Base & mesh) : Base(mesh) {
            vvelocity_ = add_property_map<vertex_descriptor, Vector>("v:velocity").first;
            for(vertex_descriptor vd : vertices()){
                put(vvelocity_, vd, Vector(0, 0, 0));
            }
        }
        Simulation_mesh(Surface_mesh&& mesh) : Base(mesh) {
            vvelocity_ = add_property_map<vertex_descriptor, Vector>("v:velocity").first;
            for(vertex_descriptor vd : vertices()){
                put(vvelocity_, vd, Vector(0, 0, 0));
            }
        }

        const Vector& velocity(vertex_descriptor v) const;
        Vector& velocity(vertex_descriptor v);
        
};

    /// returns the velocity associated to vertex `v`.
    template <class K>
    const typename K::Vector_3& Simulation_mesh<K>::velocity(vertex_descriptor v) const { return vvelocity_[v]; }

    /// returns the velocity associated to vertex `v`.
    template <class K>
    typename K::Vector_3& Simulation_mesh<K>::velocity(vertex_descriptor v) { return vvelocity_[v]; }

}