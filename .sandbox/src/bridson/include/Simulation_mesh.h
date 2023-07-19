#ifndef SIMULATION_MESH_H
#define SIMULATION_MESH_H

#include <CGAL/Surface_mesh.h>

namespace CGAL {

template <class K_>
class Simulation_mesh {

    public:
        typedef typename K_                                 K;
        typedef typename K::Point_3                         Point;
        typedef typename K::Vector_3                        Vector;
        typedef typename Surface_mesh<Point>                Base;
        typedef typename Base::Vertex_index                 Vertex_index;
        typedef typename Base::Edge_index                   Edge_index;
        typedef typename Base::Halfedge_index               Halfedge_index;
        typedef typename Base::Face_index                   Face_index;
        typedef typename Base::Face_range                   Face_range;
        typedef typename Base::Vertex_around_target_range   Vertex_around_target_range;

    private:  
        Base & mesh_;
        decltype(mesh_.add_property_map<Vertex_index, Vector>("v:velocity").first) vvelocity_;
        decltype(mesh_.add_property_map<Vertex_index, Point>("v:next_point").first) vnext_point_;

    public:
        Simulation_mesh() {}
        Simulation_mesh(Base & mesh) : mesh_{mesh} {
            vvelocity_ = mesh.add_property_map<Vertex_index, Vector>("v:velocity").first;
            vnext_point_ = mesh.add_property_map<Vertex_index, Point>("v:next_point").first;
            for(Vertex_index vd : mesh.vertices()){
                put(vvelocity_, vd, Vector(0, 0, 0));
                put(vnext_point_, vd, point(vd));
            }
        }
        Simulation_mesh(Base && mesh) : mesh_{mesh} {
            vvelocity_ = mesh.add_property_map<Vertex_index, Vector>("v:velocity").first;
            vnext_point_ = mesh.add_property_map<Vertex_index, Point>("v:next_point").first;
            for(Vertex_index vd : mesh.vertices()){
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

        Face_range faces();
        Vertex_around_target_range vertices_around_face(Halfedge_index h) const;
        Halfedge_index halfedge(Face_index f) const;

};

    /// returns the velocity associated to vertex `v`.
    template <class K>
    const typename K::Vector_3& Simulation_mesh<K>::velocity(Vertex_index v) const { return vvelocity_[v]; }

    template <class K>
    typename K::Vector_3& Simulation_mesh<K>::velocity(Vertex_index v) { return vvelocity_[v]; }

    template <class K>
    const typename K::Point_3& Simulation_mesh<K>::point(Vertex_index v) const { return mesh_.point(v); }
    
    template <class K>
    typename K::Point_3& Simulation_mesh<K>::point(Vertex_index v) { return mesh_.point(v); }

    template <class K>
    const typename K::Point_3& Simulation_mesh<K>::next_point(Vertex_index v) const { return vnext_point_[v]; }
    
    template <class K>
    typename K::Point_3& Simulation_mesh<K>::next_point(Vertex_index v) { return vnext_point_[v]; }

    template <class K>
    typename Simulation_mesh<K>::Face_range Simulation_mesh<K>::faces() { return mesh_.faces(); } 

    template <class K>
    typename Simulation_mesh<K>::Vertex_around_target_range Simulation_mesh<K>::vertices_around_face(Halfedge_index h) const
    {
      return mesh_.vertices_around_target(h);
    }

    template <class K>
    typename Simulation_mesh<K>::Halfedge_index Simulation_mesh<K>::halfedge(Face_index f) const
    {
        return mesh_.halfedge(f);
    }

}

#endif