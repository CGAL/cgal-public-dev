#ifndef CGAL_OCTREE_GRID_DOMAIN_H
#define CGAL_OCTREE_GRID_DOMAIN_H

#include "Octree_wrapper.h"

#include <array>

namespace CGAL {

    template<typename GeomTraits>
    class Octree_domain {
      public:
        typedef GeomTraits Geom_traits;
        typedef typename Geom_traits::FT FT;
        typedef typename Geom_traits::Point_3 Point_3;
        typedef typename Geom_traits::Vector_3 Vector_3;

        // octree type and handles
        typedef Octree_wrapper<Geom_traits> Octree;
        typedef typename Octree::Vertex_handle Vertex_handle;
        typedef typename Octree::Edge_handle Edge_handle;
        typedef typename Octree::Voxel_handle Voxel_handle;

        // return types of topology methods
        typedef std::array<Vertex_handle, 2> Edge_vertices;
        typedef std::array<Voxel_handle, 4> Voxels_incident_to_edge;
        typedef std::array<Vertex_handle, 8> Voxel_vertices;
        typedef std::array<Edge_handle, 12> Voxel_edges;

      public:
        Octree_domain( const Octree& octree ) : octree_( &octree ) {}

        // returns the position of vertex v in 3D space
        Point_3 position( const Vertex_handle& v ) const { return octree_->point( v ); }

        // not necessary! this should be removed
        Vector_3 gradient( const Vertex_handle& v ) const { return octree_->gradient( v ); }

        // returns the stored value of vertex v
        FT value( const Vertex_handle& v ) const { return octree_->vertex_value( v ); }

        // topology methods:

        // returns the two vertices incident to edge e
        Edge_vertices edge_vertices( const Edge_handle& e ) const { return octree_->edge_vertices( e ); }

        // returns all voxels incident to edge e (usually the 4 neighbor cells)
        Voxels_incident_to_edge voxels_incident_to_edge( const Edge_handle& e ) const { return octree_->edge_voxels( e ); }

        // returns all vertices of the cell vox (usually 8 for a cube)
        Voxel_vertices voxel_vertices( const Voxel_handle& vox ) const { return octree_->voxel_vertices( vox ); }

        // currently not used! returns all edges of the cell vox (usually 12 for a cube)
        Voxel_edges voxel_edges( const Voxel_handle& vox ) const { return octree_->voxel_edges( vox ); }

        // iterators: for examples see Dual_contouring_octree_3.h:DC_position_functor and DC_quad_functor

        // iterate over all (leaf) vertices of the octree and call the functor f on each one
        template<typename Functor>
        void iterate_vertices( Functor& f ) const {
            for( const Vertex_handle& v: octree_->leaf_vertices() ) {
                f( v );
            }
        }

        // iterate over all (leaf) edges of the octree and call the functor f on each one
        template<typename Functor>
        void iterate_edges( Functor& f ) const {
            for( const Edge_handle& e: octree_->leaf_edges() ) {
                f( e );
            }
        }

        // iterate over all (leaf) cells of the octree and call the functor f on each one
        template<typename Functor>
        void iterate_voxels( Functor& f ) const {
            for( const Voxel_handle& v: octree_->leaf_voxels() ) {
                f( v );
            }
        }

        const Octree& getOctree() const { return *octree_; }

      private:
        // reference to the Octree_wrapper
        const Octree* octree_;
    };

}    // namespace CGAL

#endif    // CGAL_OCTREE_GRID_DOMAIN_H
