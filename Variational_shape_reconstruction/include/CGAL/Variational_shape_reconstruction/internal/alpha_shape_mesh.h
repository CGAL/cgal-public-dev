// Copyright (c) 2018  Liangliang Nan. All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Liangliang Nan

#ifndef ALPHA_SHAPE_MESH_H
#define ALPHA_SHAPE_MESH_H

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Surface_mesh.h>

#include "types.h"

namespace qem {

/* A vertex class with an additional member representing its index */
template < class Gt, class VB = CGAL::Triangulation_hierarchy_vertex_base_2<Gt> >
class AS_vertex_base : public VB
{
public:
    typedef VB                                  Base;
    typedef typename Base::Vertex_handle        Vertex_handle;
    typedef typename Base::Face_handle          Face_handle;
    typedef typename Base::Point                Point;

    template < typename TDS2 >
    struct Rebind_TDS {
        typedef typename VB::template Rebind_TDS<TDS2>::Other	VB2;
        typedef AS_vertex_base<Gt, VB2>                         Other;
    };

public:
    AS_vertex_base() : Base(), index_(-1) {}
    AS_vertex_base(const Point & p) : Base(p), index_(-1) {}
    AS_vertex_base(const Point & p, Face_handle f) : Base(f, p), index_(-1) {}
    AS_vertex_base(Face_handle f) : Base(f), index_(-1) {}

    void set_index(int idx) { index_ = idx; }
    int  index() const { return index_; }

private:
    int index_;
}; // end of class AS_vertex_base


template <typename Ht>
class Alpha_shape : public CGAL::Alpha_shape_2<Ht>
{
public:
    typedef CGAL::Alpha_shape_2<Ht>					Parent_class;
    typedef typename Ht::Point_2					Point2;
    typedef typename Parent_class::Vertex_handle	Vertex_handle;

public:
    // constructs alpha shapes from the input points
    template <typename InputIterator>
    Alpha_shape(InputIterator first, InputIterator beyond)
    {
        InputIterator it = first;
        for (int id = 0; it != beyond; ++it, ++id) {
            const Point2& p = *it;
            Vertex_handle vh = Ht::insert(p);
            if (vh->index() == -1)
                vh->set_index(id);
        }

        if (Parent_class::dimension() == 2) {
            // Computes the associated _interval_face_map
            Parent_class::initialize_interval_face_map();

            // Computes the associated _interval_edge_map
            Parent_class::initialize_interval_edge_map();

            // Computes the associated _interval_vertex_map
            Parent_class::initialize_interval_vertex_map();

            // Merges the two maps
            Parent_class::initialize_alpha_spectrum();
        }
    }
}; // end of class Alpha_shape


/**
*	An Alpha Shape Mesh approximates the point covered region by a mesh representation.
*/
class Alpha_shape_mesh
{
    typedef typename Kernel::FT						FT;
    typedef typename Kernel::Point_2				Point2;
    typedef typename Kernel::Point_3				Point3;
    typedef typename Kernel::Plane_3				Plane3;
    typedef CGAL::Surface_mesh<Point3>				Mesh3;
    typedef typename Mesh3::Vertex_index			Vertex_descriptor;

    typedef CGAL::Alpha_shape_vertex_base_2<Kernel>			Avb;
    typedef AS_vertex_base<Avb>								Av;
    typedef CGAL::Triangulation_face_base_2<Kernel>			Tf;
    typedef CGAL::Alpha_shape_face_base_2<Kernel, Tf>		Af;
    typedef CGAL::Triangulation_default_data_structure_2<Kernel, Av, Af> Tds;
    typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>		Dt;
    typedef CGAL::Triangulation_hierarchy_2<Dt>				Ht;

private:
    Alpha_shape<Ht>*			alpha_shape_;
    std::vector<const Point3*>  original_points_;
    Point                       plane_center;
    Vector                      plane_base_1;
    Vector                      plane_base_2;

public:
    /// Given a set of 3D points lying on 'plane', constructs alpha shapes from the
    /// the projection of the points onto 'plane'
    template <typename InputIterator>
    Alpha_shape_mesh(InputIterator first, InputIterator beyond, const Point& center, const Vector& base_1, const Vector& base_2)
    {
        original_points_.clear();

        plane_center = center;
        plane_base_1 = base_1;
        plane_base_2 = base_2;

        std::vector<Point2> pts;
        for (InputIterator it = first; it != beyond; ++it) {
            const Point3& p = *it;
            const Point2& q = project_point_on_plane(p);
            pts.push_back(q);
            original_points_.push_back(&p);
        }

        double dist_2d = std::sqrt(CGAL::squared_distance(pts[0], pts[1]));
        double dist_3d = std::sqrt(CGAL::squared_distance(*original_points_[0], *original_points_[1]));
        //std::cout << "dist 2d: " << dist_2d << ", dist 3d: " << dist_3d << std::endl;
        //std::cout << "check 2d: " << std::sqrt(CGAL::squared_distance(pts[0], pts[2])) << " check 3d: " << std::sqrt(CGAL::squared_distance(*original_points_[0], *original_points_[2])) << std::endl;

        alpha_shape_ = new Alpha_shape<Ht>(pts.begin(), pts.end());
    }

    Point_2 project_point_on_plane(const Point& query)
    {
        Vector qc = query - plane_center;
        double x = plane_base_1 * qc;
        double y = plane_base_2 * qc;

        Point_2 proj(x, y);

        return proj;
    }

    ~Alpha_shape_mesh() 
    { 
        delete alpha_shape_; 
    }

    /// Extracts the 3D mesh representation of the alpha shapes
    bool extract_mesh(FT alpha_value, Mesh3& mesh)
    {
        alpha_shape_->set_alpha(alpha_value);

        typedef std::vector<std::size_t> Triangle;
        std::vector<Triangle>	faces;

        typedef Alpha_shape<Ht>	Alpha_shape;

        typename Alpha_shape::Finite_faces_iterator fit = alpha_shape_->finite_faces_begin();
        for (; fit != alpha_shape_->finite_faces_end(); ++fit) {
            if (alpha_shape_->classify(fit) == Alpha_shape::INTERIOR) {
                Triangle tri;
                for (int i = 0; i < 3; ++i) {
                    typename Alpha_shape::Vertex_handle vh = fit->vertex(i);
                    int idx = vh->index();
                    tri.push_back(idx);
                }
                faces.push_back(tri);
            }
        }

        if (faces.empty())
            return false;

        mesh.clear();

        std::vector<Vertex_descriptor> descriptors(original_points_.size());
        for (std::size_t i = 0; i < original_points_.size(); ++i) {
            const Point3* p = original_points_[i];
            descriptors[i] = mesh.add_vertex(*p);
        }

        for (std::size_t i = 0; i < faces.size(); ++i) {
            std::vector<Vertex_descriptor> face;
            const Triangle& tri = faces[i];
            for (std::size_t j = 0; j < tri.size(); ++j) {
                std::size_t idx = tri[j];
                face.push_back(descriptors[idx]);
            }
            mesh.add_face(face);;
        }

        return true;
    }

}; // end of class Alpha_shape_mesh

} //namespace qem

#endif