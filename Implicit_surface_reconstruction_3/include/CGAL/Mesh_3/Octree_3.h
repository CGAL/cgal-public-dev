// Copyright (c) 2007-2008  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Tong Zhao

#ifndef CGAL_MESH_3_OCTREE_3_H
#define CGAL_MESH_3_OCTREE_3_H

#include <CGAL/license/Implicit_surface_reconstruction_3.h>
#include <CGAL/bounding_box.h>

#include <stack>
#include <vector>

#endif // CGAL_MESH_3_OCTREE_3_H

namespace CGAL {

// F B U D L R
const static bool NU[8][6] = {false, true,  false, true,  false, true,
                              false, true,  false, true,  true,  false,
                              false, true,  true,  false, false, true,
                              false, true,  true,  false, true,  false,
                              true,  false, false, true,  false, true,
                              true,  false, false, true,  true,  false,
                              true,  false, true,  false, false, true,
                              true,  false, true,  false, true,  false};

enum direction{ FRONT  = 0,
                BACK   = 1,
                UP     = 2,
                DOWN   = 3,
                LEFT   = 4,
                RIGHT  = 5};

template < class Kernel,
           class PointRange >
class OctreeNode
{   
// Public types
public:

    typedef typename Kernel::FT                 FT;
    typedef typename Kernel::Point_3            Point;
    typedef typename Kernel::Vector_3           Vector;

    typedef typename PointRange::const_iterator InputIterator;
    typedef typename std::vector<InputIterator> IterList;

    typedef typename OctreeNode<Kernel>         Node;

    OctreeNode(): 
        m_isleaf(true),
        m_children(NULL),
        m_parent(NULL),
        m_level(0),
        m_length(1),
        m_size(0),
        m_index(0),
        m_steiner(false)
    {}

    ~OctreeNode(){
        unpslit();
    }

    void unsplit()
    {
        if(m_children != NULL)
        {
            for(int i = 0; i < 8; i++)
                m_children[i].unsplit();
            delete [] m_children;
            m_children = NULL;
        }

        m_isleaf = true;
    }

    void set_parent(Node *parent)
    {
        m_parent = parent;
        m_level = parent->level() + 1;
        m_length = parent->length() / 2.;
    }

    void split()
    {
        m_children = new Node[8];
        const FT offset = m_length / 4.;

        for(int i=0;i<8;i++){
            m_children[i].set_parent(this);
            m_children[i].index() = i;
        }

        m_children[0].point() = Point(m_point.x() - offset, m_point.y() + offset, m_point.z() + offset);
        m_children[1].point() = Point(m_point.x() + offset, m_point.y() + offset, m_point.z() + offset);
        m_children[2].point() = Point(m_point.x() - offset, m_point.y() + offset, m_point.z() - offset);
        m_children[3].point() = Point(m_point.x() + offset, m_point.y() + offset, m_point.z() - offset);
        m_children[4].point() = Point(m_point.x() - offset, m_point.y() - offset, m_point.z() + offset);
        m_children[5].point() = Point(m_point.x() + offset, m_point.y() - offset, m_point.z() + offset);
        m_children[6].point() = Point(m_point.x() - offset, m_point.y() - offset, m_point.z() - offset);
        m_children[7].point() = Point(m_point.x() + offset, m_point.y() - offset, m_point.z() - offset);

        m_isleaf = false;
    } 

    Node *parent() { return m_parent; }

    unsigned int& level() { return m_level; }
    const unsigned int& level() const { return m_level; }

    unsigned char& index() { return m_index; }
    const unsigned char& index() const { return m_index; }

    const Point& point() const { return m_point; }
    Point& point() { return m_point; }

    const FT& length() const { return m_length; }
    FT& length() { return m_length; }

    bool& is_leaf() { return m_isleaf; }
    const bool& is_leaf() const { return m_isleaf; }

    Node *child(const unsigned int index) {
        if(m_isleaf || index > 7)
            return NULL;
        else
            return m_children[index];
    }


private:
    size_t              m_size;
    Point               m_point;
    IterList            m_pts;
    Node*               m_children;
    bool                m_isleaf;
    bool                m_steiner;
    unsigned int        m_level;
    unsigned char       m_index;
    Node*               m_parent;
    FT                  m_length;

}; // end class OctreeNode


template < class Gt,
           class PointRange,
           class PointMap,
           class NormalMap >
class Octree
{   
// Public types
    /// Geometric traits class
    typedef Gt  Geom_traits;
    typedef typename Geom_traits::FT            FT;
    typedef typename Geom_traits::Iso_cuboid_3  Iso_cuboid;
    typedef typename OctreeNode<Gt>             Node;

public:

    Octree(
        PointRange& points, ///< input point range
        PointMap point_map, ///< property map: `value_type of InputIterator` -> `Point` (the position of an input point).
        NormalMap normal_map, ///< property map: `value_type of InputIterator` -> `Vector` (Normal)
        const FT enlarge_ratio = 1.2):
        m_ranges(points),
        m_points(point_map),
        m_normals(normal_map)
    {
        m_bounding_box = CGAL::bounding_box(points.begin(), points.end());
        m_center = midpoint((ic.min)(), (ic.max)());
        m_box_length = 2 * enlarge_ratio * squared_distance(center, (ic.max)());

        m_root.point() = m_center;
        m_root.length() = m_box_length;
        m_root.size() = points.size();
    }

    ~Octree(){
        my_root.unsplit();
    }

    template < typename CellCriteria,
               typename NormalCriteria > // or other useful criterion 
    void build(size_t max_depth, size_t max_pts_num, CellCriteria cell_criteria, NormalCriteria normal_criteria){
        
    }

    void grade(){

    }

    template <typename VertexOutputIterator>>
    void generate_balanced_pts(VertexOutputIterator vertices_out){  // std::back_inserter(new_vertices)

    }

private:

    PointMap    m_points;
    PointRange  m_ranges;
    NormalMap   m_normals;
    Iso_cuboid  m_bounding_box;
    FT          m_box_length;
    Point       m_center;
    Node        m_root;
    unsigned int m_depth;



}; // end class Octree




} // namespace CGAL