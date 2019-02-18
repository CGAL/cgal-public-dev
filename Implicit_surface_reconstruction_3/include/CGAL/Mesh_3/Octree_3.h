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
// Author(s)     : Tong Zhao, Cédric Portaneri

#ifndef CGAL_OCTREE_3_H
#define CGAL_OCTREE_3_H

#include <CGAL/license/Implicit_surface_reconstruction_3.h>
#include <CGAL/bounding_box.h>
#include <boost/iterator/transform_iterator.hpp>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>

#include <stack>
#include <vector>

#endif // CGAL_OCTREE_3_H

namespace CGAL {
	namespace OCTREE {

// F B U D L R
const static bool NU[8][6] = {false, true,  false, true,  false, true,
                              false, true,  false, true,  true,  false,
                              false, true,  true,  false, false, true,
                              false, true,  true,  false, true,  false,
                              true,  false, false, true,  false, true,
                              true,  false, false, true,  true,  false,
                              true,  false, true,  false, false, true,
                              true,  false, true,  false, true,  false};

enum direction { FRONT = 0,
                 BACK  = 1,
                 UP    = 2,
                 DOWN  = 3,
                 LEFT  = 4,
                 RIGHT = 5};

template <class Kernel,
           class PointRange>
class Octree_node
{   
// Public types
public:

    typedef typename Kernel::FT                 FT;
    typedef typename Kernel::Point_3            Point;
    typedef typename Kernel::Vector_3           Vector;

    typedef typename PointRange::const_iterator InputIterator;
    typedef typename std::list<InputIterator> IterList;

    typedef Octree_node<Kernel,PointRange> Node;

    Octree_node(): 
        m_isleaf(true),
        m_children(NULL),
        m_parent(NULL),
        m_level(0),
        m_length(1),
        m_size(0),
        m_index(0)
        //m_steiner(false)
    {}

    ~Octree_node()
	{
		unsplit();
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
        m_length = (FT)0.5 * parent->length();
    }

    void split()
    {
        m_children = new Node[8];
        const FT offset = m_length / 4.; // change

        for(int i = 0; i < 8; i++)
		{
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
    
    const size_t& size() const { return m_size; }
    size_t& size() { return m_size; }

    Node *child(const unsigned int index) 
	{
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
    bool                m_isleaf; // remove
    unsigned int        m_level;
    unsigned char       m_index; // index of current node in 8-subdivision
    Node*               m_parent;
    FT                  m_length;

}; // end class Octree_node


template < class Kernel,
           class PointRange,
           class PointMap,
           class NormalMap >
class Octree
{   
// Public types
    /// Geometric traits class
    typedef typename Kernel::FT            FT;
	typedef typename Kernel::Iso_cuboid_3  Iso_cuboid;
	typedef typename Kernel::Point_3       Point;
	typedef typename Kernel::Vector_3      Vector;
	typedef Octree_node<Kernel, PointRange> Node;

public:

    Octree(
        PointRange& pwn, ///< input point range
        PointMap point_map, ///< property map: `value_type of InputIterator` -> `Point` (the position of an input point).
        NormalMap normal_map, ///< property map: `value_type of InputIterator` -> `Vector` (Normal)
        const FT enlarge_ratio = 1.2):
        m_ranges(pwn),
        m_points(point_map),
        m_normals(normal_map)
    {
        // compute bbox
        typedef typename PointRange::value_type PointRange_t;
        boost::function<Point(PointRange_t&)> pwn_to_point = boost::bind(&PointRange_t::first, _1);
        Iso_cuboid bbox = CGAL::bounding_box(boost::make_transform_iterator(pwn.begin(), pwn_to_point), 
                                             boost::make_transform_iterator(pwn.end(), pwn_to_point));
        debug_bbox(bbox.min(), bbox.max(), "bbox");
        
        // scale bbox
        Iso_cuboid bbox_scaled = bbox.transform(Aff_transformation_3<Kernel>(SCALING, enlarge_ratio));
        Point bbox_centroid = midpoint(bbox.min(), bbox.max());
        Point bbox_scaled_centroid = midpoint(bbox_scaled.min(), bbox_scaled.max());
        Vector diff_centroid = bbox_centroid - bbox_scaled_centroid;
        bbox_scaled = bbox_scaled.transform(Aff_transformation_3<Kernel>(TRANSLATION, diff_centroid));
        debug_bbox(bbox_scaled.min(), bbox_scaled.max(), "bbox_scaled");
       
        // save octree attributes
        m_center = bbox_centroid;
        m_box_length = (FT) sqrt(2.0 * CGAL::squared_distance(bbox_centroid, bbox_scaled.max()));
        m_bounding_box = bbox_scaled.bbox();
        m_root.point() = m_center;
        m_root.length() = m_box_length;
        m_root.size() = pwn.size(); 
    }

    ~Octree()
	{
        m_root.unsplit();
    }

    template < typename CellCriteria,
               typename NormalCriteria > // or other useful criterion 
    void refine(size_t max_depth, 
		size_t max_pts_num, 
		CellCriteria cell_criteria, 
		NormalCriteria normal_criteria)
	{
        
    }

    void grade()
	{
    }

    template <typename OutputIterator> 
    void generate_points(OutputIterator out)
	{  
    }
    
    // (debugging) export bbox defined by min and max as .obj file 
    void debug_bbox(const Point &min, const Point &max, const std::string &filename) {
      std::ofstream out_file(filename+".ply");
      out_file << "ply\n" // header
                  "format ascii 1.0\n"
                  "element vertex 8\n"
                  "property float x\n"
                  "property float y\n"
                  "property float z\n"
                  "element edge 12\n"
                  "property int vertex1\n"            
                  "property int vertex2\n"
                  "end_header\n";
      for(int i = 0; i < 8; i++) { // vertices
        for(int j = 0; j < 3; j++) {
          if((i >> j) & 1)
            out_file << min[j] << " ";
          else
            out_file << max[j] << " ";
        }
        out_file << "\n";
      }
      for(int i = 0; i < 8; i++) { // edges
        int id = ((i & 3) == 3) ? i-3 : i;
        out_file << id << " " << id + 1 + (i & 1) << "\n";
      }
      for(int i = 0; i < 4; i++) {
        out_file << i << " " << i+4 << "\n";
      }
      out_file.close();
    }

private:

	Node        m_root;
	PointMap    m_points;
    PointRange  m_ranges;
    NormalMap   m_normals;
    Bbox_3      m_bounding_box;
    FT          m_box_length;
    Point       m_center;
    unsigned int m_depth;

}; // end class Octree


} // namespace OCTREE
} // namespace CGAL
