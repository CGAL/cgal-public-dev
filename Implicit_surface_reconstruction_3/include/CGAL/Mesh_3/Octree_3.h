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
  
public:
    typedef typename Kernel::FT                 FT;
    typedef typename Kernel::Point_3            Point;
    typedef typename Kernel::Vector_3           Vector;

    typedef typename PointRange::const_iterator InputIterator;
    typedef typename std::list<InputIterator> IterList;

    typedef Octree_node<Kernel,PointRange> Node;

    Octree_node(): 
        m_children(NULL),
        m_parent(NULL),
        m_index(0),
        m_depth(0)
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
    }
    
    void split() 
    {
      m_children = new Node[8];
      for(int i = 0; i < 8; i++) {
        m_children[i].set_parent(this);
        m_children[i].depth() = m_depth + 1;
        m_children[i].index() = i;
        m_children[i].half_size() = (FT) 0.5 * m_half_size;
        
        // barycenter
        FT offsets[3];
        for(int j = 0; j < 3; j++) {
          FT offset = m_children[i].half_size()[j];
          if(((i >> j) & 1) == 0)
            offsets[j] = -offset;
          else
            offsets[j] = +offset;
        }
        m_children[i].barycenter() = m_barycenter + Vector(offsets[0], offsets[1], offsets[2]);
      }
    } 

    bool is_leaf() const { return (m_children == NULL); }
    Node *children() { return m_children; }
    Node *child(const unsigned int index) 
    {
      if(m_children == NULL || index > 7)
        return NULL;
      else
        return &(m_children[index]);
    }
    
    Node *parent() { return m_parent; }
    void set_parent(Node *parent) { m_parent = parent; }

    unsigned char& index() { return m_index; }
    const unsigned char& index() const { return m_index; }
    
    unsigned int& depth() { return m_depth; }
    const unsigned int& depth() const { return m_depth; }
    
    IterList& points() { return m_points; }
    const IterList& points() const { return m_points; }
    void add_point(InputIterator point) { m_points.push_back(point); }
  
    size_t num_points() const { return m_points.size(); }
    bool is_empty() { return (m_points.size() == 0); }
 
    Point& barycenter() { return m_barycenter; }    
    const Point& barycenter() const { return m_barycenter; }

    Vector& half_size() { return m_half_size; }    
    const Vector& half_size() const { return m_half_size; }

private:
    Node*               m_children;
    Node*               m_parent;
    IterList            m_points;
    Point               m_barycenter;
    Vector              m_half_size;
    unsigned char       m_index; // NOT USED YET
    unsigned int        m_depth; // NOT USED YET

}; // end class Octree_node

template < class Kernel,
           class PointRange,
           class PointMap,
           class NormalMap >
class Octree
{   

    typedef typename Kernel::FT            FT;
    typedef typename Kernel::Iso_cuboid_3  Iso_cuboid;
    typedef typename Kernel::Point_3       Point;
    typedef typename Kernel::Vector_3      Vector;
    typedef Octree_node<Kernel, PointRange> Node;
    typedef typename PointRange::const_iterator InputIterator;
    typedef typename std::list<InputIterator> IterList;

public:

    Octree(
        PointRange& pwn, ///< input point range
        PointMap point_map, ///< property map: `value_type of InputIterator` -> `Point` (the position of an input point).
        NormalMap normal_map, ///< property map: `value_type of InputIterator` -> `Vector` (Normal)
        const FT enlarge_ratio = 1.2):
        m_ranges(pwn),
        m_points_map(point_map),
        m_normals_map(normal_map),
        m_num_nodes(0),
        m_num_leafs(0),
        m_num_empty(0)
        
    {
      // compute bbox
      typedef typename PointRange::value_type PointRange_t;
      boost::function<Point(PointRange_t&)> pwn_it_to_point_it = boost::bind(&PointRange_t::first, _1);
      Iso_cuboid bbox = CGAL::bounding_box(boost::make_transform_iterator(pwn.begin(), pwn_it_to_point_it), 
                                           boost::make_transform_iterator(pwn.end(), pwn_it_to_point_it));
      debug_bbox(bbox.min(), bbox.max(), "bbox");
      
      // scale bbox
      Iso_cuboid bbox_scaled = bbox.transform(Aff_transformation_3<Kernel>(SCALING, enlarge_ratio));
      Point bbox_centroid = midpoint(bbox.min(), bbox.max());
      Point bbox_scaled_centroid = midpoint(bbox_scaled.min(), bbox_scaled.max());
      Vector diff_centroid = bbox_centroid - bbox_scaled_centroid;
      bbox_scaled = bbox_scaled.transform(Aff_transformation_3<Kernel>(TRANSLATION, diff_centroid));
      debug_bbox(bbox_scaled.min(), bbox_scaled.max(), "bbox_scaled");
     
      // save octree attributes
      m_bounding_box = bbox_scaled.bbox();
      m_root.barycenter() = bbox_centroid;
      m_root.half_size() = Vector((m_bounding_box.xmax() - m_bounding_box.xmin()) / 2.0,
                                  (m_bounding_box.ymax() - m_bounding_box.ymin()) / 2.0,
                                  (m_bounding_box.zmax() - m_bounding_box.zmin()) / 2.0);
      for (InputIterator it = pwn.cbegin(); it != pwn.cend(); it++)
        m_root.add_point(it);
    }
 
    ~Octree()
    {
      m_root.unsplit();
    }

    // template < typename CellCriteria, typename NormalCriteria > // or other useful criterion 
    void refine(size_t max_depth, 
                size_t max_pts_num) 
                //CellCriteria cell_criteria, 
                //NormalCriteria normal_criteria)
    {
      if(max_depth < 0 || max_pts_num < 1) { 
        CGAL_TRACE_STREAM << "wrong octree refinement criteria\n"; 
        return;
      }
      refine_recurse(&m_root, max_depth, max_pts_num);
      debug_octree("octree_all_nodes", SHOW_ALL_LEAFS); // drawing all leafs is same as all nodes but cheaper
      debug_octree("octree_non_empty_nodes", SHOW_NON_EMPTY_NODES);
      debug_octree("octree_non_empty_leafs", SHOW_NON_EMPTY_LEAFS);
    }
    
    void refine_recurse(Node *node, size_t dist_to_max_depth, size_t max_pts_num) 
    {
      m_num_nodes++;
      if(node->is_empty()) m_num_empty++;
      if (dist_to_max_depth == 0 || node->num_points() <= max_pts_num) {
        m_num_leafs++;
        return;
      }
        
      node->split();
      
      // use parent barycenter to add points in child list
      for (const InputIterator &pwn_it : node->points()) {
        const Point &point = get(m_points_map, *pwn_it);
        int is_right = (node->barycenter()[0] < point[0]);
        int is_up = (node->barycenter()[1] < point[1]);
        int is_front = (node->barycenter()[2] < point[2]);
        int child_id = (is_front << 2) | (is_up << 1) | is_right; 
        node->child(child_id)->add_point(pwn_it);
      }
      
      // recursive calls
      for(int i = 0; i < 8; i++) {
        refine_recurse(node->child(i), dist_to_max_depth-1, max_pts_num);
      }
    }

    void grade()
	{
    }

    template <typename OutputIterator> 
    void generate_points(OutputIterator out)
	{  
    }
    
    // DEBUG
    void debug_header(int num_cuboids, std::ofstream &out_file) 
    {
      out_file << "ply\n" // header
                  "format ascii 1.0\n"
                  "element vertex " 
                  << num_cuboids * 8 << "\n"
                  "property float x\n"
                  "property float y\n"
                  "property float z\n"
                  "element edge " 
                  << num_cuboids * 12 << "\n"
                  "property int vertex1\n"            
                  "property int vertex2\n"
                  "end_header\n";      
    }
    
    void debug_cuboid_vertices(const Point &min, const Point &max, std::ofstream &out_file) 
    {
      for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 3; j++) {
          if(((i >> j) & 1) == 0)
            out_file << min[j] << " ";
          else
            out_file << max[j] << " ";
        }
        out_file << "\n";
      }
    }
    
    void debug_cuboid_vertices(const Point &barycenter, const Vector &half_size, std::ofstream &out_file) 
    {
      debug_cuboid_vertices(barycenter - half_size, barycenter + half_size, out_file);
    }
    
    void debug_cuboid_edges(int cuboid_id, std::ofstream &out_file) 
    {
      for(int i = 0; i < 8; i++) {
        int v1 = (cuboid_id*8) + i;
        if ((i & 3) == 3) v1-=3;
        int v2 = v1 + 1 + (i & 1);
        out_file << v1 << " " << v2 << "\n";
      }
      for(int i = 0; i < 4; i++) {
        int v1 = (cuboid_id*8) + i;
        int v2 = v1 + 4;
        out_file << v1 << " " << v2 << "\n";
      }
    }
    
    void debug_bbox(const Point &min, const Point &max, const std::string &filename) 
    {
      std::cout << "debug bbox " + filename + "\n";
      std::ofstream out_file(filename+".ply");
      debug_header(1, out_file);
      debug_cuboid_vertices(min, max, out_file);
      debug_cuboid_edges(0, out_file);
      out_file.close();
    }

    enum DebugOctreeVisuType 
    {
      SHOW_ALL_LEAFS = 0,
      SHOW_NON_EMPTY_LEAFS = 1,
      SHOW_NON_EMPTY_NODES = 2
    };
    
    void debug_octree(const std::string &filename, DebugOctreeVisuType visu_type) 
    {
      std::cout << "debug octree " + filename + "\n";
      std::ofstream out_file(filename+".ply");
      int num_cuboid_to_draw = 0;
      if(visu_type == SHOW_ALL_LEAFS) {
        num_cuboid_to_draw = m_num_leafs; 
      } else if (visu_type == SHOW_NON_EMPTY_LEAFS) {
        num_cuboid_to_draw = m_num_leafs - m_num_empty;
      } else if (visu_type == SHOW_NON_EMPTY_NODES) {
        num_cuboid_to_draw = m_num_nodes - m_num_empty;
      }
      debug_header(num_cuboid_to_draw, out_file);
      debug_octree_vertices_recursive(&m_root, out_file, visu_type);
      for (int i = 0; i < num_cuboid_to_draw; i++) {
        debug_cuboid_edges(i, out_file);
      }  
      out_file.close();
    }
    
    void debug_octree_vertices_recursive(Node* node, std::ofstream &out_file, DebugOctreeVisuType visu_type) 
    {       
      bool is_leaf = node->is_leaf();
      bool is_non_empty = !node->is_empty();
      if((visu_type == SHOW_ALL_LEAFS && is_leaf) || 
         (visu_type == SHOW_NON_EMPTY_NODES && is_non_empty) ||
         (visu_type == SHOW_NON_EMPTY_LEAFS && is_non_empty && is_leaf)) {
        debug_cuboid_vertices(node->barycenter(), node->half_size(), out_file);
      }   
      if(!is_leaf) {
        for(int i = 0; i < 8; i++) {
          debug_octree_vertices_recursive(node->child(i), out_file, visu_type);
        }    
      }
    }
    
private:
    Node        m_root;
    size_t      m_num_nodes;
    size_t      m_num_leafs;
    size_t      m_num_empty;
    PointRange  m_ranges;
    PointMap    m_points_map;
    NormalMap   m_normals_map;
    Bbox_3      m_bounding_box;
    
}; // end class Octree

} // namespace OCTREE
} // namespace CGAL
