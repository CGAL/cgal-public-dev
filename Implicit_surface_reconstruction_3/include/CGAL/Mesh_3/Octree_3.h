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
#include <queue>
#include <vector>

#endif // CGAL_OCTREE_3_H

namespace CGAL {
namespace OCTREE {

// F B U D L R
/*const static bool NU[8][6] = {false, true,  false, true,  false, true,
                              false, true,  false, true,  true,  false,
                              false, true,  true,  false, false, true,
                              false, true,  true,  false, true,  false,
                              true,  false, false, true,  false, true,
                              true,  false, false, true,  true,  false,
                              true,  false, true,  false, false, true,
                              true,  false, true,  false, true,  false};*/

template <class Kernel,
          class PointRange>
class Octree_node
{ 
  
    typedef typename Kernel::FT                 FT;
    typedef typename Kernel::Point_3            Point;
    typedef typename Kernel::Vector_3           Vector;
    typedef typename PointRange::const_iterator InputIterator;
    typedef typename std::list<InputIterator> IterList;
    typedef Octree_node<Kernel,PointRange> Node;
  
public:
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
    
    bool is_sibling(Node *neighbor) {
      return (m_parent == neighbor->parent());
    }

    // dir: LEFT = 0, RIGHT = 1, DOWN = 2, UP = 3, BACK = 4, FRONT= 5
    Node* find_greater_or_equal_neighbor(int dir)
    {
      if(m_parent == NULL) return NULL;
      unsigned int axis_dir = dir & 1;  // 0, 1, 0, 1, 0, 1
      unsigned int axis_bit = dir >> 1; // 0, 0, 1, 1, 2, 2
      unsigned int offset_idx_dir = 1;  // -1, 1, -2, 2, -4, 4
      offset_idx_dir <<= axis_bit;
      if(!axis_dir) offset_idx_dir = -offset_idx_dir;
      for(int i = 0; i < 8; i++) {
        if(((i >> axis_bit) & 1) != axis_dir && m_parent->child(i) == this) { // is 'this' an opposite 'dir' child?
          return m_parent->child(i+offset_idx_dir); // return 'dir' sibling child
        }
      }
      Node *parent_neighbor = m_parent->find_greater_or_equal_neighbor(dir);
      if (parent_neighbor == NULL || parent_neighbor->is_leaf()) return parent_neighbor;
      for(int i = 0; i < 8; i++) { 
        if(((i >> axis_bit) & 1) == axis_dir && m_parent->child(i) == this) { // 'this' is guarantedd to be a 'dir' child
          return parent_neighbor->child(i-offset_idx_dir);  // return opposite 'dir' neighbor child
        }
      }
      return NULL;
    }

    // dir: LEFT = 0, RIGHT = 1, DOWN = 2, UP = 3, BACK = 4, FRONT= 5    
    std::list<Node *> find_smaller_neighbors(Node* ge_neighbor, int dir) 
    {
      std::list<Node *> le_neighbors;
      unsigned int axis_dir = dir & 1;  // 0, 1, 0, 1, 0, 1
      unsigned int axis_bit = dir >> 1; // 0, 0, 1, 1, 2, 2
      std::queue<Node *> possible_neighbors;
      if(ge_neighbor != NULL) possible_neighbors.push(ge_neighbor);
      while (!possible_neighbors.empty()) {
        Node* node = possible_neighbors.front();
        if(node->is_leaf()) {
          le_neighbors.push_back(node);
        }
        else {
          for(int i = 0; i < 8; i++) {
            if(((i >> axis_bit) & 1) != axis_dir) { // add to queue the opposite 'dir' neighbor child
              possible_neighbors.push(node->child(i));    
            }
          }
        }
        possible_neighbors.pop();
      }
      return le_neighbors;
    }
    
    bool is_balanced() 
    { 
      for(int dir = 0; dir < 6; dir++) {
        Node* ge_neighbor = find_greater_or_equal_neighbor(dir);
        std::list<Node *> neighbors = find_smaller_neighbors(ge_neighbor, dir);
        for(Node* neighbor : neighbors) {
          if(neighbor != NULL && !is_sibling(neighbor) 
             && (std::abs((int)m_depth - (int)neighbor->depth()) > 1)) {
            return false;
          }
        }
      }
      return true;
    }
    
    std::list<Node *> find_unbalanced_neighbors_to_split() 
    {
      std::list<Node *> neighbors_to_split;
      for(int dir = 0; dir < 6; dir++) {
        Node* neighbor = find_greater_or_equal_neighbor(dir);
        if(neighbor != NULL && !is_sibling(neighbor) && neighbor->is_leaf() 
           && ((m_depth - neighbor->depth()) > 1)) {
          neighbors_to_split.push_back(neighbor);
        }
      }      
      return neighbors_to_split;
    }

private:
    Node*               m_children;
    Node*               m_parent;
    IterList            m_points;
    Point               m_barycenter;
    Vector              m_half_size;
    unsigned int        m_depth; 
    unsigned char       m_index; // NOT USED YET

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
        PointRange& pwn, 
        PointMap point_map,
        NormalMap normal_map,
        const FT enlarge_ratio = 1.2):
        m_ranges(pwn),
        m_points_map(point_map),
        m_normals_map(normal_map)
        
    {
      // compute bbox
      typedef typename PointRange::value_type PointRange_t;
      boost::function<Point(PointRange_t&)> pwn_it_to_point_it = boost::bind(&PointRange_t::first, _1);
      Iso_cuboid bbox = CGAL::bounding_box(boost::make_transform_iterator(pwn.begin(), pwn_it_to_point_it), 
                                           boost::make_transform_iterator(pwn.end(), pwn_it_to_point_it));
      dump_bbox(bbox.min(), bbox.max(), "bbox");
      
      // scale bbox
      Iso_cuboid bbox_scaled = bbox.transform(Aff_transformation_3<Kernel>(SCALING, enlarge_ratio));
      Point bbox_centroid = midpoint(bbox.min(), bbox.max());
      Point bbox_scaled_centroid = midpoint(bbox_scaled.min(), bbox_scaled.max());
      Vector diff_centroid = bbox_centroid - bbox_scaled_centroid;
      bbox_scaled = bbox_scaled.transform(Aff_transformation_3<Kernel>(TRANSLATION, diff_centroid));
      dump_bbox(bbox_scaled.min(), bbox_scaled.max(), "bbox_scaled");
     
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
      dump_octree("octree_all_nodes", SHOW_ALL_LEAFS); // drawing all leafs is same as all nodes but cheaper
      dump_octree("octree_non_empty_nodes", SHOW_NON_EMPTY_NODES);
      dump_octree("octree_non_empty_leafs", SHOW_NON_EMPTY_LEAFS);
    }
    
    void refine_recurse(Node *node, size_t dist_to_max_depth, size_t max_pts_num) 
    {
      if(dist_to_max_depth == 0 || node->num_points() <= max_pts_num)
        return;
    
      node->split();
      reassign_points(node);
      for(int i = 0; i < 8; i++) {
        refine_recurse(node->child(i), dist_to_max_depth-1, max_pts_num);
      }
    }
    
    void reassign_points(Node *node){
      for (const InputIterator &pwn_it : node->points()) {
        const Point &point = get(m_points_map, *pwn_it);
        int is_right = (node->barycenter()[0] < point[0]);
        int is_up = (node->barycenter()[1] < point[1]);
        int is_front = (node->barycenter()[2] < point[2]);
        int child_id = (is_front << 2) | (is_up << 1) | is_right; 
        node->child(child_id)->add_point(pwn_it);
      }
    }

    void grade()
	{
      std::queue<Node *> leaf_nodes;
      fill_leaf_queue(&m_root, leaf_nodes);
      while(!leaf_nodes.empty()) {
        Node *node = leaf_nodes.front();
        leaf_nodes.pop();
        if(!node->is_leaf()) continue;
        std::list<Node *> neighbors_to_split = node->find_unbalanced_neighbors_to_split();
        if(!neighbors_to_split.empty()) leaf_nodes.push(node);
        for(Node* neighbor : neighbors_to_split) {
          neighbor->split();
          reassign_points(neighbor);
          for(int i = 0; i < 8; i++) {
            Node* neighbor_child = neighbor->child(i);
            leaf_nodes.push(neighbor_child);
          }
        }
      }
      dump_octree("balanced_octree_all_nodes", SHOW_ALL_LEAFS);
      dump_octree("balanced_octree_non_empty_nodes", SHOW_NON_EMPTY_NODES);
      dump_octree("balanced_octree_non_empty_leafs", SHOW_NON_EMPTY_LEAFS);
      (debug_grading(&m_root)) ? std::cout << "octree correctly balanced!\n" : 
                                 std::cerr << "Error: octree not correctly balanced!\n"; 
    }
    
    void fill_leaf_queue(Node *node, std::queue<Node *> &queue) {
      if (node->is_leaf()) queue.push(node);
      else {
        for(int i = 0; i < 8; i++) {
          fill_leaf_queue(node->child(i), queue);
        }     
      }   
    }

    template <typename OutputIterator> 
    void generate_points(OutputIterator out)
	{  
    }
    
    // DEBUG functions
    void dump_header(int num_cuboids, std::ofstream &out_file) 
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
    
    void dump_cuboid_vertices(const Point &min, const Point &max, std::ofstream &out_file) 
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
    
    void dump_cuboid_vertices(const Point &barycenter, const Vector &half_size, std::ofstream &out_file) 
    {
      dump_cuboid_vertices(barycenter - half_size, barycenter + half_size, out_file);
    }
    
    void dump_cuboid_edges(int cuboid_id, std::ofstream &out_file) 
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
    
    void dump_bbox(const Point &min, const Point &max, const std::string &filename) 
    {
      std::cout << "dump bbox " + filename + "\n";
      std::ofstream out_file(filename+".ply");
      dump_header(1, out_file);
      dump_cuboid_vertices(min, max, out_file);
      dump_cuboid_edges(0, out_file);
      out_file.close();
    }

    enum DebugOctreeVisuType 
    {
      SHOW_ALL_LEAFS = 0,
      SHOW_NON_EMPTY_LEAFS = 1,
      SHOW_NON_EMPTY_NODES = 2
    };
    
    void dump_octree(const std::string &filename, DebugOctreeVisuType visu_type) 
    {
      std::cout << "dump octree " + filename + "\n";
      std::ofstream out_file(filename+".ply");
      int num_cuboid_to_draw = 0; 
      get_num_cuboid_to_draw(&m_root, num_cuboid_to_draw, visu_type);
      dump_header(num_cuboid_to_draw, out_file);
      dump_octree_vertices_recursive(&m_root, out_file, visu_type);
      for (int i = 0; i < num_cuboid_to_draw; i++) {
        dump_cuboid_edges(i, out_file);
      }  
      out_file.close();
    }
    
    void get_num_cuboid_to_draw(Node* node, int &num_cuboid_to_draw, DebugOctreeVisuType visu_type) 
    {
      bool is_leaf = node->is_leaf();
      bool is_non_empty = !node->is_empty();
      if((visu_type == SHOW_ALL_LEAFS && is_leaf) || 
         (visu_type == SHOW_NON_EMPTY_NODES && is_non_empty) ||
         (visu_type == SHOW_NON_EMPTY_LEAFS && is_non_empty && is_leaf)) {
        num_cuboid_to_draw++;     
      }
      if(!is_leaf) {
        for(int i = 0; i < 8; i++) {
          get_num_cuboid_to_draw(node->child(i), num_cuboid_to_draw, visu_type);
        }  
      }
    }
    
    void dump_octree_vertices_recursive(Node* node, std::ofstream &out_file, DebugOctreeVisuType visu_type) 
    {       
      bool is_leaf = node->is_leaf();
      bool is_non_empty = !node->is_empty();
      if((visu_type == SHOW_ALL_LEAFS && is_leaf) || 
         (visu_type == SHOW_NON_EMPTY_NODES && is_non_empty) ||
         (visu_type == SHOW_NON_EMPTY_LEAFS && is_non_empty && is_leaf)) {
        dump_cuboid_vertices(node->barycenter(), node->half_size(), out_file);
      }   
      if(!is_leaf) {
        for(int i = 0; i < 8; i++) {
          dump_octree_vertices_recursive(node->child(i), out_file, visu_type);
        }    
      }
    }
    
    bool debug_grading(Node* node) {
      if(node->is_leaf()) {
        return node->is_balanced();
      }
      for(int i = 0; i < 8; i++) {
        if(!debug_grading(node->child(i)))
          return false;
      }
      return true;
    }
    // end DEBUG functions
    
private:
    Node        m_root;         /* root node of the octree */
    PointRange  m_ranges;       /* input point range */
    PointMap    m_points_map;   /* property map: `value_type of InputIterator` -> `Point` (Position) */
    NormalMap   m_normals_map;  /* property map: `value_type of InputIterator` -> `Vector` (Normal) */
    Bbox_3      m_bounding_box; /* octree bounding box */
    
}; // end class Octree

} // namespace OCTREE
} // namespace CGAL
