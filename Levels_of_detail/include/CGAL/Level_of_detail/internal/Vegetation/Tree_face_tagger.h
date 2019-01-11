#ifndef CGAL_LEVEL_OF_DETAIL_TREE_FACE_TAGGER_H
#define CGAL_LEVEL_OF_DETAIL_TREE_FACE_TAGGER_H

// STL includes.
#include <utility>

// CGAL includes.
#include <CGAL/Level_of_detail/internal/Vegetation/Tree.h>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits, class InputRange, class PointMap, class Triangulation>
class Tree_face_tagger {
			
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Triangle_2 Triangle_2;
  typedef typename InputRange::const_iterator const_iterator;
  typedef Tree<GeomTraits, InputRange, PointMap> Tree_item;

  typedef typename Triangulation::Face_handle Face_handle;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Edge Edge;

  using Triangulation_vertex_handle  = typename Triangulation::Vertex_handle;
  
private:

  std::vector<Tree_item>& m_trees;
  PointMap m_point_map;
  Triangulation& m_triangulation;

public:

  Tree_face_tagger(std::vector<Tree_item>& input, PointMap point_map,
                   Triangulation &triangulation)
    : m_trees (input), m_point_map (point_map), m_triangulation(triangulation)
  {
  }

  void tag ()
  {
    for (typename Triangulation::Finite_faces_iterator it = m_triangulation.finite_faces_begin();
         it != m_triangulation.finite_faces_end(); ++ it)
    {
      if (it->info().visibility_label() == Visibility_label::VEGETATION &&
          it->info().group_number() == -1)
        tag_tree (it);
    }
  }


private:

  void tag_tree (Face_handle face)
  {
    Point_2 reference = center (face);
    for (std::size_t i = 0; i < m_trees.size(); ++ i)
    {
      Tree_item& tree = m_trees[i];

      Point_2 center = tree.center_2();
      FT radius = tree.radius();

      if (CGAL::squared_distance (reference, center) < radius * radius)
      {
        flood (face, i);
        break;
      }
    }
    CGAL_assertion (face->info().group_number() != -1);    
  }

  void flood (Face_handle face, int idx)
  {
    std::queue<Face_handle> todo;
    todo.push(face);

    while (!todo.empty())
    {
      Face_handle fh = todo.front();
      todo.pop();

      fh->info().group_number() = idx;
      
      for (std::size_t i = 0; i < 3; ++ i)
      {
        if (m_triangulation.is_constrained (std::make_pair(fh, i)))
          continue;
        
        Face_handle candidate = fh->neighbor(i);

        if (m_triangulation.is_infinite(candidate))
          continue;
        if (candidate->info().group_number() != -1)
          continue;

        todo.push(candidate);
      }
    }
  }

  Point_2 center (Face_handle face) const
  {
    const Point_2& a = face->vertex(0)->point();
    const Point_2& b = face->vertex(1)->point();
    const Point_2& c = face->vertex(2)->point();
    return Point_2 ((a.x() + b.x() + c.x()) / 3.,
                    (a.y() + b.y() + c.y()) / 3.);
  }
    

};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TREE_FACE_TAGGER_H
