#ifndef CGAL_LEVEL_OF_DETAIL_WALL_TO_TRIANGLES_H
#define CGAL_LEVEL_OF_DETAIL_WALL_TO_TRIANGLES_H

// STL includes.
#include <vector>
#include <fstream>

// CGAL includes
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// LOD includes.
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits, class InputTriangulation>
class Wall_to_triangles
{

public:
  using Kernel        = GeomTraits;

  using Point_2 = typename Kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Plane_3 = typename Kernel::Plane_3;
  using FT = typename Kernel::FT;

  typedef CGAL::Triangulation_vertex_base_with_info_2<Point_3, Kernel> Vertex_base;
  typedef CGAL::Triangulation_face_base_2<Kernel> Face_base;
  typedef CGAL::Triangulation_data_structure_2<Vertex_base, Face_base> TDS;
  typedef CGAL::Delaunay_triangulation_2<Kernel, TDS> Triangulation;

private:
  
  const InputTriangulation& m_input_triangulation;
  
public:

  Wall_to_triangles (const InputTriangulation& triangulation)
    : m_input_triangulation (triangulation)
  {
  }

  void compute (const typename InputTriangulation::Edge& e,
                std::vector<Triangle_3>& triangles)
  {
    typename InputTriangulation::Face_handle
      f0 = e.first,
      f1 = e.first->neighbor(e.second);

    typename InputTriangulation::Vertex_handle
      va = f0->vertex ((e.second + 1) % 3),
      vb = f0->vertex ((e.second + 2) % 3);

    int va0 = f0->index(va);
    int va1 = f1->index(va);
    int vb0 = f0->index(vb);
    int vb1 = f1->index(vb);

    std::set<FT> za;
    get_heights_of_vertex (va, za);
    std::set<FT> zb;
    get_heights_of_vertex (vb, zb);

    FT za0 = internal::point_3<Point_3>(f0, va0).z();
    FT za1 = internal::point_3<Point_3>(f1, va1).z();
    FT zb0 = internal::point_3<Point_3>(f0, vb0).z();
    FT zb1 = internal::point_3<Point_3>(f1, vb1).z();
    
    FT za_min = (std::min)(za0, za1);
    FT za_max = (std::max)(za0, za1);
    FT zb_min = (std::min)(zb0, zb1);
    FT zb_max = (std::max)(zb0, zb1);

    cpp11::array<Point_3, 3> corners = {{ Point_3 (va->point().x(), va->point().y(), za_min),
                                          Point_3 (va->point().x(), va->point().y(), za_max),
                                          Point_3 (vb->point().x(), vb->point().y(), zb_max) }};

    if (za0 > za1 || zb0 > zb1) // orient normal consistently with 3D mesh
      std::swap (corners[1], corners[2]);

    Vector_3 normal = CGAL::normal (corners[0], corners[1], corners[2]);
    Plane_3 plane (corners[0], normal);

    Triangulation triangulation;

    for (typename std::set<FT>::iterator it = za.begin(); it != za.end(); ++ it)
    {
      FT height = *it;
      if (height < za_min || height > za_max)
        continue;

      Point_3 p = Point_3 (va->point().x(), va->point().y(), height);
      typename Triangulation::Vertex_handle
        v = triangulation.insert (plane.to_2d (p));
      v->info() = p;
    }

    for (typename std::set<FT>::iterator it = zb.begin(); it != zb.end(); ++ it)
    {
      FT height = *it;
      if (height < zb_min || height > zb_max)
        continue;

      Point_3 p = Point_3 (vb->point().x(), vb->point().y(), height);
      typename Triangulation::Vertex_handle
        v = triangulation.insert (plane.to_2d (p));
      v->info() = p;
    }

    for (typename Triangulation::Finite_faces_iterator it = triangulation.finite_faces_begin();
         it != triangulation.finite_faces_end(); ++ it)
    {
      triangles.push_back (Triangle_3 (it->vertex(0)->info(),
                                       it->vertex(1)->info(),
                                       it->vertex(2)->info()));
    }
  }

private:

  void get_heights_of_vertex (typename InputTriangulation::Vertex_handle v, std::set<FT>& heights)
  {
    typename InputTriangulation::Face_circulator
      start = m_input_triangulation.incident_faces(v),
      circ = start;

    do
    {
      if (!m_input_triangulation.is_infinite(circ))
        heights.insert (internal::point_3<Point_3>(circ, circ->index(v)).z());

      ++ circ;
    }
    while (circ != start);
  }

};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_WALL_TO_TRIANGLES_H
