#include "bench_classes.h"

// cgal triangulation
void Bench_triangulation_voronoi::op()
{
  Triangulation_voronoi vd;
  vd.insert(_sites.begin(), _sites.end());
  
  if (m_verbose_level > 0)
  {
    std::cout << "# of vertices: " << vd.number_of_vertices() << std::endl;
    std::cout << "# of halfedges: " << vd.number_of_halfedges() << std::endl;
    std::cout << "# of faces: " << vd.number_of_faces() << std::endl;
    
    if (m_verbose_level > 1) 
    {
      Triangulation_voronoi::Vertex_iterator it = vd.vertices_begin();
      for (; it != vd.vertices_end(); ++it)
      {
        std::cout << it->point() << std::endl;
      }
    }
  }
}

