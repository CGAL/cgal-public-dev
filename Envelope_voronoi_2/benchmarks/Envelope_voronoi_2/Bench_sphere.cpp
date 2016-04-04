#include "bench_classes.h"
#include <CGAL/Envelope_voronoi_2.h>

// Voronoi diagram on a sphere
void Bench_sphere::op()
{
  Envelope_sphere_2 diagram;
  CGAL::voronoi_2(_sites.begin(), _sites.end(), diagram);
  
  if (m_verbose_level > 0)
  {
    std::cout << "# of vertices: " << diagram.number_of_vertices() << std::endl;
    std::cout << "# of halfedges: " << diagram.number_of_halfedges() << std::endl;
    std::cout << "# of faces: " << diagram.number_of_faces() << std::endl;
    
    if (m_verbose_level > 1) 
    {
      Envelope_sphere_2::Vertex_iterator it = diagram.vertices_begin();
      for (; it != diagram.vertices_end(); ++it)
      {
        std::cout << it->point() << std::endl;
      }
    }
  }
}

