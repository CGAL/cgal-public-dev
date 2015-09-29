#include <CGAL/EXACUS_apollonius_traits_2.h>
#include "bench_classes.h"

// EXACUS apollonius
void Bench_EXACUS_apollonius::op()
{
  EXACUS_apollonius_diagram diagram;
  CGAL::lower_envelope_3(_sites.begin(), _sites.end(), diagram);
  
  if (m_verbose_level > 0)
  {
    std::cout << "# of vertices: " << diagram.number_of_vertices() << std::endl;
    std::cout << "# of halfedges: " << diagram.number_of_halfedges() << std::endl;
    std::cout << "# of faces: " << diagram.number_of_faces() << std::endl;
    
    if (m_verbose_level > 1) 
    {
      EXACUS_apollonius_diagram::Vertex_iterator it = 
        diagram.vertices_begin();
      for (; it != diagram.vertices_end(); ++it)
      {
        std::cout << it->point() << std::endl;
      }
    }
  }
}

