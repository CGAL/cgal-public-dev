#include "bench_classes.h"

// CGAL apollonius benchmark
void Bench_apollonius::op()
{
  Apollonius_graph ag;
  ag.insert(_sites.begin(), _sites.end());
  
  if (m_verbose_level > 0)
  {
    std::cout << "# of vertices: " <<  ag.number_of_vertices() << std::endl;
    std::cout << "# of faces: " << ag.number_of_faces() << std::endl;
  }
}

