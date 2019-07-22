#include "include/isr_test_types.h"

#include <CGAL/Polygon_mesh_processing/connected_components.h>

typedef Mesh::Halfedge_index halfedge_descriptor;

size_t nb_boundaries(const Mesh &m)
{
    size_t nb = 0;
    std::set<halfedge_descriptor> he_set;

    BOOST_FOREACH(halfedge_descriptor hd, m.halfedges())
    {
      if(m.is_border(hd) && (he_set.find(hd) == he_set.end()))
      {
        nb++;
        halfedge_descriptor curr = hd;
        do
        {
          curr = m.next(curr);
          he_set.insert(curr);
        }
        while(curr != hd);
      }
    }
    return nb;
}