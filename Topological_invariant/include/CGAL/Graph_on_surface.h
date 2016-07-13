#ifndef GRAPH_ON_SURFACE_H
#define GRAPH_ON_SURFACE_H



#include <CGAL/Compact_container.h>

namespace CGAL
{
    template<class Surface>
    class Graph_on_surface : public Compact_container_base {
    public:
        Graph_on_surface(Surface& surface):
        mSurface(surface)
        {
        }
        
    private:
        Surface& mSurface;
    };
};

#endif
