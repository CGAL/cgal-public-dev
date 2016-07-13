#ifndef TOPOLOGICAL_SURFACE_H
#define TOPOLOGICAL_SURFACE_H

#include <CGAL/Generalized_map.h>
#include <CGAL/Generalized_map_constructors.h>
#include <CGAL/Path.h>
#include <queue>
#include <map>
#include <list>

#include <CGAL/Basic_surface.h>
#include <CGAL/Basic_surface_items.h>

/* TODO :
 * graph : bonne verif ou pas?
 * 
 * order range
 * order range test
 * 
 * verif add arc hole
 * verif path
 * verif doc 
 * 
 * comment march generic LLC
 * 
 * with point
 * surface file save
 * surface file load
 * path file save
 * path file load
 * 
 * sew
 * 
 * Test découpe poly non cononique
 * 
 * 
 */

namespace CGAL
{
    template<
            typename Items_ = CGAL::Generalized_map_min_items<2> , 
            typename Alloc_ = CGAL_ALLOCATOR(int)
            >
    class Topological_surface : public Basic_surface<Basic_surface_items<Items_>, Alloc_>
    {
    public:
    };
}

#endif //TOPOLOGICAL_SURFACE_H