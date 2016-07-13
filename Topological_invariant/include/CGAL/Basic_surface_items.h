#ifndef BASIC_SURFACE_ITEMS_H
#define BASIC_SURFACE_ITEMS_H

#include <CGAL/Generalized_map.h>

namespace CGAL
{
    template<class GMap_items, class Surface>
    struct Basic_surface_GMap_item
    {
        template < class GMap >
        struct Dart_wrapper
        {
            typedef typename GMap_items::template Dart_wrapper<GMap> Old_wrapper;
            struct Dart : public Old_wrapper::Dart{
                typedef typename Path<Surface>::Arc_occurence_handle Arc_occurence_handle;
                Arc_occurence_handle  one_arc;
            };
            typedef typename Old_wrapper::Attributes Attributes;
        };
    };
    
    template<class GMap_items_>
    struct Basic_surface_items
    {
        typedef GMap_items_ GMap_items;
        
        template<class Surface_>
        struct Wrapper
        {
            typedef Surface_ Surface;
            
            typedef Generalized_map<2, Basic_surface_GMap_item<GMap_items, Surface> > GMap;
        };
    };
}

#endif
