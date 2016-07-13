#ifndef BASIC_SURFACE_ITEMS_WITH_POINT_H
#define BASIC_SURFACE_ITEMS_WITH_POINT_H

//#include <CGAL/Linear_cell_complex.h>

namespace CGAL
{
    /*template<class GMap_items, class Surface>
    struct Basic_surface_GMap_item_with_point
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
    
    template<class GMap_items_, class Alloc_>
    struct Basic_surface_items_with_point
    {
        typedef GMap_items_ GMap_items;
        typedef Alloc_ Alloc;
        
        template<class Surface_>
        struct Wrapper
        {
            typedef Surface_ Surface;
            
            typedef Linear_cell_complex<2, 3,  CGAL::Linear_cell_complex_traits<d2>, Basic_surface_GMap_item_with_point<GMap_items, Surface>, Alloc, Generalized_map> GMap;
        };
    };*/
}

#endif
