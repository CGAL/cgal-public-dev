#ifndef SCENE_SWEEPT_VOLUME_3_ITEM_CONFIG_H
#define SCENE_SWEEPT_VOLUME_3_ITEM_CONFIG_H

#ifdef scene_swept_volume_3_item_EXPORTS
#  define SCENE_SWEEPT_VOLUME_3_ITEM_EXPORT Q_DECL_EXPORT
#else
#  define SCENE_SWEEPT_VOLUME_3_ITEM_EXPORT Q_DECL_IMPORT
#endif

//#include <SV/Swept_volume_with_bvh_3.h> 
//#include <SV/Swept_volume_with_delaunay_3.h> 
//#include <SV/Swept_volume_with_octree_3.h> 
#include <SV/Swept_volume_with_vhull_3.h> 


//typedef Swept_volume_with_bvh_3      Swept_volume_3; 
//typedef Swept_volume_with_delaunay_3 Swept_volume_3; 
//typedef Swept_volume_with_octree_3   Swept_volume_3; 
typedef SV::Swept_volume_with_vhull_3<>    Swept_volume_3; 



#include "Scene_swept_volume_3_item.h"


#endif // SCENE_SWEEPT_VOLUME_3_ITEM_CONFIG_H
