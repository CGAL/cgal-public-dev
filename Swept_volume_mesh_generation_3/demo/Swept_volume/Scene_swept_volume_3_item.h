#ifndef SCENE_SWEEPT_VOLUME_3_ITEM_H
#define SCENE_SWEEPT_VOLUME_3_ITEM_H

#include "Scene_swept_volume_3_item_config.h"
#include <CGAL_demo/Scene_item_with_display_list.h>
#include <iostream>

// This class represents a swept_volume_3 in the OpenGL scene
class SCENE_SWEEPT_VOLUME_3_ITEM_EXPORT Scene_swept_volume_3_item 
  : public Scene_item_with_display_list {
  Q_OBJECT
public:  
  Scene_swept_volume_3_item();
//   Scene_swept_volume_3_item(const Scene_swept_volume_3_item&);
  Scene_swept_volume_3_item(const Swept_volume_3& p);
  Scene_swept_volume_3_item(Swept_volume_3* const p);
  ~Scene_swept_volume_3_item();

  Scene_swept_volume_3_item* clone() const;
  
  // IO
  bool load(std::istream& in);
  bool save(std::ostream& out) const;

  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const { return true; }
  // Points/Wireframe/Flat/Gouraud OpenGL drawing in a display list
  virtual void direct_draw() const;

  // Get wrapped swept_volume_3
  Swept_volume_3*       swept_volume_3();
  const Swept_volume_3* swept_volume_3() const;

  // Get dimensions
  bool isFinite() const { return true; }
  bool isEmpty() const;
  Bbox bbox() const;

private:
  Swept_volume_3* m_swept_volume_3;

}; // end class Scene_swept_volume_3_item

#endif // SCENE_SWEEPT_VOLUME_3_ITEM_H
