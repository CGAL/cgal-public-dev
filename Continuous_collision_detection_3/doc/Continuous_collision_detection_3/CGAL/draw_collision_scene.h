
namespace CGAL {

/// \ingroup PkgCollisions3Functions
/// @{

/*!
  Draws all the meshes in the scene. If `draw_next` is true, the meshes will be drawn with their points at their "next" location, i.e., at the end of the time period under consideration, and any triangles that collide will be colored white.
*/
template <class CollisionScene>
void draw_collision_scene( CollisionScene& scene, bool draw_next=false );

/// @}

}
