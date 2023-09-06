/*!
\ingroup PkgCollisions3Concepts
\cgalConcept

The concept `CollisionSceneUpdateFunctor` defines the requirements for any functor used to update a collision scene
with the method `CGAL::Collision_scene_3<Kernel>::update_state`.

\cgalHasModel `CGAL::Collisions::internal::Translate_functor<CollisionScene>`
\cgalHasModel `CGAL::Collisions::internal::Swap_current_next_functor<CollisionScene>`
\cgalHasModel `CGAL::Collisions::internal::Contraction_functor<CollisionScene>`

*/

template <class CollisionScene>
class CollisionSceneUpdateFunctor {
public:

/*!
  Applies the user-defined update to the vertex characterized by the `Scene_vertex_index` argument.
*/
void operator() (CollisionScene::Mesh* mesh, CollisionScene::Scene_vertex_index svi);

}; /* end CollisionSceneUpdateFunctor */

