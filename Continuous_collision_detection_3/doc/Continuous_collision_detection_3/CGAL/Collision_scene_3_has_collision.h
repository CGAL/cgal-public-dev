
namespace CGAL{


/// \ingroup PkgCollisions3Predicates
/// @{

/*!
  Returns true if any pair of triangle trajectories in the scene collides.
*/
template <class K>
bool has_collision(
    const Collision_scene<K>& scene
);


/*!
  Returns true if the pair of triangle trajectories contained by the `Collision_candidate` collides.
*/
template <class CollisionCandidate>
bool candidate_has_collision(
    const CollisionCandidate& candidate
);


/// @}


/// \ingroup PkgCollisions3Functions
/// @{


/*!
  \brief Returns a vector of `Collision_candidate` objects, each containing a pair of triangle trajectories that were confirmed to collide.
*/
template <class K>
std::vector<Collision_candidate<typename Collision_scene<K>::Trajectory>>
get_collisions(
    const Collision_scene<K>& scene
);

/*!
  Returns a vector of `Collision_candidate` objects, each containing a pair of triangle trajectories that were flagged as possibly colliding.
*/
template <class K>
auto get_collision_candidates(
  const Collision_scene<K>& scene
) -> std::vector< Collision_candidate<typename Collision_scene<K>::Trajectory> >;


/// @}

} // end CGAL
