
namespace CGAL{

/// \ingroup do_collide_grp
/// @{

/*!
    \brief Returns true if a collision occurs between any of the collision meshes provided

    \details This function computes efficiently computes the occurence of a collision by only considering candidate pairs of triangle trajectories, the bounding isocuboids of which intersect.
*/
template <class K>
bool do_collide(
    std::vector< Collision_mesh<K> >& meshes
);

/*!
    \brief Returns true if a collision occurs between the two collision meshes provided

    \details This function computes efficiently computes the occurence of a collision by only considering candidate pairs of triangle trajectories, the bounding isocuboids of which intersect.
*/
template <class K>
bool do_collide(
    Collision_mesh<K>& mesh_1,
    Collision_mesh<K>& mesh_2
);

/// @}

} // end CGAL
#endif