
namespace CGAL{

/// \ingroup do_collide_grp
/// @{

/*!
    \brief Returns true if a collision occurs between the two segment trajectories provided

    \details This function returns true for an odd-numbered ray-intersection parity of the computed boundary and the provided ray.
*/
template <class K>
bool do_collide(
    const Segment_3_trajectory<K>& s0,
    const Segment_3_trajectory<K>& s1,
    typename K::Ray_3       test_ray
);

/*!
    \brief Returns true if a collision occurs between the two segment trajectories provided

    \details This function returns true for an odd-numbered ray-intersection parity of the computed boundary and a random ray.
*/
template <class K>
bool do_collide(
    const Segment_3_trajectory<K>& s0,
    const Segment_3_trajectory<K>& s1
);

/// @}


} // end CGAL


