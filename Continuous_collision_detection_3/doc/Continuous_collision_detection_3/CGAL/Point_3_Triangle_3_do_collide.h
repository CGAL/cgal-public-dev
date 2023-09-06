
namespace CGAL{

/// \ingroup do_collide_grp
/// @{

/*!
    \brief Returns true if a collision occurs between the point trajectory and triangle trajectory provided

    \details This function returns true for an odd-numbered ray-intersection parity of the computed boundary and the provided ray.
*/
template <class K>
bool do_collide(
    const Point_3_trajectory<K>& p,
    const Triangle_3_trajectory<K>& t,
    typename K::Ray_3       test_ray
);

/*!
    \brief Returns true if a collision occurs between the point trajectory and triangle trajectory provided

    \details This function returns true for an odd-numbered ray-intersection parity of the computed boundary and a random ray.
*/
template <class K>
bool do_collide(
    const Point_3_trajectory<K>& p,
    const Triangle_3_trajectory<K>& t
);

/// @}


} // end CGAL


