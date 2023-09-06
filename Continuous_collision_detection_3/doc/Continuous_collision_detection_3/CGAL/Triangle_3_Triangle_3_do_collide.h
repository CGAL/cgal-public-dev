
namespace CGAL{

/// \defgroup do_collide_grp CGAL::do_collide()
/// \ingroup PkgCollisions3Predicates
/// @{

  /*!
     \brief Returns true if a collision occurs between the two triangle trajectories provided

     \details This function determines whether a collision occurs by checking all possible segment-segment and point-triangle collisions that can be inferred from the provided triangle trajectories.
  */
  template <class K>
  bool do_collide(
      const Triangle_3_trajectory<K>& t0,
      const Triangle_3_trajectory<K>& t1
  );

/// @}

}


