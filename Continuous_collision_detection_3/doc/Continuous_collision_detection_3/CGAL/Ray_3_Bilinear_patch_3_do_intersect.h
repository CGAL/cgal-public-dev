
namespace CGAL {

/// \ingroup PkgCollisions3Predicates
/// @{

/*!
  Returns true if the number of ray-bilinear-patch intersections is odd
*/
template <class Kernel>
bool do_intersect_odd_parity(
  const typename CGAL::BilinearPatchC3<Kernel> &bp,
  const typename Kernel::Ray_3 &r
);

/// @}

}
