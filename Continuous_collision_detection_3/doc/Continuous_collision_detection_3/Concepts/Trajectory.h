/*!
\ingroup PkgCollisions3Concepts
\cgalConcept

The concept `Trajectory` provides a data structure for 
the current and next positions of the points in a geometric
object.

\cgalHasModel `CGAL::Triangle_3_trajectory<Kernel>`
\cgalHasModel `CGAL::Segment_3_trajectory<Kernel>`
\cgalHasModel `CGAL::Point_3_trajectory<Kernel>`

*/
template <class Kernel>
class Trajectory {
public:

  /*!
    The type of object whose trajectory is stored.
   */
  typedef unspecified_type Base_object;

  /*!
    Returns the `Base_object` at the start of its trajectory.
   */
  Base_object current();

  /*!
    Returns the `Base_object` at the end of its trajectory.
   */
  Base_object next();

}; /* end Trajectory */
