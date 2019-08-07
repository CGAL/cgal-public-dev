/*!
 * \ingroup PkgKDOPTreeConcepts
 * \cgalConcept
 *
 * The concept 'KDOPKdop' provides types and functions for computing k-dops in the class 'CGAL::KDOP_tree<KDOPTraits>'. The k-dop is represented as a set of support heights in the prescribed k directions.
 *
 * \cgalHasModel 'CGAL::KDOP_kdop<typename GeomTraits, unsigned int N>'
 *
 * \sa `CGAL::KDOP_tree<KDOPTraits>`
 *
 */


class KDOPKdop
  {
  public:
    /// \name Types
    /// @{

    /*!
     * Type of the direction.
     */
    typedef unspecified_type Direction_type;

    /*!
     * Type of the vector of k directions.
     */
    typedef std::vector< Direction_type > Vec_direction;

    /*!
     * Type of support heights of a k-DOP.
     */
    typedef unspecified_type Array_height;

    /// Type of support heights for rays
    typedef unspecified_type Array_height_ray;

    /// @}

    /// \name Functions
    /// @{

    /// Set support heights
    unspecified_type set_support_heights(const Array_height& support_heights);

    /*!
     * Return support heights
     */
    inline const Array_height& support_heights();

    /// Return support heights of a ray
    inline const Array_height_ray& support_heights_ray();

    /// @}

    /// \name Overlap detection
    /// @{

    /*!
     * A functor object to check if two k-dops overlap by comparing support heights of the two k-dops. Provides the operators:
     * - bool operator () (const Array_height& support_heights, const unspecified_type & triangle), which returns true if it is possible that the triangle intersects the k-DOP
     * - bool operator () (const Array_height& support_heights, const unspecified_type & ray), which returns true if it is possible that the ray intersects the k-DOP
     * - bool operator () (const Array_height& support_heights, const unspecified_type & squared_radius), which returns true if it is possible that the sphere intersects the k-DOP
     */
    typedef unspecified_type Do_overlap;

    /// @}

    /// \name Operators
    /// @{

    /// Return true if the query is possible to intersect the k-DOP.
    Do_overlap do_overlap_object();

    ///@}

  };
