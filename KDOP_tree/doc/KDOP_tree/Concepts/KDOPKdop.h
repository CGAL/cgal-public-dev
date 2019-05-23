/*!
 * \ingroup PkgKDOPTreeConcepts
 * \cgalConcept
 *
 * The concept 'KDOPKdop' provides types and functions for computing k-dops in the class 'CGAL::KDOP_tree<KDOPTraits>'.
 *
 * \cgalHasModel 'CGAL::KDOP_kdop<unsigned int T>'
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
     * Value type of k-dop directions.
     */
    typedef unspecified_type Vec_direction;

    /*!
     * Value type of k-dop support heights.
     */
    typedef unspecified_type Vec_height;

    /// @}

    /// \name Constructors
    /// @{

    /*!
     * Default constructor, the directions are inferred from the prescribed number of directions N.
     */
    KDOP_kdop() { }

    /*!
     * Constructor with directions given, the dimension should agree with the prescribed number of directions N.
     */
    KDOP_kdop(Vec_direction vector_direction)
      : vector_direction_(vector_direction)
    {}

    /// @}

    //TODO some flexibility can be provided for users to choose whether to use user-defined directions or pre-computed directions.

    /// \name Functions
    /// @{

    /*!
     * Inline function to compute the minimum support height.
     */
    inline double min_height() const { return *std::min_element( vector_height_.begin(), vector_height_.end() ); }

    /*!
     * Inline function to compute the maximum support height
     */
    inline double max_height() const { return *std::max_element( vector_height_.begin(), vector_height_.end() ); }

    /// @}

  };
