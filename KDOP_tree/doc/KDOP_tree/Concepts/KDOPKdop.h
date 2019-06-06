/*!
 * \ingroup PkgKDOPTreeConcepts
 * \cgalConcept
 *
 * The concept 'KDOPKdop' provides types and functions for computing k-dops in the class 'CGAL::KDOP_tree<KDOPTraits>'. The k-dop is represented as a set of support heights in the prescribed k directions.
 *
 * \cgalHasModel 'CGAL::KDOP_kdop<unsigned int N>'
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
     * Value type of the direction vector.
     */
    typedef unspecified_type Direction_type;

    /*!
     * Value type of the vector of k directions.
     */
    typedef std::vector< Direction_type > Vec_direction;

    /*!
     * Value type of k-dop support heights.
     */
    typedef unspecified_type Vec_height;

    /// @}

    /// \name Constructors
    /// K-dop can be constructed by either prescribing the number of directions N
    /// or giving explicitly the directions. If a fixed number N is given, the
    /// default directions corresponding to N are considered.
    /// @{

    /*!
     * Default constructor, the directions are inferred from the prescribed
     * number of directions N. Cartesian axes are considered
     * when N = 6; Cartesian axes + diagonal directions are considered when
     * N = 14, etc.
     * \todo Define default directions for some selected numbers N.
     */
    KDOP_kdop() { }

    /*!
     * Constructor with directions given, the dimension should agree with the
     * prescribed number of directions N. The default directions are not active
     * in this case.
     */
    KDOP_kdop(const Vec_direction& vector_direction)
      : vector_direction_(vector_direction)
    {}

    /// @}

    //TODO some flexibility can be provided for users to choose whether to use user-defined directions or pre-computed directions.

    /// \name Functions
    /// @{

    /*!
     * Return support heights
     */
    inline Vec_height give_support_heights() const;

    /*!
     * Return the minimum support height.
     */
    inline double min_height() const;

    /*!
     * Return the maximum support height.
     */
    inline double max_height() const;

    /*!
     * Add a new direction to vector_direction.
     */
    void add_direction(Direction_type new_direction);

    /*!
     * Check if two k-dops overlap by comparing support heights of the two k-dops.
     */
    inline bool do_overlap(const KDOP& kdop1, const KDOP& kdop2);

    /// @}

  };
