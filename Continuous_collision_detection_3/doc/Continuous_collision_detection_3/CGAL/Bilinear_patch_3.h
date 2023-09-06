
namespace CGAL {

/*!
\ingroup PkgCollisions3Classes

The class `CGAL::BilinearPatchC3` provides a bilinear-patch implementation for cartesian
and related kernels. The patch is defined by the bilinear interpolation between four points
in space, \f$ { x_0, x_1, x_2, x_3 \in \mathbb{R}^3 }\f$. Any point on the surface can be 
expressed as \f$ { p(u,v) = (1-u)(1-v)x_0 + u(1-v)x_1 + (1-u)vx_2 + uvx_3  }\f$. 

*/
template <class Kernel>
class BilinearPatchC3
{

  public: 
    /// \name Types
    /// @{

    /*!
    the underlying kernel
    */
    typedef Kernel                            R;

    /*!
    the `FieldNumberType` used
    */
    typedef typename R::FT                    FT;

    /*!
    the `Point_3` type
    */
    typedef typename R::Point_3               Point_3;

    /*!
    the `Vector_3` type
    */
    typedef typename R::Vector_3              Vector_3;

    /*!
    the `Plane_3` type
    */
    typedef typename R::Plane_3               Plane_3;

    /*!
    the `Triangle_3` type
    */
    typedef typename R::Triangle_3            Triangle_3;

    /*!
    the `Segment_3` type
    */
    typedef typename R::Segment_3             Segment_3;

    /*!
    the `Tetrahedron_3` type
    */
    typedef typename R::Tetrahedron_3         Tetrahedron_3;

    /// @}

    // ===================================================================

    /// \name Creation
    /// @{

    /*!
    creates an empty bilinear patch
    */
    BilinearPatchC3() {}


    /*!
    Creates a bilinear patch for which pqsrp is a valid cycle around the boundary
    */
    BilinearPatchC3(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s);

    /// @}

    // ===================================================================

    /// \name Operators
    /// @{

    /*!
    returns a point corresponding to the parametric coordinates \f$ { (u,v) } \f$ on the bilinear patch.
    */
    Point_3 operator()(const FT& u, const FT& v) const;

    /*!
    returns true if the vertices are the same and are ordered in a way that produces an equivalent cycle
    */
    bool  operator==(const BilinearPatchC3& bp) const;

    /*!
    returns true if `operator==()` is false
    */
    bool  operator!=(const BilinearPatchC3& bp) const;
    
    /// @}

    // ===================================================================

    /// \name Predicates
    /// @{

    /*!
    returns the orientation of the point with respect to the bilinear patch as determined by `signed_scaled_patch_distance()`
    */
    ::CGAL::Orientation orientation(const Point_3& p) const;

    /*!
    returns true if the signed distance from the bilinear is patch is zero
    */
    bool  has_on(const Point_3& p) const;

    /*!
    returns true if the bilinear patch is equivalent to a collection of segments or a point
    */
    bool  is_degenerate() const;

    /*!
    returns true if the all corners of the bilinear patch are coplanar
    */
    bool  is_planar() const;

    /// @}

    // ===================================================================

    /// \name Methods
    /// @{

    /*!
    returns the `Point_3` object corresponding to the i-th index
    */
    const Point_3 & vertex(int i) const;
    
    /*!
    returns the `Point_3` object corresponding to the i-th index
    */
    const Point_3 & operator[](int i) const;

    /*!
    returns the `Tetrahedron_3` object whose vertices coincide with the vertices of the bilinear patch
    */
    const Tetrahedron_3 & tetrahedron() const;

    /*!
    returns a signed, scaled distance between the point and the bilinear patch, which can be used for orientation and `has_on()`
    */
    FT signed_scaled_patch_distance(const Point_3& x) const;

    /// @}
};

} //namespace CGAL

