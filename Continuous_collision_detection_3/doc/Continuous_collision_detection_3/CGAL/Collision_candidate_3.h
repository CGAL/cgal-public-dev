
namespace CGAL{

/*!
\ingroup PkgCollisions3Classes

The class `Collision_candidate` is a container for a pair of primitives that have been identified for potential collision.
*/
template <class Primitive>
class Collision_candidate : public std::pair<const Primitive*, const Primitive*> {

  public:
    /// \name Types
    /// @{

    /*!
    the type from which `Collision_candidate` derives
    */
    using Base = std::pair<const Primitive*, const Primitive*>;

    /*!
    the kernel on which the primitives are built
    */
    using K     = typename Primitive::K;

    /*!
    the index type used to identify the primitives
    */
    using Index = typename Primitive::Index;

    /// @}

    /// \name Creation
    /// @{

    /*!
    constructs a candidate from individual pointers to primitives
    */
    Collision_candidate(const Primitive* first, const Primitive* second) : Base(first, second) {}

    /*!
    constructs a candidate from a pair of primitive pointers
    */
    Collision_candidate(const Base& index_pair) : Base(index_pair) {}

    /// @}
};



}
