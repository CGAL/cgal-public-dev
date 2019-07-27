#ifndef CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS
#define CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS

// #include <CGAL/license/Shape_regularization.h>

#include <map>
#include <vector>

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits>
  struct Parallel_groups {

  public:
    using Traits = GeomTraits;
    using FT = typename GeomTraits::FT;
    using Groups_type = std::map <FT, std::vector<std::size_t>>;

    Parallel_groups() {}

    void set_parallel_groups(const Groups_type & parallel_groups) {
      m_parallel_groups = parallel_groups;
    }

    Groups_type parallel_groups() {
      return m_parallel_groups;
    }


  private:
    Groups_type m_parallel_groups;

  };

} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS