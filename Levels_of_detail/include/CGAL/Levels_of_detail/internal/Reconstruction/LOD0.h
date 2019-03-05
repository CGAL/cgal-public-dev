#ifndef CGAL_LEVELS_OF_DETAIL_0_H
#define CGAL_LEVELS_OF_DETAIL_0_H

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class LOD0 {

  public:
    using Data_structure = DataStructure;
    using Traits = typename Data_structure::Traits;

    LOD0(const Data_structure& data_structure) :
    m_data(data_structure)
    { }

    void reconstruct() {


    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void return_result(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) {

      
    }

  private:
    const Data_structure& m_data;

  }; // LOD0

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_0_H
