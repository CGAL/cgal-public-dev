#ifndef CGAL_LEVELS_OF_DETAIL_VEGETATION_H
#define CGAL_LEVELS_OF_DETAIL_VEGETATION_H

// STL includes.
#include <vector>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Vegetation {

  public:
    using Data_structure = DataStructure;

    using Traits = typename DataStructure::Traits;
    using FT = typename Traits::FT;
    
    Vegetation(Data_structure& data_structure) :
    m_data(data_structure)
    { }
    
    void detect_tree_boundaries() {

      if (m_data.verbose) 
        std::cout << "- Detecting tree boundaries" 
        << std::endl;
    }

  private:
    Data_structure& m_data;

  }; // Vegetation

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_VEGETATION_H
