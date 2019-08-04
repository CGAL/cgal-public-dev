#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_CONDITIONS_ANGLES_2
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_CONDITIONS_ANGLES_2

// #include <CGAL/license/Shape_regularization.h>

#include <CGAL/Shape_regularization/internal/Segment_data_2.h>
#include <map>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<
    typename GeomTraits>
  class Conditions_angles_2 {

    public:
      using Traits = GeomTraits;
      using FT = typename GeomTraits::FT;
      using Segment_data = typename internal::Segment_data_2<Traits>;

      Conditions_angles_2() :
      m_moe(FT(1) / FT(4)) {}


     FT reference(const Segment_data & seg_data, const FT suffix) const {
      FT val = seg_data.m_orientation + suffix; 

      if (val < FT(0)) val += FT(180); 
      else if (val > FT(180)) val -= FT(180);

      return val;
     }

    int group_index(const FT val, const FT val_j, const int g_index) const {
      int g_j = -1;
      for (int k = -1; k <= 1; ++k) {  
        if (CGAL::abs(val_j - val + static_cast<FT>(k) * FT(180)) < m_moe) {  
          g_j = g_index;
          break;  
        }
      }
      return g_j;
    }

    FT get_margin_of_error() const {
      return m_moe;
    }

    void set_margin_of_error(const FT max_bound) {
      CGAL_precondition(max_bound > 0);
      m_moe = max_bound / FT(100);
    }

    private:
      FT m_moe;

    };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_CONDITIONS_ANGLES_2
