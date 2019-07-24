#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_DATA_2
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_DATA_2

// #include <CGAL/license/Shape_regularization.h>

#include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<
    typename GeomTraits>
  struct Segment_data_2 {

  public:
    using Traits = GeomTraits;
    using Segment = typename GeomTraits::Segment_2;
    using Point = typename GeomTraits::Point_2;
    using Vector  = typename GeomTraits::Vector_2;
    using FT = typename GeomTraits::FT;

    const std::size_t  m_index;
    const Segment& m_segment;
    const FT m_angle;
    Vector  m_direction;
    FT      m_orientation;
    Point   m_barycentre;
    FT      m_length;
    Point   m_reference_coordinates;
    FT      m_a;
    FT      m_b;
    FT      m_c;


    Segment_data_2(
      const Segment& segment,
      const std::size_t index,
      const FT angle = FT(-1000)):
    m_segment(segment),
    m_index(index),
    m_angle(angle) {

      m_barycentre = compute_middle_point(m_segment.source(), m_segment.target());
      m_length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(m_segment.squared_length())));

      if(m_angle != FT(-1000)) {
        const FT x = static_cast<FT>(cos(CGAL::to_double(m_angle * static_cast<FT>(CGAL_PI) / FT(180))));
        const FT y = static_cast<FT>(sin(CGAL::to_double(m_angle * static_cast<FT>(CGAL_PI) / FT(180))));

        m_direction = Vector(x, y);
        const Vector v_ort = Vector(-m_direction.y(), m_direction.x());
        
        m_a = v_ort.x();
        m_b = v_ort.y();
        m_c = -m_a * m_barycentre.x() - m_b * m_barycentre.y();
      }
      else {
        m_direction = compute_direction(m_segment);
        m_orientation = compute_orientation(m_direction);
        m_a = -static_cast<FT>(sin(CGAL::to_double(m_orientation) * CGAL_PI / 180.0));
        m_b =  static_cast<FT>(cos(CGAL::to_double(m_orientation) * CGAL_PI / 180.0));
        m_c = -m_a * m_barycentre.x() - m_b * m_barycentre.y();

      }
      std::cout << "Segment #" << m_index << ": direction = " << m_direction << "; a = " << m_a 
                << "; b = " << m_b << "; c = " << m_c << ";" << std::endl;

    }

  private:
    // FT      m_difference;

  };

} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_DATA_2
