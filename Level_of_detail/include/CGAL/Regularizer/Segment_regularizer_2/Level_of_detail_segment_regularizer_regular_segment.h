#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_H

// STL includes.
#include <map>
#include <list>
#include <cmath>
#include <cassert>

// CGAL includes.
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_regular_segment {

        public:
            typedef KernelTraits Kernel;
            
            using FT      = typename Kernel::FT;
            using Point   = typename Kernel::Point_2;
            using Segment = typename Kernel::Segment_2;
            using Vector  = typename Kernel::Vector_2;

            Level_of_detail_segment_regularizer_regular_segment() : 
            m_orientation(-FT(1)), m_length(-FT(1)), m_is_defined(false) { }

            Level_of_detail_segment_regularizer_regular_segment(const size_t index, const Segment &segment) : 
            m_index(index), m_segment(segment), m_orientation(-FT(1)), m_length(-FT(1)), m_is_defined(true) { 

                compute_orientation();
                compute_length();
                compute_barycentre();
            }

            size_t get_index() const {
                return m_index;
            }

            const Segment &get() const {
                assert(m_is_defined);
                return m_segment;
            }

            FT get_orientation() const {
                assert(m_is_defined);
                return m_orientation;
            }

            FT get_length() const {
                assert(m_is_defined);
                assert(m_length >= FT(0));

                return m_length;
            }

            const Point &get_barycentre() const {
                assert(m_is_defined);
                return m_barycentre;
            }

            const Vector &get_direction() const {
                assert(m_is_defined);
                return m_direction;
            }

            void set_orientation(const FT new_orientation, const FT a, const FT b, const FT c, const Vector &direction) {
                
                // Rotate segment by an angle new_orientation around its barycentre.
                m_orientation = new_orientation;

                // The position of the m_barycentre remains unchanged,
                // which is not the case for the ends. The abscissa of the points are
                // computed according to the new angle new_orientation, however we want to make
                // sure that the final points belong to the line defined by equation
                // _ax + _by + _c = 0, so we are going to use this equation to find
                // the ordinates of the points.
                Vector final_direction = direction;
                if (final_direction.y() < FT(0) || (final_direction.y() == FT(0) && final_direction.x() < FT(0))) 
                    final_direction = -final_direction;

                FT x1, y1, x2, y2;
                if (CGAL::abs(final_direction.x()) > CGAL::abs(final_direction.y())) {
                    
                    x1 = m_barycentre.x() - m_length * final_direction.x() / FT(2);
                    x2 = m_barycentre.x() + m_length * final_direction.x() / FT(2);

                    y1 = (-c - a * x1) / b;
                    y2 = (-c - a * x2) / b;

                } else {
                    y1 = m_barycentre.y() - m_length * final_direction.y() / FT(2);
                    y2 = m_barycentre.y() + m_length * final_direction.y() / FT(2);

                    x1 = (-c - b * y1) / a;
                    x2 = (-c - b * y2) / a;
                }

                const Point source = Point(x1, y1);
                const Point target = Point(x2, y2);

                m_segment = Segment(source, target);
                compute_length();
            }

            void set_reference_coordinates(const Point &new_reference_coordinates) {
                m_reference_coordinates = new_reference_coordinates;
            }

            const Point &get_reference_coordinates() const {
                return m_reference_coordinates;
            }

        private:
            size_t  m_index;
            Segment m_segment;
            FT      m_orientation;
            FT      m_length;
            Point   m_barycentre;
            Vector  m_direction;
            bool    m_is_defined;

            Point m_reference_coordinates;

            void compute_orientation() {
                compute_direction();

                const FT atan = static_cast<FT>(std::atan2(CGAL::to_double(m_direction.y()), CGAL::to_double(m_direction.x())));
                m_orientation = atan * FT(180) / static_cast<FT>(CGAL_PI);

                if (m_orientation < FT(0)) m_orientation += FT(180);
            }

            void compute_length() {
                m_length = static_cast<FT>(CGAL::sqrt(CGAL::to_double(m_segment.squared_length())));
            }

            void compute_barycentre() {
                const FT half = FT(1) / FT(2);

                const Point &source = m_segment.source();
                const Point &target = m_segment.target();

                const FT x = half * (source.x() + target.x());
                const FT y = half * (source.y() + target.y());

                m_barycentre = Point(x, y);
            }

            void compute_direction() {
                
                m_direction = m_segment.to_vector();
                if (m_direction.y() < FT(0) || (m_direction.y() == FT(0) && m_direction.x() < FT(0))) m_direction = -m_direction;
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_H