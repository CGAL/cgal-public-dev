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

        template<class Kernel>
        class Level_of_detail_segment_regularizer_tree_parallel_segments_node;

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_regular_segment {

        public:
            typedef KernelTraits Kernel;
            
            using FT      = typename Kernel::FT;
            using Point   = typename Kernel::Point_2;
            using Segment = typename Kernel::Segment_2;
            using Vector  = typename Kernel::Vector_2;

            using Parallel_segments_tree_node = Level_of_detail_segment_regularizer_tree_parallel_segments_node<Kernel>;

            Level_of_detail_segment_regularizer_regular_segment() : 
            m_orientation(-FT(1)), m_difference(-FT(1)), m_length(-FT(1)), m_is_defined(false) { }

            Level_of_detail_segment_regularizer_regular_segment(const size_t index, const Segment &segment) : 
            m_index(index), m_segment(segment), m_orientation(-FT(1)), m_difference(-FT(1)), m_length(-FT(1)), m_is_defined(true) { 

                compute_orientation();
                compute_difference();

                compute_length();
                compute_barycentre();
                compute_line_coefficients();

                compute_reference_coordinates();
                parallel_node = NULL;
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

            FT get_difference() const {
                assert(m_is_defined);
                return m_difference;
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

            const Point &get_reference_coordinates() const {
                return m_reference_coordinates;
            }

            FT get_a() const {
                return m_a;
            }

            FT get_b() const {
                return m_b;
            }

            FT get_c() const {
                return m_c;
            }

            void set_reference_coordinates(const Point &new_reference_coordinates) {
                m_reference_coordinates = new_reference_coordinates;
            }

            void set_orientation(const FT new_orientation, const FT a, const FT b, const FT c, const Vector &direction) {
                
                // Rotate segment by an angle new_orientation around its barycentre.
                m_orientation = new_orientation;

                // We update the equation of the support line.
                m_a = a;
                m_b = b;
                m_c = c;

                // The position of the m_barycentre remains unchanged,
                // which is not the case for the ends. The abscissa of the points are
                // computed according to the new angle new_orientation, however we want to make
                // sure that the final points belong to the line defined by equation
                // _ax + _by + _c = 0, so we are going to use this equation to find
                // the ordinates of the points.
                m_direction = direction;
                if (m_direction.y() < FT(0) || (m_direction.y() == FT(0) && m_direction.x() < FT(0))) 
                    m_direction = -m_direction;

                FT x1, y1, x2, y2;
                if (CGAL::abs(m_direction.x()) > CGAL::abs(m_direction.y())) {
                    
                    x1 = m_barycentre.x() - m_length * m_direction.x() / FT(2);
                    x2 = m_barycentre.x() + m_length * m_direction.x() / FT(2);

                    y1 = (-m_c - m_a * x1) / m_b;
                    y2 = (-m_c - m_a * x2) / m_b;

                } else {
                    y1 = m_barycentre.y() - m_length * m_direction.y() / FT(2);
                    y2 = m_barycentre.y() + m_length * m_direction.y() / FT(2);

                    x1 = (-m_c - m_b * y1) / m_a;
                    x2 = (-m_c - m_b * y2) / m_a;
                }

                const Point source = Point(x1, y1);
                const Point target = Point(x2, y2);

                m_segment = Segment(source, target);
                compute_length();
            }

            void set_difference(const FT new_difference, const FT a, const FT b, const FT c, const Vector &direction) {

                // We translate the segment by the distance new_difference in the direction of the normal vector.
                m_difference = new_difference;

                // We update the equation of the support line.
                m_a = a;
                m_b = b;
                m_c = c;

                m_direction = direction;
                if (m_direction.y() < FT(0) || (m_direction.y() == FT(0) && m_direction.x() < FT(0))) m_direction = -m_direction;

                Vector final_normal = Vector(-m_direction.y(), m_direction.x());
                FT bx, by, x1, x2, y1, y2;

                const Point &source = m_segment.source();
                const Point &target = m_segment.target();

                if (CGAL::abs(m_direction.x()) > CGAL::abs(m_direction.y())) {
                    bx = m_barycentre.x() + m_difference * final_normal.x();

                    x1 = source.x() + m_difference * final_normal.x();
                    x2 = target.x() + m_difference * final_normal.x();

                    by = (-m_c - m_a * m_barycentre.x()) / m_b;

                    y1 = (-m_c - m_a * x1) / m_b;
                    y2 = (-m_c - m_a * x2) / m_b;

                } else {

                    by = m_barycentre.y() + m_difference * final_normal.y();
                    
                    y1 = source.y() + m_difference * final_normal.y();
                    y2 = target.y() + m_difference * final_normal.y();

                    bx = (-m_c - m_b * m_barycentre.y()) / m_a;

                    x1 = (-m_c - m_b * y1) / m_a;
                    x2 = (-m_c - m_b * y2) / m_a;
                }

                const Point new_source = Point(x1, y1);
                const Point new_target = Point(x2, y2);

                m_segment    = Segment(new_source, new_target);
                m_barycentre = Point(bx, by);

                compute_length();
            }

            void set_difference(const FT new_difference) {
                
                m_difference        = new_difference;
                Vector final_normal = Vector(-m_direction.y(), m_direction.x());

                const Point &source = m_segment.source();
                const Point &target = m_segment.target();

                Point new_source = Point(source.x() + m_difference * final_normal.x(), source.y() + m_difference * final_normal.y());
                Point new_target = Point(target.x() + m_difference * final_normal.x(), target.y() + m_difference * final_normal.y());

                const FT bx = (new_source.x() + new_target.x()) / FT(2);
                const FT by = (new_source.y() + new_target.y()) / FT(2);

                m_segment    = Segment(new_source, new_target);
                m_barycentre = Point(bx, by);

                m_c = -m_a * bx - m_b * by;
                compute_length();
            }

            Parallel_segments_tree_node *parallel_node;

        private:
            size_t  m_index;
            Segment m_segment;
            FT      m_orientation;
            FT      m_difference;
            FT      m_length;
            Point   m_barycentre;
            Vector  m_direction;
            bool    m_is_defined;
            
            FT m_a, m_b, m_c;
            Point m_reference_coordinates;

            void compute_orientation() {
                compute_direction();

                const FT atan = static_cast<FT>(std::atan2(CGAL::to_double(m_direction.y()), CGAL::to_double(m_direction.x())));
                m_orientation = atan * FT(180) / static_cast<FT>(CGAL_PI);

                if (m_orientation < FT(0)) m_orientation += FT(180);
            }

            void compute_difference() {
                m_difference = FT(0);
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

            void compute_line_coefficients() {
	            
                m_a = -static_cast<FT>(sin(CGAL::to_double(m_orientation) * CGAL_PI / 180.0));
	            m_b =  static_cast<FT>(cos(CGAL::to_double(m_orientation) * CGAL_PI / 180.0));
	            m_c = -m_a * m_barycentre.x() - m_b * m_barycentre.y();
            }

            void compute_reference_coordinates() {
                m_reference_coordinates = m_barycentre;
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_REGULAR_SEGMENT_H