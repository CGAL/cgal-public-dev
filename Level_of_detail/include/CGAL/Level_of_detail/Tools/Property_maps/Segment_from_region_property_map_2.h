#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_FROM_REGION_PROPERTY_MAP_2_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_FROM_REGION_PROPERTY_MAP_2_H

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Fitting/Segment_to_points_fitter.h>

namespace CGAL {

	namespace Level_of_detail {

		namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel, class Elements, class PointMap, class Output>
		class Segment_from_region_property_map_2 {

		public:
			using Kernel 	= InputKernel;
			using Range     = Elements;
			using Point_map = PointMap;
			using Ranges    = std::list<Range>;

			using Line_2 	= typename Kernel::Line_2;
			using Point_2 	= typename Kernel::Point_2;
			using Segment_2 = typename Kernel::Segment_2;

			using KeyType 	= Range;
			using ValueType = Segment_2;

			using Segment_to_points_fitter = LOD::Segment_to_points_fitter<Kernel>;

			Segment_from_region_property_map_2(const Ranges &input, const Point_map &point_map) :
			m_input(input),
			m_point_map(point_map) 
			{ }

			inline const Point_map& point_map() const {
				return m_point_map;
			}

            using key_type   = KeyType;
            using value_type = ValueType;
            using reference  = const value_type&;
			using category   = boost::lvalue_property_map_tag;
            
            using Self = Segment_from_region_property_map_2<Kernel, Elements, Point_map, Output>;

			value_type operator[](key_type &key) const { 
				return get(this, key);
			}

            friend value_type get(const Self &self, const key_type &key) {
				
				ValueType segment;
				self.create_segment_2(key, segment);
				return segment;
            }

			friend void put(const Self&, Output &output, const value_type &value) {
				output.push_back(value);
            }

			inline void create_segment_2(const Range &points, Segment_2 &segment) const {
				
				const Segment_to_points_fitter segment_to_points_fitter;
				segment_to_points_fitter.fit_segment_2(points, m_point_map, segment);
			}

		private:
			const Ranges 	&m_input;
			const Point_map &m_point_map;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_FROM_REGION_PROPERTY_MAP_2_H