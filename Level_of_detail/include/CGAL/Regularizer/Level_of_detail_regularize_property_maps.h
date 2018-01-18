#ifndef CGAL_LEVEL_OF_DETAIL_REGULARIZE_PROPERTY_MAPS_H
#define CGAL_LEVEL_OF_DETAIL_REGULARIZE_PROPERTY_MAPS_H

namespace CGAL {

	namespace LOD {

		template <typename Kernel>
		class Point_to_shape_index_map {
			
		private:
			typedef typename Kernel::Point_3 Point_3;
			typedef typename Kernel::Plane_3 Shape;

			boost::shared_ptr< std::vector<int> > m_indices;

		public:
			typedef std::size_t key_type;
			typedef int value_type;
			typedef value_type reference;
			typedef boost::readable_property_map_tag category;

			Point_to_shape_index_map(const std::vector< std::pair<Point_3, int> > &points) : m_indices(new std::vector<int>(points.size(), -1)) {
				for (std::size_t i = 0; i < points.size(); ++i) (*m_indices)[i] = points[i].second;
			}

			inline friend value_type get(const Point_to_shape_index_map &pm, const key_type &k) {
				return (*(pm.m_indices))[k];
			}
		};

	} // namespace LOD

} // namespace CGAL

#endif // CGAL_LEVEL_OF_DETAIL_REGULARIZE_PROPERTY_MAPS_H