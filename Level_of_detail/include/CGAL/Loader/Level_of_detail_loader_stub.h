#ifndef CGAL_LEVEL_OF_DETAIL_LOADER_STUB_H
#define CGAL_LEVEL_OF_DETAIL_LOADER_STUB_H

// New CGAL includes.
#include <CGAL/Loader/Level_of_detail_loader.h>

namespace CGAL {

	namespace LOD {

		enum class Mock_data_type { BASIC };

		template<class KernelTraits, class OutputContainer>
		class Level_of_detail_loader_stub : public Level_of_detail_loader<KernelTraits, OutputContainer> {
		
		public:
			typedef KernelTraits 			  Traits;
			typedef typename Traits::Point_3  Point;
			typedef typename Traits::Vector_3 Vector;
			typedef OutputContainer 		  Container;

			void get_data(const std::string &, Container &input) const override {

				switch(m_mock_data_type) {
					
					case Mock_data_type::BASIC:
						get_mock_basic(input);
						return;

					default:	
						get_mock_default(input);	
						return;
				}
			}

			void insert_point(const Point &newPoint, Container &input) const {
				input.insert(newPoint);
			}

			void get_mock_basic(Container &input) const {
				
				input.clear();

				const Point   point_0 =  Point(0, 0, 0);
				const Vector normal_0 = Vector(1, 1, 1);

				input.add_normal_map();

				input.insert(point_0, normal_0);
			}

			void get_mock_default(Container &input) const {

				input.clear();

				const Point new_point = Point(0, 0, 0);

				input.insert(new_point);
			}

			inline void set_mock_data_type(const Mock_data_type mock_data_type) {

				m_mock_data_type = mock_data_type;
			}

		private:
			Mock_data_type m_mock_data_type;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_LOADER_STUB_H