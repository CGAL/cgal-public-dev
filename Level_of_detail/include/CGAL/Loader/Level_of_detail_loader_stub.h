#ifndef CGAL_LEVEL_OF_DETAIL_LOADER_STUB_H
#define CGAL_LEVEL_OF_DETAIL_LOADER_STUB_H

// CGAL includes.
#include <CGAL/array.h>
#include <CGAL/Point_set_3.h>

// New CGAL includes.
#include <CGAL/Loader/Level_of_detail_loader.h>

namespace CGAL {

	namespace LOD {

		enum class Mock_data_type { BASIC, DEFAULT };

		template<class KernelTraits, class OutputContainer>
		class Level_of_detail_loader_stub : public Level_of_detail_loader<KernelTraits, OutputContainer> {
		
		public:
			typedef KernelTraits 			  Traits;
			typedef typename Traits::Point_3  Point;
			typedef typename Traits::Vector_3 Normal;
			typedef OutputContainer 		  Container;

			typedef unsigned char 				Type;
			typedef CGAL::cpp11::array<Type, 3> Color;
			typedef int 						Label;
			typedef typename Traits::Plane_3 	Plane;
			typedef int 						Index;

			typedef typename Container:: template Property_map<Color> Color_map;
			typedef typename Container:: template Property_map<Label> Label_map;
			typedef typename Container:: template Property_map<Plane> Plane_map;
			typedef typename Container:: template Property_map<Index> Index_map;

			typedef typename Container::iterator Iter;

			Level_of_detail_loader_stub() : m_mock_data_type(Mock_data_type::BASIC) { }

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

			// For points see the geogebra file ggb/basic.ggb.
			// Here I also use double as my primary scalar type.
			void get_mock_basic(Container &input) const {
				
				input.clear();

				Point  point, zero(0, 0, 0);
				Plane  plane, default_plane(zero, zero, zero);
				Normal normal;

				input.add_normal_map();

				const Color red   = {{255, 0, 0}};
				const Color green = {{0, 255, 0}};
				const Color blue  = {{0, 0, 255}};
				const Color black = {{0, 0,   0}};

				const Label ground     = 0;
				const Label facade     = 1;
				const Label roof       = 2;
				const Label vegetation = 3; 

				bool success = false;
				Color_map colors; Label_map labels; Plane_map planes; Index_map indices;

				boost::tie(colors, success)  = input. template add_property_map<Color>("color", black);
				assert(success);

				boost::tie(labels, success)  = input. template add_property_map<Label>("label", vegetation);
				assert(success);

				boost::tie(planes, success)  = input. template add_property_map<Plane>("plane", default_plane);
				assert(success);

				boost::tie(indices, success) = input. template add_property_map<Index>("index", -1);
				assert(success);

				// All normals are not oriented!
				// All normals are not normalized!

				// (Facade 1) Add first vertical facade:
				normal = Normal(0.0, 0.22, 0.0); plane = Plane(0.0, 0.22, 0.0, 0.0);
				
				point = Point(0.2, 0.0, 0.3); Iter it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 0;
				point = Point(0.4, 0.0, 0.6); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 0;
				point = Point(0.8, 0.0, 0.1); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 0;

				// (Facade 2) Add second vertical facade || to the first one:
				normal = Normal(0.0, -0.2, 0.0); plane = Plane(0.0, -0.2, 0.0, 0.2);
				
				point = Point(0.3, 1.0, 0.4); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 1;
				point = Point(0.7, 1.0, 0.8); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 1;
				point = Point(0.5, 1.0, 0.1); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 1;

				// (Facade 3) Add third vertical facade:
				normal = Normal(-0.16, 0.0, 0.0); plane = Plane(-0.16, 0.0, 0.0, 0.0); 
				  
				point = Point(0.0, 0.1, 0.9); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 2;
				point = Point(0.0, 0.6, 0.6); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 2;
				point = Point(0.0, 0.9, 0.1); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 2;

				// (Facade 4) Add fourth nearly vertical facade (0.05 away from vertical) opposite to the third one:
				normal = Normal(0.36, 0.0, 0.03); plane = Plane(0.36, 0.0, 0.03, 0.387);
				  
				point = Point(1.0 , 0.2, 0.9); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 3;
				point = Point(1.05, 0.5, 0.3); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 3;
				point = Point(1.0 , 0.8, 0.9); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 3;

				// (Roof) Add a roof above the four facades defined before.
				normal = Normal(0.0, 0.0, 0.06); plane = Plane(0.0, 0.0, 0.06, -0.0594);
				  
				point = Point(0.1, 0.1, 0.99); it = input.insert(point, normal); colors[*it] = red; labels[*it] = roof; planes[*it] = plane; indices[*it] = 4;
				point = Point(0.3, 0.4, 0.99); it = input.insert(point, normal); colors[*it] = red; labels[*it] = roof; planes[*it] = plane; indices[*it] = 4;
				point = Point(0.7, 0.7, 0.99); it = input.insert(point, normal); colors[*it] = red; labels[*it] = roof; planes[*it] = plane; indices[*it] = 4;

				// (Ground) Add the ground below the building above:
				normal = Normal(0.0, 0.0, 2.56); plane = Plane(0.0, 0.0, 2.56, 0.0);
				
				point = Point(-0.1, -0.1, 0.0); it = input.insert(point, normal); colors[*it] = black; labels[*it] = ground; planes[*it] = plane; indices[*it] = 5;
				point = Point( 1.5, -0.1, 0.0); it = input.insert(point, normal); colors[*it] = black; labels[*it] = ground; planes[*it] = plane; indices[*it] = 5;
				point = Point( 0.5,  1.5, 0.0); it = input.insert(point, normal); colors[*it] = black; labels[*it] = ground; planes[*it] = plane; indices[*it] = 5;

				// (Vegetation) Add vegetation in the top right corner of the ground above:
				point = Point(1.40, 1.40,  0.0); normal = Normal( 0.004,  0.0,   1.0);
				it = input.insert(point, normal); colors[*it] = green; labels[*it] = vegetation; planes[*it] = default_plane; indices[*it] = -1;

				point = Point(1.45, 1.38, 0.10); normal = Normal(-0.046,  0.02,  0.9);
				it = input.insert(point, normal); colors[*it] = green; labels[*it] = vegetation; planes[*it] = default_plane; indices[*it] = -1;

				point = Point(1.37, 1.42, 0.20); normal = Normal( 0.034, -0.02,  0.8);
				it = input.insert(point, normal); colors[*it] = green; labels[*it] = vegetation; planes[*it] = default_plane; indices[*it] = -1;
				
				point = Point(1.43, 1.43, 0.05); normal = Normal(-0.026, -0.03, 0.95);
				it = input.insert(point, normal); colors[*it] = green; labels[*it] = vegetation; planes[*it] = default_plane; indices[*it] = -1;

				point = Point(1.37, 1.37, 0.15); normal = Normal( 0.034,  0.03, 0.85);
				it = input.insert(point, normal); colors[*it] = green; labels[*it] = vegetation; planes[*it] = default_plane; indices[*it] = -1;
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