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
			typedef typename Traits::Vector_3 Vector;
			typedef OutputContainer 		  Container;

			typedef CGAL::cpp11::array<unsigned char, 3> Color;
			typedef int Label;
			typedef typename Traits::Plane_3 Plane;

			typedef typename Container:: template Property_map<Color> Color_map;
			typedef typename Container:: template Property_map<Label> Label_map;
			typedef typename Container:: template Property_map<Plane> Plane_map;

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

				Point point0, point1, point2, point;
				Point z(0,0,0);
				Vector normal;

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
				Color_map color; Label_map label; Plane_map plane;

				boost::tie(color, success) = input. template add_property_map<Color>("color", black);
				assert(success);

				boost::tie(label, success) = input. template add_property_map<Label>("label", ground);
				assert(success);

				boost::tie(plane, success) = input. template add_property_map<Plane>("plane", Plane(z, z, z));
				assert(success);

				// All normals are not oriented!

				// (Facade 1) Add first vertical facade:
				normal = Vector(0.0, 0.22, 0.0);
				point0 =   Point(0.2, 0.0, 0.3); point1 = Point(0.4, 0.0, 0.6); point2 = Point(0.8, 0.0, 0.1);

				typename Container::iterator 
				it = input.insert(point0, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point1, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point2, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);

				// (Facade 2) Add second vertical facade || to the first one:
				normal = Vector(0.0, -0.2, 0.0);
				point0 =   Point(0.3, 1.0, 0.4); point1 = Point(0.7, 1.0, 0.8); point2 = Point(0.5, 1.0, 0.1);

				it = input.insert(point0, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point1, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point2, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);

				// (Facade 3) Add third vertical facade:
				normal = Vector(-0.16, 0.0, 0.0);
				point0 =    Point(0.0, 0.1, 0.9); point1 = Point(0.0, 0.6, 0.6); point2 = Point(0.0, 0.9, 0.1); 

				it = input.insert(point0, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point1, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point2, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);

				// (Facade 4) Add fourth nearly vertical facade (0.05 away from vertical) opposite to the third one:
				normal = Vector(0.36, 0.0, 0.03);
				point0 =   Point(1.0 , 0.2, 0.9); point1 = Point(1.05, 0.5, 0.3); point2 = Point(1.0 , 0.8, 0.9); 

				it = input.insert(point0, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point1, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point2, normal); color[*it] = blue; label[*it] = facade; plane[*it] = Plane(point0, point1, point2);

				// (Roof) Add a roof above the four facades defined before.
				normal = Vector(0.0, 0.0, 0.06);
				point0 =  Point(0.1, 0.1, 0.99); point1 = Point(0.3, 0.4, 0.99); point2 = Point(0.7, 0.7, 0.99); 

				it = input.insert(point0, normal); color[*it] = red; label[*it] = roof; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point1, normal); color[*it] = red; label[*it] = roof; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point2, normal); color[*it] = red; label[*it] = roof; plane[*it] = Plane(point0, point1, point2);

				// (Ground) Add the ground below the building above:
				normal = Vector(0.0, 0.0, -2.56);
				point0 =  Point(-0.1, -0.1, 0.0); point1 = Point( 1.5, -0.1, 0.0); point2 = Point( 0.5,  1.5, 0.0); 

				it = input.insert(point0, normal); color[*it] = black; label[*it] = ground; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point1, normal); color[*it] = black; label[*it] = ground; plane[*it] = Plane(point0, point1, point2);
				it = input.insert(point2, normal); color[*it] = black; label[*it] = ground; plane[*it] = Plane(point0, point1, point2);

				// (Vegetation) Add vegetation in the top right corner of the ground above:
				point = Point(1.40, 1.40,  0.0); normal = Vector( 0.004,  0.0,   1.0);
				it = input.insert(point, normal); color[*it] = green; label[*it] = vegetation; plane[*it] = Plane(z, z, z);

				point = Point(1.45, 1.38, 0.10); normal = Vector(-0.046,  0.02,  0.9);
				it = input.insert(point, normal); color[*it] = green; label[*it] = vegetation; plane[*it] = Plane(z, z, z);

				point = Point(1.37, 1.42, 0.20); normal = Vector( 0.034, -0.02,  0.8);
				it = input.insert(point, normal); color[*it] = green; label[*it] = vegetation; plane[*it] = Plane(z, z, z);
				
				point = Point(1.43, 1.43, 0.05); normal = Vector(-0.026, -0.03, 0.95);
				it = input.insert(point, normal); color[*it] = green; label[*it] = vegetation; plane[*it] = Plane(z, z, z);

				point = Point(1.37, 1.37, 0.15); normal = Vector( 0.034,  0.03, 0.85);
				it = input.insert(point, normal); color[*it] = green; label[*it] = vegetation; plane[*it] = Plane(z, z, z);
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