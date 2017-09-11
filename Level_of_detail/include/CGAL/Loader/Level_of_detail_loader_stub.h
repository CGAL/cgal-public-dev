#ifndef CGAL_LEVEL_OF_DETAIL_LOADER_STUB_H
#define CGAL_LEVEL_OF_DETAIL_LOADER_STUB_H

// CGAL includes.
#include <CGAL/array.h>

// Boost includes.
#include <boost/tuple/tuple.hpp>

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
			
			typedef int Label;
			typedef int Types;
			typedef int Index;

			// typedef typename Traits::Plane_3 Plane;

			typedef typename Container:: template Property_map<Color> Color_map;
			typedef typename Container:: template Property_map<Label> Label_map;
			typedef typename Container:: template Property_map<Types> Types_map;
			typedef typename Container:: template Property_map<Index> Index_map;

			// typedef typename Container:: template Property_map<Plane> Plane_map;

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

				Point  point; // , zero(0, 0, 0);
				// Plane  plane, default_plane(zero, zero, zero);
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
				Color_map colors; Label_map labels; Types_map types; Index_map indices; // Plane_map planes;

				boost::tie(colors, success)  = input. template add_property_map<Color>("color", black);
				assert(success);

				boost::tie(labels, success)  = input. template add_property_map<Label>("label", -1);
				assert(success);

				boost::tie(types, success)   = input. template add_property_map<Types>("types", -1);
				assert(success);

				boost::tie(indices, success) = input. template add_property_map<Index>("index", -1);
				assert(success);

				// boost::tie(planes, success)  = input. template add_property_map<Plane>("plane", default_plane);
				// assert(success);

				// All normals are not oriented!
				// All normals are not normalized!
				// It has 38 points and 11 planes.

				// The first building - xy aligned.

				// (Facade 1) Add first vertical facade:
				normal = Normal(0.0, 0.22, 0.0); // plane = Plane(0.0, 0.22, 0.0, 0.0);
				
				point = Point(0.2, 0.0, 0.3); Iter it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 0;
				point = Point(0.4, 0.0, 0.6);      it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 0;
				point = Point(0.8, 0.0, 0.1);      it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 0;

				// (Facade 2) Add second vertical facade || to the first one:
				normal = Normal(0.0, -0.2, 0.0); // plane = Plane(0.0, -0.2, 0.0, 0.2);
				
				point = Point(0.3, 1.0, 0.4); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 1;
				point = Point(0.7, 1.0, 0.8); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 1;
				point = Point(0.5, 1.0, 0.1); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 1;

				// (Facade 3) Add third vertical facade:
				normal = Normal(-0.16, 0.0, 0.0); // plane = Plane(-0.16, 0.0, 0.0, 0.0); 
				  
				point = Point(0.0, 0.1, 0.9); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 2;
				point = Point(0.0, 0.6, 0.6); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 2;
				point = Point(0.0, 0.9, 0.1); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 2;

				// (Facade 4) Add fourth nearly vertical facade (0.05 away from vertical) opposite to the third one:
				normal = Normal(0.36, 0.0, 0.03); // plane = Plane(0.36, 0.0, 0.03, 0.387);
				  
				point = Point(1.0 , 0.2, 0.9); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 3;
				point = Point(1.05, 0.5, 0.3); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 3;
				point = Point(1.0 , 0.8, 0.9); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 3;

				// (Roof) Add a horizontal roof above the four facades defined before.
				normal = Normal(0.0, 0.0, 0.06); // plane = Plane(0.0, 0.0, 0.06, -0.0594);
				  
				point = Point(0.1, 0.1, 0.99); it = input.insert(point, normal); colors[*it] = red; labels[*it] = roof; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 4;
				point = Point(0.3, 0.4, 0.99); it = input.insert(point, normal); colors[*it] = red; labels[*it] = roof; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 4;
				point = Point(0.7, 0.7, 0.99); it = input.insert(point, normal); colors[*it] = red; labels[*it] = roof; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 4;

				// (Ground) Add the ground below the building above:
				normal = Normal(0.0, 0.0, 2.56); // plane = Plane(0.0, 0.0, 2.56, 0.0);
				
				point = Point(-0.1, -0.1, 0.0); it = input.insert(point, normal); colors[*it] = black; labels[*it] = ground; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 5;
				point = Point( 1.5, -0.1, 0.0); it = input.insert(point, normal); colors[*it] = black; labels[*it] = ground; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 5;
				point = Point( 0.5,  1.5, 0.0); it = input.insert(point, normal); colors[*it] = black; labels[*it] = ground; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 5;

				// (Vegetation) Add vegetation in the top right corner of the ground above:
				point = Point(1.40, 1.40,  0.0); normal = Normal( 0.004,  0.0,   1.0);
				it = input.insert(point, normal); colors[*it] = green; labels[*it] = vegetation; types[*it] = -1; /* planes[*it] = default_plane; */ indices[*it] = -1;

				point = Point(1.45, 1.38, 0.10); normal = Normal(-0.046,  0.02,  0.9);
				it = input.insert(point, normal); colors[*it] = green; labels[*it] = vegetation; types[*it] = -1; /* planes[*it] = default_plane; */ indices[*it] = -1;

				point = Point(1.37, 1.42, 0.20); normal = Normal( 0.034, -0.02,  0.8);
				it = input.insert(point, normal); colors[*it] = green; labels[*it] = vegetation; types[*it] = -1; /* planes[*it] = default_plane; */ indices[*it] = -1;
				
				point = Point(1.43, 1.43, 0.05); normal = Normal(-0.026, -0.03, 0.95);
				it = input.insert(point, normal); colors[*it] = green; labels[*it] = vegetation; types[*it] = -1; /* planes[*it] = default_plane; */ indices[*it] = -1;

				point = Point(1.37, 1.37, 0.15); normal = Normal( 0.034,  0.03, 0.85);
				it = input.insert(point, normal); colors[*it] = green; labels[*it] = vegetation; types[*it] = -1; /* planes[*it] = default_plane; */ indices[*it] = -1;

				// Rotated building - not xy aligned.

				// (Facade 1) :
				normal = Normal(0.202854, 0.12171, -3.984e-7); // plane = Plane(-509172.0, -305497.0, 1.0, -631371.0);
				
				point = Point(-0.95172, -0.48047, 0.2); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 6;
				point = Point(-0.89196, -0.58007, 0.8); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 6;
				point = Point(-0.73891, -0.83516, 0.3); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 6;

				// (Facade 2) Follows the first one and includes the exact corners of the building:
				normal = Normal(-0.056698, 0.099226, -3.0e-6); // plane = Plane(18899.3, -33075.3, 1.0, -16538.4);
				
				point = Point(-0.7    , -0.9    , 0.1); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 7;
				point = Point(-0.39871, -0.72783, 0.5); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 7;
				point = Point(0.0     , -0.5    , 0.7); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 7;

				// (Facade 3) Follows the second one:
				normal = Normal(0.145876, 0.175048, -9.419e-7); // plane = Plane(-154874.0, -185846.0, 1.0, -92923.7);
				
				point = Point(-0.05418, -0.45485, 0.7); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 8;
				point = Point(-0.30487, -0.24594, 0.3); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 8;
				point = Point(-0.4918 , -0.09016, 0.7); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 8;

				// (Facade 4) Follows the third one, must be rejected by the regularizer:

				// normal = Normal(0.26, -0.26, -0.03016); plane = Plane(-8.62069, 8.62069, 1.0, -5.40241); // does not reject
				// point = Point(-0.76504, -0.24044, 0.58); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; planes[*it] = plane; indices[*it] = 9; // does not reject, middle point

				normal = Normal(0.26, -0.26, -0.099744); // plane = Plane(-2.60667, 2.60667, 1.0, -1.794);
				
				point = Point(-0.6    , 0.0     , 0.23); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 9;
				point = Point(-0.64175, -0.29111, 0.88); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 9;
				point = Point(-1.0    , -0.4    , 0.23); it = input.insert(point, normal); colors[*it] = blue; labels[*it] = facade; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 9;

				// (Roof) This one is not horizontal :
				normal = Normal(0.053372, 0.052358, 0.120235); // plane = Plane(0.443898, 0.435464, 1.0, -0.683335);
				  
				point = Point(-0.79882, -0.37218, 1.2); it = input.insert(point, normal); colors[*it] = red; labels[*it] = roof; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 10;
				point = Point(-0.53703, -0.63904, 1.2); it = input.insert(point, normal); colors[*it] = red; labels[*it] = roof; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 10;
				point = Point(-0.45596, -0.26240, 1.0); it = input.insert(point, normal); colors[*it] = red; labels[*it] = roof; types[*it] = 0; /* planes[*it] = plane; */ indices[*it] = 10;
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