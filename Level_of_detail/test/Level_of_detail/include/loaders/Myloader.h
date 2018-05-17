#ifndef CGAL_LEVEL_OF_DETAIL_MYLOADER_H
#define CGAL_LEVEL_OF_DETAIL_MYLOADER_H

// STL includes.
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/array.h>

namespace CGAL {

	namespace Level_of_detail {

		template<class KernelTraits, class OutputContainer>
		class Myloader {
		
		public:
			typedef KernelTraits    Traits;
			typedef OutputContainer Container;

			typedef typename Traits::FT       FT;
			typedef typename Traits::Point_3  Point;
			typedef typename Traits::Vector_3 Normal;

			typedef unsigned char 				Type;
			typedef CGAL::cpp11::array<Type, 3> Color;

			typedef int Label;
			typedef int Types;
			typedef int Index;

			typedef typename Container:: template Property_map<Color> Color_map;
			typedef typename Container:: template Property_map<Label> Label_map;
			typedef typename Container:: template Property_map<Types> Types_map;
			typedef typename Container:: template Property_map<Index> Index_map;

			typedef typename Container::iterator Iterator;

			Myloader(Traits traits = Traits())
			: m_traits(traits) { }

			void get_data(const std::string &filePath, Container &input) const {
            	
				// We use this to fix wrong order of the given input data in the file.
            	// std::ofstream saver((filePath + ".new").c_str(), std::ios_base::out); saver.precision(20);

				std::ifstream loader(filePath.c_str(), std::ios_base::in);
            	if (!loader) {

                	std::cerr << std::endl << std::endl << "ERROR: Error loading file with LOD data!" << std::endl << std::endl;
                	exit(EXIT_FAILURE);
            	}

				Color_map colors; Label_map labels; Types_map types; Index_map indices;
				set_default_properties(input, colors, labels, types, indices);

            	std::string tmp;
            	std::getline(loader, tmp);
            	std::getline(loader, tmp);
            	std::getline(loader, tmp);
            	std::getline(loader, tmp);
            	std::getline(loader, tmp);
            	std::getline(loader, tmp);
            	std::getline(loader, tmp);
            	std::getline(loader, tmp);

            	size_t num_points;
            	loader >> tmp >> tmp >> num_points;

            	for (size_t i = 0; i < 13; ++i) std::getline(loader, tmp);

				FT x, y, z, nx, ny, nz;
				int r, g, b, la;

            	for (size_t i = 0; i < num_points; ++i) {
            		loader >> x >> y >> z >> nx >> ny >> nz >> r >> g >> b >> tmp >> la;
            		
					// We use this to fix wrong order of the given input data in the file.
					// saver << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << " " << r << " " << g << " " << b << " " << tmp << " " << la << std::endl;

					Iterator it = input.insert(Point(x, y, z), Normal(nx, ny, nz));

					 colors[*it] = {{static_cast<Type>(r), static_cast<Type>(g), static_cast<Type>(b)}};
					 labels[*it] = la;
					  types[*it] = -1; 
					indices[*it] = -1;
            	}
            	loader.close();
            	
				// We use this to fix wrong order of the given input data in the file.
            	// saver.close(); exit(0);
			}

		private:
			Traits m_traits;

			void set_default_properties(Container &input, Color_map &colors, Label_map &labels, Types_map &types, Index_map &indices) const {

				bool success = false;
				input.add_normal_map();

				boost::tie(colors , success) = input. template add_property_map<Color>("color", {{0, 0, 0}});
				assert(success);

				boost::tie(labels , success) = input. template add_property_map<Label>("label", -1);
				assert(success);

				boost::tie(types  , success) = input. template add_property_map<Types>("types", -1);
				assert(success);

				boost::tie(indices, success) = input. template add_property_map<Index>("index", -1);
				assert(success);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_MYLOADER_H