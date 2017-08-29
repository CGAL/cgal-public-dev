#ifndef CGAL_MYLOG_H
#define CGAL_MYLOG_H

// STL includes.
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/array.h>

namespace CGAL {

	namespace LOD {

		class Mylog {

		public:
			std::string state() const {
				return "ok";
			}

			std::string data() const {
				return out.str();
			}

			void clear() {
				out.str(std::string());
			}

			template<typename T>
			void append(const T &token) {
				out << token;
			}

			bool save(const std::string &fileName, const std::string &extension = ".log", const std::string path = "/Users/danisimo/Documents/pipeline/logs/") const {

				const std::string finalPath = path + fileName + extension;
				std::ofstream file(finalPath.c_str(), std::ios_base::out);

				if (!file) {
					std::cerr << "\nERROR: Error saving log file with the name " << fileName << "\n" << std::endl;
					return false;
				}

				file << data() << std::endl;
				file.close();

				return true;
			}

			template<class Traits, class Container>
			void save_ply(const Container &input, 
						  const std::string &fileName,
						  const bool withExtraProperties = false) {

				using Type  = unsigned char;
				using Color = CGAL::cpp11::array<Type, 3>;
				using Label = int;
				using Plane = typename Traits::Plane_3;
				using Index = int;

				using Color_map = typename Container:: template Property_map<Color>;
				using Label_map = typename Container:: template Property_map<Label>;
				using Plane_map = typename Container:: template Property_map<Plane>;
				using Index_map = typename Container:: template Property_map<Index>;

				typedef typename Container::const_iterator Iter;

				clear();

				Color_map colors;
				Label_map labels;
				Plane_map planes;
				Index_map indices;

				out << 
				"ply\n"                  << 
				"format ascii 1.0\n"     << 
				"element vertex "        << input.number_of_points() << "\n" << 
				"property double x\n"    << 
				"property double y\n"    << 
				"property double z\n"    <<
				"property double nx\n"   <<
				"property double ny\n"   <<
				"property double nz\n"   <<
				"property uchar red\n"   << 
				"property uchar green\n" <<
				"property uchar blue\n";

				boost::tie(colors,  boost::tuples::ignore) = input. template property_map<Color>("color");

				if (withExtraProperties) {

					out << 
					"property int label\n"  <<
					"property int index\n"  <<
					"property double planea\n" <<
					"property double planeb\n" <<
					"property double planec\n" <<
					"property double planed\n" <<
					"end_header\n";
					
					boost::tie(labels,  boost::tuples::ignore) = input. template property_map<Label>("label");
					boost::tie(indices, boost::tuples::ignore) = input. template property_map<Index>("index");
					boost::tie(planes,  boost::tuples::ignore) = input. template property_map<Plane>("plane");
				
				} else out << "end_header\n";

				for (Iter it = input.begin(); it != input.end(); ++it) {

					// if (static_cast<int>(*it) % 3 == 0) out << "\n"; // remove if not needed

					out.precision(10);

					out 
					<< input.point(*it)  << " " 
					<< input.normal(*it) << " "
					<< static_cast<int>(colors[*it][0]) << " " 
					<< static_cast<int>(colors[*it][1]) << " " 
					<< static_cast<int>(colors[*it][2]);

					if (withExtraProperties) out << " " <<  labels[*it] << " " << indices[*it] << " " <<  planes[*it];
					
					out << "\n";
				}
				save(fileName, ".log");
			}

			std::stringstream out;
		};
	}
}

#endif // CGAL_MYLOG_H