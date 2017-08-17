#ifndef CGAL_MYLOG_H
#define CGAL_MYLOG_H

// STL includes.
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

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

			bool save(const std::string &fileName) const {

				const std::string path = "/Users/danisimo/Documents/pipeline/logs/";
				const std::string finalPath = path + fileName + ".log";

				std::ofstream file(finalPath.c_str(), std::ios_base::out);

				if (!file) {
					std::cerr << "\nERROR: Error saving log file with the name " << fileName << "\n" << std::endl;
					return false;
				}

				file << data() << std::endl;
				file.close();

				return true;
			}

			std::stringstream out;
		};
	}
}

#endif // CGAL_MYLOG_H