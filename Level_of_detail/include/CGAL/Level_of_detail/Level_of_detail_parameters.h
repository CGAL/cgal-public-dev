#ifndef CGAL_LEVEL_OF_DETAIL_PARAMETERS_H
#define CGAL_LEVEL_OF_DETAIL_PARAMETERS_H

// STL include.
#include <string>

namespace CGAL {

	namespace Level_of_detail {

		template<typename FT>
		struct Level_of_detail_parameters {

		public:
			Level_of_detail_parameters() :
			m_path_to_input("default_path"),
			m_verbose(true)
			{ }

			//////////////////////////////////
			// Functions to be not documented:

			inline std::string& path_to_input() {
				return m_path_to_input;
			}

			inline const std::string& path_to_input() const {
				return m_path_to_input;
			}

			//////////////////////////////////
			// Functions to be documented:

			inline bool& verbose() {
				return m_verbose;
			}

			inline const bool& verbose() const {
				return m_verbose;
			}

			// to be added!

		private:
			std::string m_path_to_input;
			bool 		m_verbose;

			// add here all LOD parameters that can be provided by the user!
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PARAMETERS_H