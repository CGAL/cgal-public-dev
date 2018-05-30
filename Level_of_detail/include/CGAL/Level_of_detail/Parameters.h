#ifndef CGAL_LEVEL_OF_DETAIL_PARAMETERS_H
#define CGAL_LEVEL_OF_DETAIL_PARAMETERS_H

// STL include.
#include <string>

namespace CGAL {

	namespace Level_of_detail {

		template<typename FT>
		struct Parameters {

		public:
			Parameters() :
			m_path_to_input("default_path"),
			m_verbose(true),
			m_scale(FT(5)), 
			m_alpha(m_scale)
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

			inline FT& scale() {
				return m_scale;
			}

			inline const FT& scale() const {
				return m_scale;
			}

			inline FT& alpha() {
				return m_alpha;
			}

			inline const FT& alpha() const {
				return m_alpha;
			}

			void update_scale_dependent() {
				m_alpha = m_scale;
			}

		private:
			std::string m_path_to_input;
			bool 		m_verbose;

			FT m_scale;
			FT m_alpha;
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PARAMETERS_H