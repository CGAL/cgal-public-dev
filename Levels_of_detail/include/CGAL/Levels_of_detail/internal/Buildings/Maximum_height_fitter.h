#ifndef CGAL_LEVEL_OF_DETAIL_MAXIMUM_HEIGHT_FITTER_H
#define CGAL_LEVEL_OF_DETAIL_MAXIMUM_HEIGHT_FITTER_H

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel>
		class Maximum_height_fitter {

		public:
			using Kernel = InputKernel;
            using FT     = typename Kernel::FT;

			Maximum_height_fitter() : 
            m_max_value(FT(0))
            { }

			void add_height(const FT value) {
        m_max_value = (std::max) (value, m_max_value);
			}

			FT get_result() const {
				return m_max_value;
			}

			void clear() {
				m_max_value = FT(0);
			}

		private:
			FT m_max_value;
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_MAXIMUM_HEIGHT_FITTER_H
