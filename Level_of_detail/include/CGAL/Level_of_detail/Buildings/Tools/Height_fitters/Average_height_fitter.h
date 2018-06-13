#ifndef CGAL_LEVEL_OF_DETAIL_AVERAGE_HEIGHT_FITTER_H
#define CGAL_LEVEL_OF_DETAIL_AVERAGE_HEIGHT_FITTER_H

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel>
		class Average_height_fitter {

		public:
			using Kernel = InputKernel;
            using FT     = typename Kernel::FT;

			Average_height_fitter() : 
            m_num_values(FT(0)), 
            m_sum_height(FT(0)) 
            { }

			void add_height(const FT value) {
				CGAL_precondition(value >= FT(0));

				m_sum_height += value;
				m_num_values += FT(1);
			}

			FT get_result() const {

				CGAL_precondition(m_num_values != FT(0));
				return m_sum_height / m_num_values;
			}

			void clear() {
				m_num_values = FT(0);
				m_sum_height = FT(0);
			}

		private:
			FT m_num_values;
			FT m_sum_height;
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_AVERAGE_HEIGHT_FITTER_H