#ifndef CGAL_LEVEL_OF_DETAIL_FIXED_10_METERS_HEIGHT_FITTER_H
#define CGAL_LEVEL_OF_DETAIL_FIXED_10_METERS_HEIGHT_FITTER_H

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel>
		class Fixed_10_meters_height_fitter {

		public:
			using Kernel = InputKernel;
            using FT     = typename Kernel::FT;

			void add_height(const FT ) { }

			FT get_result() const {
				return FT(10);
			}

			void clear() { }
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_FIXED_10_METERS_HEIGHT_FITTER_H