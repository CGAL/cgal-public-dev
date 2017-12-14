#ifndef CGAL_LEVEL_OF_DETAIL_DISTORTION_H
#define CGAL_LEVEL_OF_DETAIL_DISTORTION_H

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class LodReconstruction>
		class Level_of_detail_distortion {

		public:
			typedef KernelTraits   	  Kernel;
			typedef InputContainer 	  Container;
			typedef LodReconstruction LODS;

			typedef typename Kernel::FT FT;

			Level_of_detail_distortion(const Container &input, const LODS &lods) : m_input(input), m_lods(lods), m_distortion(-FT(2)) { }

			void estimate() {
				
			}

			FT get() const {
				return m_distortion;
			}

		private:
			const Container &m_input;
			const LODS 		&m_lods;

			FT m_distortion;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_DISTORTION_H