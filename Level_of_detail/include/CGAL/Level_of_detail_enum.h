#ifndef CGAL_LEVEL_OF_DETAIL_ENUM_H
#define CGAL_LEVEL_OF_DETAIL_ENUM_H

namespace CGAL {

	namespace LOD {

		// Visibility labels.
		enum class Visibility_label { 
			IN,  	// inside a triangle
			OUT, 	// outside a triangle
			UNKNOWN // 50/50 inside outside
		};

		// Visibility methods.
		enum class Visibility_method {
			CLASSIFICATION,
			SAMPLING
		};

		// Visibility approaches.
		enum class Visibility_approach {
			POINT_BASED,
			FACE_BASED
		};

		// Face info class.
		template<typename FT>
		class My_face_info {

		public:
			FT in = FT(1) / FT(2);
		};

		// Visibility samplers.
		enum class Visibility_sampler {
			REGULAR,
			RANDOM
		};

	} // LOD

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ENUM_H