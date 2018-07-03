#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_CONSISTENCY_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_CONSISTENCY_H

// STL includes.
#include <vector>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        template<class InputTriangulation>
		class Visibility_consistency {

		public:
            using Triangulation = InputTriangulation;

            using Visibility_label  = LOD::Visibility_label;
            using Visibility_labels = std::vector<Visibility_label>;

            using Triangulation_face_handle    = typename Triangulation::Face_handle;
            using Triangulation_faces_iterator = typename Triangulation::Finite_faces_iterator;

            void make_consistent(Triangulation &triangulation) const {

                const size_t num_faces = triangulation.number_of_faces();
				Visibility_labels visibility_labels(num_faces);

				size_t i = 0;
				for (Triangulation_faces_iterator tf_it = triangulation.finite_faces_begin(); tf_it != triangulation.finite_faces_end(); ++tf_it, ++i) {
					const Visibility_label label = tf_it->info().visibility_label();

					const Triangulation_face_handle &fh1 = tf_it->neighbor(0);
					const Triangulation_face_handle &fh2 = tf_it->neighbor(1);
					const Triangulation_face_handle &fh3 = tf_it->neighbor(2);

					const Visibility_label label1 = fh1->info().visibility_label();
					const Visibility_label label2 = fh2->info().visibility_label();
					const Visibility_label label3 = fh3->info().visibility_label();

					if ((label == Visibility_label::INSIDE  && label1 == Visibility_label::OUTSIDE && label2 == Visibility_label::OUTSIDE && label3 == Visibility_label::OUTSIDE) ||
						(label == Visibility_label::OUTSIDE && label1 == Visibility_label::INSIDE  && label2 == Visibility_label::INSIDE  && label3 == Visibility_label::INSIDE) ){
						
						visibility_labels[i] = label1;
						continue;
					}
					visibility_labels[i] = label;
				}

                i = 0;
				for (Triangulation_faces_iterator tf_it = triangulation.finite_faces_begin(); tf_it != triangulation.finite_faces_end(); ++tf_it, ++i)
					tf_it->info().visibility_label() = visibility_labels[i];
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_CONSISTENCY_H