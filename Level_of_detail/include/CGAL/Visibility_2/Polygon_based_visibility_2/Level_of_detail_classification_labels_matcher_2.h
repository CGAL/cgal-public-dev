#ifndef CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_LABELS_MATCHER_2_H
#define CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_LABELS_MATCHER_2_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <utility>
#include <cassert>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class VisibilityOutput>
		class Level_of_detail_classification_labels_matcher_2 {

        public:
            typedef KernelTraits     Kernel;
            typedef VisibilityOutput Visibility_output;

            using FT          = typename Kernel::FT;
            using Point_label = int;
            using Visibility  = Visibility_output;

            Level_of_detail_classification_labels_matcher_2() { }

            void add_visibility(const size_t container_index, const Point_label point_label, Visibility &visibility) const {

				const Point_label ground 	 = 0;
				const Point_label facade 	 = 1;
				const Point_label roof 		 = 2;
				const Point_label vegetation = 3;

				switch (point_label) {

					case ground:
						set_outside(container_index, visibility);
						break;

					case facade:
						set_alike(container_index, visibility);
						break;

					case roof:
						set_inside(container_index, visibility);
						break;

					case vegetation:
						set_outside(container_index, visibility);
						break;

					default:
                        set_unknown(container_index, visibility);
                        break;
				}
			}

			FT match_label(const Point_label point_label) {
				
				const Point_label ground 	 = 0;
				const Point_label facade 	 = 1;
				const Point_label roof 		 = 2;
				const Point_label vegetation = 3;

				switch (point_label) {

					case ground:
						return FT(0);
						
					case facade:
						return FT(1) / FT(2);

					case roof:
						return FT(1);

					case vegetation:
						return FT(0);

					default:
                        return FT(1) / FT(2);
				}
			}

        private:
            inline void set_inside(const size_t container_index, Visibility &visibility) const {
				visibility[container_index].first  += FT(1);
			}

			inline void set_outside(const size_t container_index, Visibility &visibility) const {
				visibility[container_index].second += FT(1);
			}

            void set_alike(const size_t container_index, Visibility &visibility) const {
                visibility[container_index].first  += FT(1) / FT(2);
                visibility[container_index].second += FT(1) / FT(2);
            }

			void set_unknown(const size_t container_index, Visibility &visibility) const {
				visibility[container_index].first  += FT(1) / FT(2);
                visibility[container_index].second += FT(1) / FT(2);
			}
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_LABELS_MATCHER_2_H