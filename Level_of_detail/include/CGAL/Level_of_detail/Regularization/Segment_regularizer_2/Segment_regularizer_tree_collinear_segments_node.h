#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_COLLINEAR_SEGMENTS_NODE_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_COLLINEAR_SEGMENTS_NODE_H

// STL includes.
#include <map>
#include <list>
#include <memory>

// LOD includes.
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Regular_segment.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        template<class InputKernel>
		class Segment_regularizer_tree_collinear_segments_node {

        public:
            using Kernel = InputKernel;
            using FT     = typename Kernel::FT;

            using Regular_segment    = LOD::Regular_segment<Kernel>;
            using Collinear_segments = std::list<Regular_segment *>;

            Segment_regularizer_tree_collinear_segments_node() { 
                allocate_memory();
            }

            ~Segment_regularizer_tree_collinear_segments_node() { 
                deallocate_memory();
            }

            inline const Collinear_segments &get_collinear_segments() const {
                return m_collinear_segments;
            }

            inline Collinear_segments &get_collinear_segments() {
                return m_collinear_segments;
            }

            inline void add(Regular_segment *segment_pointer) {
                m_collinear_segments.push_back(segment_pointer);
            }

        private:
            Collinear_segments m_collinear_segments;

            inline void allocate_memory() {
                m_collinear_segments = Collinear_segments();
            }

            inline void deallocate_memory() {
                m_collinear_segments.clear();
            }
		};

	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_COLLINEAR_SEGMENTS_NODE_H