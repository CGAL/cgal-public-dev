#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_DELAUNAY_NEIGHBOURS_GRAPH_BUILDER_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_DELAUNAY_NEIGHBOURS_GRAPH_BUILDER_H

// STL includes.
#include <map>
#include <list>
#include <cassert>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_neighbours_graph_data.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_delaunay_neighbours_graph_builder {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;

            using Regular_segment  = CGAL::LOD::Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Regular_segments = std::vector<Regular_segment>;
            using Orientations     = std::vector<FT>;
            
            using Parameters            = CGAL::LOD::Level_of_detail_segment_regularizer_parameters<Kernel>;
            using Neighbours_graph_data = CGAL::LOD::Level_of_detail_segment_regularizer_neighbours_graph_data<Kernel>;
            
            Level_of_detail_segment_regularizer_delaunay_neighbours_graph_builder(
                const Regular_segments &segments, 
                const Orientations     &orientations, 
                const Parameters       &parameters) : 
                m_segments(segments), m_orientations(orientations), m_parameters(parameters) { }

            void build_graph(Neighbours_graph_data &graph_data) const {
                graph_data.clear();
            
                // to be added!
            }

        private:
            const Regular_segments &m_segments;
            const Orientations     &m_orientations;
            const Parameters       &m_parameters;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_DELAUNAY_NEIGHBOURS_GRAPH_BUILDER_H