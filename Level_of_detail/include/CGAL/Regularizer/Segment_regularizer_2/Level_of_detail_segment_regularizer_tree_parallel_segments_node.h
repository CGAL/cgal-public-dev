#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_PARALLEL_SEGMENTS_NODE_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_PARALLEL_SEGMENTS_NODE_H

// STL includes.
#include <map>
#include <list>
#include <memory>
#include <cassert>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_tree_collinear_segments_node.h>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_tree_parallel_segments_node {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;

            using Regular_segment   = Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Parallel_segments = std::list<Regular_segment *>;

            using Parallel_segments_iterator       = typename Parallel_segments::iterator;
            using Parallel_segments_const_iterator = typename Parallel_segments::const_iterator;

            using Collinear_segments_tree_node = Level_of_detail_segment_regularizer_tree_collinear_segments_node<Kernel>;
            using Collinear_segments           = std::map<FT, Collinear_segments_tree_node *>;
            
            using Collinear_segments_iterator       = typename Collinear_segments::iterator;
            using Collinear_segments_const_iterator = typename Collinear_segments::const_iterator;

            using List_iterator = typename std::list<Regular_segment *>::iterator;

            Level_of_detail_segment_regularizer_tree_parallel_segments_node() { 
                allocate_memory();
            }

            ~Level_of_detail_segment_regularizer_tree_parallel_segments_node() { 
                deallocate_memory();
            }

            inline const Parallel_segments &get_parallel_segments() const {
                return m_parallel_segments;
            }

            inline Parallel_segments &get_parallel_segments() {
                return m_parallel_segments;
            }

            inline const Collinear_segments &get_collinear_segments() const {
                return m_collinear_segments;
            }

            inline Collinear_segments &get_collinear_segments() {
                return m_collinear_segments;
            }

            void add(Regular_segment *segment_pointer) {
                
                m_parallel_segments.push_back(segment_pointer);
                segment_pointer->parallel_node = this;
            }

            inline void delete_parallel_segments() {
                m_parallel_segments.clear();
            }

            void delete_collinear_segments() {
                for (Collinear_segments_iterator it_m = m_collinear_segments.begin(); it_m != m_collinear_segments.end(); ++it_m)
                    delete it_m->second;
                m_collinear_segments.clear();
            }

            void create_collinear_node(const FT ordinate) {
                if (m_collinear_segments.find(ordinate) == m_collinear_segments.end())
                    m_collinear_segments[ordinate] = new Collinear_segments_tree_node();
            }

            void assign_to_collinear_node(const FT ordinate, Regular_segment* segment_pointer) {
                if (m_collinear_segments.find(ordinate) != m_collinear_segments.end())
                    m_collinear_segments[ordinate]->add(segment_pointer);
            }

        private:
            Parallel_segments  m_parallel_segments;
            Collinear_segments m_collinear_segments;

            void allocate_memory() {
                m_parallel_segments  = Parallel_segments();
                m_collinear_segments = Collinear_segments();
            }

            void deallocate_memory() {
                
                delete_references_to_parallel_node();

                delete_collinear_segments();
                delete_parallel_segments();
            }

            void delete_references_to_parallel_node() {
                
                for (Collinear_segments_iterator it_m = m_collinear_segments.begin(); it_m != m_collinear_segments.end(); ++it_m) {
                    Collinear_segments_tree_node* node = it_m->second;

                    for (List_iterator it_s = node->get_collinear_segments().begin(); it_s != node->get_collinear_segments().end(); ++it_s)
                        (*it_s)->parallel_node = NULL;
                }

                for (List_iterator it_s = m_parallel_segments.begin(); it_s != m_parallel_segments.end(); ++it_s)
                    (*it_s)->parallel_node = NULL;
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_TREE_PARALLEL_SEGMENTS_NODE_H