#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <memory>
#include <utility>

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Segment_regularizer_parameters.h>

#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Angles_regularizer.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Ordinates_regularizer.h>

#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Regular_segment.h>
#include <CGAL/Level_of_detail/Regularization/Segment_regularizer_2/Regular_segment_property_map.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

        template<class InputKernel>
        class Segment_regularizer_2 {

        public:
            using Kernel = InputKernel;

            using FT      = typename Kernel::FT;
            using Point   = typename Kernel::Point_2;
            using Segment = typename Kernel::Segment_2;
            using Line    = typename Kernel::Line_2;

            using Regular_segment     = LOD::Regular_segment<Kernel>;
            using Regular_segments    = std::vector<Regular_segment *>;
            using In_regular_segments = std::vector<Regular_segment>;
            
            using Angles_regularizer    = LOD::Angles_regularizer<Kernel>;
            using Ordinates_regularizer = LOD::Ordinates_regularizer<Kernel>;

            using Parameters = LOD::Segment_regularizer_parameters<FT>;
            using Tree       = typename Angles_regularizer::Tree;

            Segment_regularizer_2(const Parameters &parameters) : 
            m_parameters(parameters)
            { }

            ~Segment_regularizer_2() {
                m_input_segments.clear();
            }

            template<class Elements, class Segment_map, class Output>
            void regularize(const Elements &elements, const Segment_map &segment_map, Output &output) {
                if (elements.size() == 0) return;

                // Copy input segments in the internal storage.
                create_input(elements, segment_map);
                
                // Make segments parallel and orthogonal.
                if (m_parameters.optimize_angles()) 
                    apply_angles_regularization();

                // Make segments collinear.
                if (m_parameters.optimize_angles() && m_parameters.optimize_ordinates()) 
                    apply_ordinates_regularization();

                // Return all resulting segments.
                create_output(segment_map, output);
            }

        private:            
            const Parameters &m_parameters;
            
            Regular_segments    m_input_segments;
            In_regular_segments m_in_segments;

            std::shared_ptr<Angles_regularizer>    m_angles_regularizer_ptr;
            std::shared_ptr<Ordinates_regularizer> m_ordinates_regularizer_ptr;
            
            template<class Elements, class Segment_map>
            void create_input(const Elements &elements, const Segment_map &segment_map) {
                
                using Const_elements_iterator = typename Elements::const_iterator;
                CGAL_precondition(elements.size() > 0);

                m_in_segments.clear();
                m_in_segments.resize(elements.size());

                size_t i = 0;
                for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it, ++i) {
                    
                    const Segment &segment = get(segment_map, *ce_it);
                    m_in_segments[i] = Regular_segment(i, segment);
                }

                CGAL_precondition(i == m_in_segments.size());

                m_input_segments.clear();
                m_input_segments.resize(m_in_segments.size());

                for (size_t i = 0; i < m_in_segments.size(); ++i)
                    m_input_segments[i] = &m_in_segments[i];
            }

            void apply_angles_regularization() {

                m_angles_regularizer_ptr = std::make_shared<Angles_regularizer>(m_input_segments, m_parameters);
                m_angles_regularizer_ptr->regularize();
            }

            void apply_ordinates_regularization() {
                Tree *tree_pointer = m_angles_regularizer_ptr->get_tree_pointer();

                m_ordinates_regularizer_ptr = std::make_shared<Ordinates_regularizer>(m_input_segments, tree_pointer, m_parameters);
                m_ordinates_regularizer_ptr->regularize();
            }

            template<class Segment_map, class Output>
            void create_output(const Segment_map &segment_map, Output &output) {

                for (size_t i = 0; i < m_input_segments.size(); ++i)
                    put(segment_map, output, m_input_segments[i]->get());
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_2_H