#ifndef CGAL_LEVEL_OF_DETAIL_SCORES_BASED_SORTING_H
#define CGAL_LEVEL_OF_DETAIL_SCORES_BASED_SORTING_H

namespace CGAL {

	namespace Level_of_detail {

		template<class ScoresContainer>
		class Scores_based_sorting {

        public:
            using Scores_container   = ScoresContainer;
            using Element_identifier = typename Scores_container::Element_identifier;

            Scores_based_sorting(const Scores_container &scores_container) :
            m_scores_container(scores_container) 
            { }

            bool operator() (const Element_identifier &i, const Element_identifier &j) const {
			    return m_scores_container.scores().at(i) > m_scores_container.scores().at(j);
			}

        private:
            const Scores_container &m_scores_container;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_SCORES_BASED_SORTING_H