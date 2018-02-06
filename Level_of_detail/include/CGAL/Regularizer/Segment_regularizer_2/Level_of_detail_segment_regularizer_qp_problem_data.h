#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_QP_PROBLEM_DATA_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_QP_PROBLEM_DATA_H

// Eigen includes.
#include <eigen3/Eigen/SparseCore>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits, class NeighboursGraphData>
		class Level_of_detail_segment_regularizer_qp_problem_data {

        public:
            typedef KernelTraits        Kernel;
            typedef NeighboursGraphData Neighbours_graph_data;

            using Mus       = typename Neighbours_graph_data::Mus;
            using Targets   = typename Neighbours_graph_data::Targets;
            using Relations = typename Neighbours_graph_data::Relations;

            using FT = typename Kernel::FT;
            
            using Mus_matrix       = Eigen::SparseMatrix<FT>;
            using Targets_matrix   = Eigen::SparseMatrix<FT>;
            using Relations_matrix = Eigen::SparseMatrix<int>;

            Level_of_detail_segment_regularizer_qp_problem_data() : m_mus_matrix(), m_targets_matrix(), m_relations_matrix() { }

            inline Mus_matrix &get_mus_matrix() {
                return m_mus_matrix;
            }

            inline const Mus_matrix &get_mus_matrix() const {
                return m_mus_matrix;
            }

            inline Targets_matrix &get_targets_matrix() {
                return m_targets_matrix;
            }

            inline const Targets_matrix &get_targets_matrix() const {
                return m_targets_matrix;
            }

            inline Relations_matrix &get_relations_matrix() {
                return m_relations_matrix;
            }

            inline const Relations_matrix &get_relations_matrix() const {
                return m_relations_matrix;
            }

            void clear() {

                m_mus_matrix.resize(0, 0);
                m_mus_matrix.data().squeeze();

                m_targets_matrix.resize(0, 0);
                m_targets_matrix.data().squeeze();

                m_relations_matrix.resize(0, 0);
                m_relations_matrix.data().squeeze();
            }

            inline bool filled() const {
                return 
                (m_mus_matrix.innerSize()       != 0 || m_mus_matrix.outerSize()       != 0) &&
                (m_targets_matrix.innerSize()   != 0 || m_targets_matrix.outerSize()   != 0) &&
                (m_relations_matrix.innerSize() != 0 || m_relations_matrix.outerSize() != 0);
            }

            void set_from(const Neighbours_graph_data &data, const size_t global_size) {
                
                clear();
                
                set_mus_matrix_from(data, global_size);
                set_targets_matrix_from(data, global_size);
                set_relations_matrix_from(data, global_size);
            }

        private:
            Mus_matrix       m_mus_matrix;
            Targets_matrix   m_targets_matrix;
            Relations_matrix m_relations_matrix;

            void set_mus_matrix_from(const Neighbours_graph_data &data, const size_t global_size) {
                const Mus &mus = data.get_mus();

                m_mus_matrix.resize(global_size, global_size);
                m_mus_matrix.setFromTriplets(mus.begin(), mus.end());

                m_mus_matrix.makeCompressed();
            }

            void set_targets_matrix_from(const Neighbours_graph_data &data, const size_t global_size) {
                const Targets &targets = data.get_targets();

                m_targets_matrix.resize(global_size, global_size);
                m_targets_matrix.setFromTriplets(targets.begin(), targets.end());

                m_targets_matrix.makeCompressed();
            }

            void set_relations_matrix_from(const Neighbours_graph_data &data, const size_t global_size) {
                const Relations &relations = data.get_relations();

                m_relations_matrix.resize(global_size, global_size);
                m_relations_matrix.setFromTriplets(relations.begin(), relations.end());
                
                m_relations_matrix.makeCompressed();
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_QP_PROBLEM_DATA_H