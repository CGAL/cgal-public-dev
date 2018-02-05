#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_QP_SOLVER_DATA_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_QP_SOLVER_DATA_H

// Eigen includes.
#include <eigen3/Eigen/SparseCore>

namespace CGAL {

	namespace LOD {

        template<class KernelTraits>
		class Level_of_detail_segment_regularizer_qp_solver_data {

        public:
            typedef KernelTraits Kernel;
            using FT = typename Kernel::FT;
            
            using Mu_matrix        = Eigen::SparseMatrix<FT>;
            using Targets_matrix   = Eigen::SparseMatrix<FT>;
            using Relations_matrix = Eigen::SparseMatrix<int>;

            Level_of_detail_segment_regularizer_qp_solver_data() : m_mu_matrix(), m_targets_matrix(), m_relations_matrix() { }

            inline Mu_matrix &get_mu_matrix() {
                return m_mu_matrix;
            }

            inline const Mu_matrix &get_mu_matrix() const {
                return m_mu_matrix;
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

                m_mu_matrix.resize(0, 0);
                m_mu_matrix.data().squeeze();

                m_targets_matrix.resize(0, 0);
                m_targets_matrix.data().squeeze();

                m_relations_matrix.resize(0, 0);
                m_relations_matrix.data().squeeze();
            }

            bool filled() {
                return 
                (m_mu_matrix.innerSize()        != 0 || m_mu_matrix.outerSize()        != 0) &&
                (m_targets_matrix.innerSize()   != 0 || m_targets_matrix.outerSize()   != 0) &&
                (m_relations_matrix.innerSize() != 0 || m_relations_matrix.outerSize() != 0);
            }

        private:
            Mu_matrix        m_mu_matrix;
            Targets_matrix   m_targets_matrix;
            Relations_matrix m_relations_matrix;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_QP_SOLVER_DATA_H