#ifndef CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_OOQP_PROBLEM_H
#define CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_OOQP_PROBLEM_H

// STL includes.
#include <map>
#include <list>
#include <cassert>

// CGAL includes.
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_parameters.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_regular_segment.h>

// OOQP solver.
#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "GondzioSolver.h"
#include "QpGenSparseMa27.h"

namespace CGAL {

	namespace LOD {

        template<class KernelTraits, class QPProblemData>
		class Level_of_detail_segment_regularizer_ooqp_problem {

        public:
            typedef KernelTraits  Kernel;
            typedef QPProblemData QP_problem_data;

            using FT     = typename Kernel::FT;
            using Bounds = std::vector<FT>;
            
            using Solution_ft_type     = std::vector<FT>;
            using Solution_double_type = std::vector<double>;

            using Mus_matrix       = typename QP_problem_data::Mus_matrix;
            using Targets_matrix   = typename QP_problem_data::Targets_matrix;
            using Relations_matrix = typename QP_problem_data::Relations_matrix;

            using Mus_iterator     = typename Mus_matrix::InnerIterator;
            using Targets_iterator = typename Targets_matrix::InnerIterator;

            using Parameters = CGAL::LOD::Level_of_detail_segment_regularizer_parameters<Kernel>;

            using Regular_segment  = CGAL::LOD::Level_of_detail_segment_regularizer_regular_segment<Kernel>;
            using Regular_segments = std::vector<Regular_segment>;

            Level_of_detail_segment_regularizer_ooqp_problem(const Bounds &bounds, const QP_problem_data &qp_data, const Parameters &parameters, const Regular_segments &segments) 
            : m_bounds(bounds), m_qp_data(qp_data), m_parameters(parameters), m_segments(segments) { }

            void solve(Solution_ft_type &solution_ft_type) {

                // Allocate objective function.
                int *krowQ, *jcolQ;
                double *c, *dQ;
                allocate_objective_function(&krowQ, &jcolQ, &dQ, &c);


                // Allocate bounds.
                double *xlow, *xupp;
                char *ixlow, *ixupp;
                allocate_bounds(&xlow, &ixlow, &xupp, &ixupp);


                // Allocate inequality constraints.
                int *krowC, *jcolC;
                double *dC, *clow, *cupp;
                char *iclow, *icupp;
                allocate_inequality_constraints(&krowC, &jcolC, &dC, &clow, &iclow, &cupp, &icupp);


                // Some extra input data.
                int krowA = 0;

                const int num_individuals = static_cast<int>(m_qp_data.get_number_of_individuals());
                const int num_variables   = static_cast<int>(m_qp_data.get_number_of_variables());

                const Mus_matrix &mus_matrix = m_qp_data.get_mus_matrix();


                // Set up the OOQP solver.
                QpGenSparseMa27   *qp = new QpGenSparseMa27(num_variables, 0, 2 * mus_matrix.nonZeros(), num_individuals, 0, 6 * mus_matrix.nonZeros());
                QpGenData       *prob = (QpGenData *)qp->makeData(c, krowQ, jcolQ, dQ, xlow, ixlow, xupp, ixupp, &krowA, NULL, NULL, NULL, krowC, jcolC, dC, clow, iclow, cupp, icupp);
                QpGenVars       *vars = (QpGenVars *)qp->makeVariables(prob);
                QpGenResiduals *resid = (QpGenResiduals *)qp->makeResiduals(prob);
                GondzioSolver *solver = new GondzioSolver(qp, prob);


                // Solve.
                const int status = solver->solve(prob, vars, resid);

                if (status == 0) {
                    
                    Solution_double_type solution_double_type = Solution_double_type(num_variables, 0.0);
                    vars->x->copyIntoArray(solution_double_type.data());
                    
                    cast_to_ft(solution_double_type, solution_ft_type);

                } else {

                    std::string m_status;
                    switch (status) {

                        case 1:  m_status = "NOT_FINISHED"; break;
                        case 2:  m_status = "MAX_ITS_EXCEEDED"; break;
                        case 3:  m_status = "INFEASIBLE"; break;
                        case 4:  m_status = "UNKNOWN"; break;
                        default: m_status = "OTHER_STATUS"; break;
                    }
                    std::cerr << std::endl << "WARNING : solver status is " << m_status << " (code " << status << ")" << std::endl << std::endl;
                }


                // Deallocate memory.
                delete solver;
                delete resid;
                delete vars;
                delete prob;
                delete qp;

                deallocate_inequality_constraints(krowC, jcolC, dC, clow, iclow, cupp, icupp);
                deallocate_bounds(xlow, ixlow, xupp, ixupp);
                deallocate_objective_function(krowQ, jcolQ, dQ, c);
            }

        private:
            const Bounds           &m_bounds;
            const QP_problem_data  &m_qp_data;
            const Parameters       &m_parameters;
            const Regular_segments &m_segments;

            void allocate_objective_function(int **p_rowQ, int **p_colQ, double **p_dQ, double **p_c) {

                // Some dimensions and parameters.
                const int num_individuals = static_cast<int>(m_qp_data.get_number_of_individuals());
                const int num_variables   = static_cast<int>(m_qp_data.get_number_of_variables());
                
                const double lambda = CGAL::to_double(m_parameters.get_lambda());

                // Allocate quadratic term.
                *p_rowQ =    new int[num_variables + 1];
                *p_colQ =    new int[num_individuals];
                *p_dQ   = new double[num_individuals];

                int *rowQ = *p_rowQ;
                int *colQ = *p_colQ;

                double *dQ = *p_dQ;
                const double weight = 100000.0;

                for (int i = 0; i <= num_variables; ++i) {
                    rowQ[i] = (i < num_individuals ? i : num_individuals);

                    if (i < num_individuals) {
                        colQ[i] = i;

                        const double M = get_bound(i);
                        const double quadratic_term = weight * 2.0 * (1.0 - lambda) / (M * M * num_individuals);
                        dQ[i] = quadratic_term;
                    }
                }

                // Allocate linear term.
                *p_c      = new double[num_variables];
                double* c = *p_c;

                const double M = get_maximum_bound();
                for (int i = 0; i < num_variables; ++i) {

                    const double linear_term = weight * lambda / (4.0 * M * (num_variables - num_individuals));
                    const double c_i = (i < num_individuals ? 0.0 : linear_term);
                    c[i] = c_i;
                }
            }

            double get_maximum_bound() const {
                assert(m_bounds.size() > 0);

                FT max_bound = -FT(1000000000000);
                for (size_t i = 0; i < m_bounds.size(); ++i)
                    max_bound = CGAL::max(max_bound, m_bounds[i]);

                return CGAL::to_double(max_bound);
            }

            double get_bound(const size_t segment_index) const {
                return CGAL::to_double(m_bounds[segment_index]);
            }

            void allocate_bounds(double **p_xlow, char **p_ixlow, double **p_xupp, char **p_ixupp) {

                // Some dimensions.
                const int num_individuals = static_cast<int>(m_qp_data.get_number_of_individuals());
                const int num_variables   = static_cast<int>(m_qp_data.get_number_of_variables());

                // Define bounds for every variable.
                *p_xlow = new double[num_variables];
                *p_xupp = new double[num_variables];

                *p_ixlow = new char[num_variables];
                *p_ixupp = new char[num_variables];

                double* xlow = *p_xlow;
                double* xupp = *p_xupp;

                char* ixlow = *p_ixlow;
                char* ixupp = *p_ixupp;

                for (int i = 0; i < num_variables; ++i) {
                    if (i < num_individuals) {

                        xlow[i] = -get_bound(i);
                        xupp[i] =  get_bound(i);

                        ixlow[i] = 1;
                        ixupp[i] = 1;

                    } else {

                        xlow[i] = 0.0;
                        xupp[i] = 0.0;

                        ixlow[i] = 0;
                        ixupp[i] = 0;
                    }
                }
            }

            void allocate_inequality_constraints(int **p_rowC, int **p_colC, double **p_dC, double **p_clow, char **p_iclow, double **p_cupp, char **p_icupp) {
                
                // Some required data and dimensions.
                const int num_individuals = static_cast<int>(m_qp_data.get_number_of_individuals());
                
                const Mus_matrix         &mus_matrix = m_qp_data.get_mus_matrix();
                const Targets_matrix &targets_matrix = m_qp_data.get_targets_matrix();

                // Allocate the sparse matrix that corresponds to the constraints.
                *p_rowC =    new int[2 * mus_matrix.nonZeros() + 1];
                *p_colC =    new int[6 * mus_matrix.nonZeros()];
                *p_dC   = new double[6 * mus_matrix.nonZeros()];

                int  *rowC = *p_rowC;
                int  *colC = *p_colC;
                double *dC = *p_dC;

                for (int i = 0; i <= 2 * mus_matrix.nonZeros(); ++i) rowC[i] = 3 * i;

                int p = 0;
                for (int k = 0; k < mus_matrix.outerSize(); ++k) {
                    for (Mus_iterator it(mus_matrix, k); it; ++it) {

                        const int i = it.row();
                        const int j = it.col();

                        const double mu_ij = CGAL::to_double(it.value());

                        colC[6 * p + 0] = i;
                        colC[6 * p + 1] = j;
                        colC[6 * p + 2] = num_individuals + p;
                        colC[6 * p + 3] = i;
                        colC[6 * p + 4] = j;
                        colC[6 * p + 5] = num_individuals + p;

                        dC[6 * p + 0] = -2.0 * mu_ij;
                        dC[6 * p + 1] =  2.0 * mu_ij;
                        dC[6 * p + 2] = -1.0;
                        dC[6 * p + 3] =  2.0 * mu_ij;
                        dC[6 * p + 4] = -2.0 * mu_ij;
                        dC[6 * p + 5] = -1.0;

                        ++p;
                    }
                }
                
                // Define clow <= Cx <= cupp.
                *p_clow = new double[2 * mus_matrix.nonZeros()];
                *p_cupp = new double[2 * mus_matrix.nonZeros()];

                *p_iclow = new char[2 * mus_matrix.nonZeros()];
                *p_icupp = new char[2 * mus_matrix.nonZeros()];

                double *clow = *p_clow;
                double *cupp = *p_cupp;

                char *iclow = *p_iclow;
                char *icupp = *p_icupp;

                p = 0;
                for (int k = 0; k < mus_matrix.outerSize(); ++k) {

                    Mus_iterator     it_mus(mus_matrix, k);
                    Targets_iterator it_targets(targets_matrix, k);

                    while (it_mus && it_targets) {

                        iclow[2 * p + 0] = 0;
                        iclow[2 * p + 1] = 0;

                        clow[2 * p + 0] = 0.0;
                        clow[2 * p + 1] = 0.0;

                        icupp[2 * p + 0] = 1;
                        icupp[2 * p + 1] = 1;

                        cupp[2 * p + 0] = -2.0 * CGAL::to_double(it_mus.value()) * CGAL::to_double(it_targets.value());
                        cupp[2 * p + 1] =  2.0 * CGAL::to_double(it_mus.value()) * CGAL::to_double(it_targets.value());

                        ++it_mus;
                        ++it_targets;
                        ++p;
                    }
                }
            }

            void cast_to_ft(const Solution_double_type &solution_double_type, Solution_ft_type &solution_ft_type) {
                
                solution_ft_type.clear();
                solution_ft_type.resize(solution_double_type.size(), FT(0));

                for (size_t i = 0; i < solution_double_type.size(); ++i)
                    solution_ft_type[i] = static_cast<FT>(solution_double_type[i]);
            }

            void deallocate_objective_function(int *rowQ, int *colQ, double *dQ, double *c) {
                delete[] rowQ;
                delete[] colQ;
                delete[] dQ;
                delete[] c;
            }

            void deallocate_bounds(double *xlow, char *ixlow, double *xupp, char *ixupp) {
                delete[] xlow;
                delete[] ixlow;
                delete[] xupp;
                delete[] ixupp;
            }

            void deallocate_inequality_constraints(int *rowC, int *colC, double *dC, double *clow, char *iclow, double *cupp, char *icupp) {
                delete[] rowC;
                delete[] colC;
                delete[] dC;
                delete[] clow;
                delete[] iclow;
                delete[] cupp;
                delete[] icupp;
            }
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_SEGMENT_REGULARIZER_OOQP_PROBLEM_H