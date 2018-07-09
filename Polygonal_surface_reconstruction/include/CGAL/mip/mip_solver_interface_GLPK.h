#include <CGAL/mip/mip_solver.h>
#include "3rd_glpk/glpk.h"


namespace CGAL {

	/// \cond SKIP_IN_MANUAL

	namespace internal {

		// infer "bound type" (required by GLPK) from the bounds values
		int bound_type(double lb, double ub) {
			if (lb <= -Variable::infinity() && ub >= Variable::infinity())
				return GLP_FR;		// free (unbounded) variable

			else if (lb > -Variable::infinity() && ub >= Variable::infinity())
				return GLP_LO;		// variable with lower bound

			else if (lb <= -Variable::infinity() && ub < Variable::infinity())
				return GLP_UP;		// variable with upper bound

			else {// lb > -Variable::infinity() && ub < Variable::infinity()
				if (lb == ub)
					return GLP_FX;	// fixed variable
				else
					return GLP_DB;  // double-bounded variable
			}
		}
	}
	/// \endcond

	bool MIP_solver::solve(const MIP_model* model) {
		try {
			error_message_.clear();

			if (!model->is_valid(true)) {
				error_message_ = "model is not valid.";
				return false;
			}

			glp_prob* lp = glp_create_prob();
			if (!lp) {
				error_message_ = "failed creating a LP model";
				return false;
			}

			std::size_t num_integer_variables = 0;

			// create variables
			const std::vector<Variable*>& variables = model->variables();

			// this suppresses many annoying warnings: "conversion from 'size_t' to 'int', possible loss of data"
			int num_variables = static_cast<int>(variables.size());

			glp_add_cols(lp, num_variables);
			for (int i = 0; i < num_variables; ++i) {
				const Variable* var = variables[i];
				glp_set_col_name(lp, i + 1, var->name().c_str());

				if (var->variable_type() == Variable::INTEGER) {
					glp_set_col_kind(lp, i + 1, GLP_IV);	// glpk uses 1-based arrays
					++num_integer_variables;
				}
				else if (var->variable_type() == Variable::BINARY) {
					glp_set_col_kind(lp, i + 1, GLP_BV);	// glpk uses 1-based arrays
					++num_integer_variables;
				}
				else 
					glp_set_col_kind(lp, i + 1, GLP_CV);	// continuous variable

				double lb, ub;
				var->get_bounds(lb, ub);

				int type = internal::bound_type(lb, ub);
				glp_set_col_bnds(lp, i + 1, type, lb, ub);
			}

			// Add constraints

			const std::vector<Linear_constraint*>& constraints = model->constraints();

			// this suppresses many annoying warnings: "conversion from 'size_t' to 'int', possible loss of data"
			int num_constraints = static_cast<int>(constraints.size());
			glp_add_rows(lp, num_constraints);

			for (int i = 0; i < num_constraints; ++i) {
				const Linear_constraint* c = constraints[i];
				const std::unordered_map<const Variable*, double>& coeffs = c->coefficients();
				std::unordered_map<const Variable*, double>::const_iterator cur = coeffs.begin();

				std::vector<int>	indices(coeffs.size() + 1, 0);		// glpk uses 1-based arrays
				std::vector<double> coefficients(coeffs.size() + 1, 0.0);  // glpk uses 1-based arrays
				std::size_t idx = 1; // glpk uses 1-based arrays
				for (; cur != coeffs.end(); ++cur) {
					int var_idx = cur->first->index();
					double coeff = cur->second;

					indices[idx] = var_idx + 1;	 // glpk uses 1-based arrays
					coefficients[idx] = coeff;
					++idx;
				}

				glp_set_mat_row(lp, i + 1, static_cast<int>(coeffs.size()), indices.data(), coefficients.data());

				double lb, ub;
				c->get_bounds(lb, ub);

				int type = internal::bound_type(lb, ub);
				glp_set_row_bnds(lp, i + 1, type, lb, ub);

				glp_set_row_name(lp, i + 1, c->name().c_str());
			}

			// set objective 

			// determine the coefficient of each variable in the objective function
			const Linear_objective * objective = model->objective();
			const std::unordered_map<const Variable*, double>& obj_coeffs = objective->coefficients();
			std::unordered_map<const Variable*, double>::const_iterator cur = obj_coeffs.begin();
			for (; cur != obj_coeffs.end(); ++cur) {
				int var_idx = cur->first->index();
				double coeff = cur->second;
				glp_set_obj_coef(lp, var_idx + 1, coeff); // glpk uses 1-based arrays
			}

			// Set objective function sense
			bool minimize = (objective->sense() == Linear_objective::MINIMIZE);
			glp_set_obj_dir(lp, minimize ? GLP_MIN : GLP_MAX);
			int msg_level = GLP_MSG_ERR;
			int status = -1;
			if (num_integer_variables == 0) { // continuous problem
				glp_smcp parm;
				glp_init_smcp(&parm);
				parm.msg_lev = msg_level;
				status = glp_simplex(lp, &parm);
			}
			else { // solve as MIP problem
				glp_iocp parm;
				glp_init_iocp(&parm);
				parm.msg_lev = msg_level;
				parm.presolve = GLP_ON;
				// The routine glp_intopt is a driver to the MIP solver based on the branch-and-cut method,
				// which is a hybrid of branch-and-bound and cutting plane methods.
				status = glp_intopt(lp, &parm);
			}

			switch (status) {
			case 0: {
				if (num_integer_variables == 0) { // continuous problem
					result_.resize(num_variables);
					for (int i = 0; i < num_variables; ++i) {
						result_[i] = glp_get_col_prim(lp, i + 1);	// glpk uses 1-based arrays
					}
				}
				else { // MIP problem
					result_.resize(num_variables);
					for (int i = 0; i < num_variables; ++i) {
						double x = glp_mip_col_val(lp, i + 1);		// glpk uses 1-based arrays
						Variable* v = variables[i];
						v->set_solution_value(x);
						if (v->variable_type() != Variable::CONTINUOUS)
							result_[i] = static_cast<int>(std::round(x));
					}
				}
				break;
			}

			case GLP_EBOUND:
				error_message_ =
					"Unable to start the search, because some double-bounded variables have incorrect"
					"bounds or some integer variables have non - integer(fractional) bounds.";
				break;

			case GLP_EROOT:
				error_message_ =
					"Unable to start the search, because optimal basis for initial LP relaxation is not"
					"provided. (This code may appear only if the presolver is disabled.)";
				break;

			case GLP_ENOPFS:
				error_message_ =
					"Unable to start the search, because LP relaxation of the MIP problem instance has"
					"no primal feasible solution. (This code may appear only if the presolver is enabled.)";
				break;

			case GLP_ENODFS:
				error_message_ =
					"Unable to start the search, because LP relaxation of the MIP problem instance has"
					"no dual feasible solution.In other word, this code means that if the LP relaxation"
					"has at least one primal feasible solution, its optimal solution is unbounded, so if the"
					"MIP problem has at least one integer feasible solution, its(integer) optimal solution"
					"is also unbounded. (This code may appear only if the presolver is enabled.)";
				break;

			case GLP_EFAIL:
				error_message_ = 
					"The search was prematurely terminated due to the solver failure.";
				break;

			case GLP_EMIPGAP:
				error_message_ =
					"The search was prematurely terminated, because the relative mip gap tolerance has been reached.";
				break;

			case GLP_ETMLIM:
				error_message_ =
					"The search was prematurely terminated, because the time limit has been exceeded.";
				break;

			case GLP_ESTOP:
				error_message_ =
					"The search was prematurely terminated by application. (This code may appear only"
					"if the advanced solver interface is used.)";
				break;

			default:
				error_message_ =
					"optimization was stopped with status code " + std::to_string(status);
				break;
			}

			glp_delete_prob(lp);

			return (status == 0);
		}
		catch (std::exception e) {
			error_message_ = e.what();
		}
		catch (...) {
			error_message_ = "exception during optimization.";
		}

		return false;
	}

}