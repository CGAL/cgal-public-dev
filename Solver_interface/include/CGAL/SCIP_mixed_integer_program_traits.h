// Copyright (c) 2018  Liangliang Nan
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Liangliang Nan

#ifndef CGAL_GLPK_MIXED_INTEGER_PROGRAM_TRAITS_H
#define CGAL_GLPK_MIXED_INTEGER_PROGRAM_TRAITS_H

#include <CGAL/Mixed_integer_program_traits.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

namespace CGAL {

	/// \ingroup PkgSolver
	///
	/// The class `SCIP_mixed_integer_program_traits` provides an interface for 
	/// formulating and solving (constrained) mixed integer programs (It can 
	/// also be used for general linear programs) using \ref thirdpartySCIP.
	///
	/// \tparam FT Number type
	///
	/// \cgalModels `MixedIntegerProgramTraits`

	template <typename FT>
	class SCIP_mixed_integer_program_traits : public Mixed_integer_program_traits<FT>
	{
	public:
		typedef CGAL::Variable<FT>						Variable;
		typedef CGAL::Linear_constraint<FT>				Linear_constraint;
		typedef CGAL::Linear_objective<FT>				Linear_objective;
		typedef typename Linear_objective::Sense		Sense;
		typedef typename Variable::Variable_type		Variable_type;

	public:

		/// Solves the program. Returns false if fails.
		virtual bool solve();
	};


	//////////////////////////////////////////////////////////////////////////

	// implementation

	template<typename FT>
	bool SCIP_mixed_integer_program_traits<FT>::solve() {

		error_message_.clear();

		Scip* scip = 0;
		SCIP_CALL(SCIPcreate(&scip));
		SCIP_CALL(SCIPincludeDefaultPlugins(scip));

		// disable scip output to stdout
		SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scip), TRUE);

		// use wall clock time because getting CPU user seconds
		// involves calling times() which is very expensive
		SCIP_CALL(SCIPsetIntParam(scip, "timing/clocktype", SCIP_CLOCKTYPE_WALL));

		// create empty problem 
		SCIP_CALL(SCIPcreateProbBasic(scip, "Polygonal_surface_reconstruction"));

		// create variables
		std::vector<SCIP_VAR*> scip_variables;
		for (std::size_t i = 0; i < variables_.size(); ++i) {
			const Variable* var = variables_[i];
			SCIP_VAR* v = 0;

			double lb, ub;
			var->get_bounds(lb, ub);

			switch (var->variable_type())
			{
			case Variable::CONTINUOUS:
				SCIP_CALL(SCIPcreateVar(scip, &v, var->name().c_str(), lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, 0, 0, 0, 0, 0));
				break;
			case Variable::INTEGER:
				SCIP_CALL(SCIPcreateVar(scip, &v, var->name().c_str(), lb, ub, 0.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE, 0, 0, 0, 0, 0));
				break;
			case Variable::BINARY:
				SCIP_CALL(SCIPcreateVar(scip, &v, var->name().c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, 0, 0, 0, 0, 0));
				break;
			}
			// add the SCIP_VAR object to the scip problem
			SCIP_CALL(SCIPaddVar(scip, v));

			// storing the SCIP_VAR pointer for later access
			scip_variables.push_back(v);
		}

		// Add constraints

		std::vector<SCIP_CONS*> scip_constraints;
		for (std::size_t i = 0; i < constraints_.size(); ++i) {
			const Linear_constraint* c = constraints_[i];
			const std::unordered_map<const Variable*, double>& coeffs = c->coefficients();
			std::unordered_map<const Variable*, double>::const_iterator cur = coeffs.begin();

			std::vector<SCIP_VAR*>	cstr_variables(coeffs.size());
			std::vector<double>		cstr_values(coeffs.size());
			std::size_t idx = 0;
			for (; cur != coeffs.end(); ++cur) {
				std::size_t var_idx = cur->first->index();
				double coeff = cur->second;
				cstr_variables[idx] = scip_variables[var_idx];
				cstr_values[idx] = coeff;
				++idx;
			}

			// create SCIP_CONS object
			SCIP_CONS* cons = 0;
			const std::string& name = c->name();

			double lb, ub;
			c->get_bounds(lb, ub);

			SCIP_CALL(SCIPcreateConsLinear(scip, &cons, name.c_str(), int(coeffs.size()), cstr_variables.data(), cstr_values.data(), lb, ub, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
			SCIP_CALL(SCIPaddCons(scip, cons));			// add the constraint to scip

			// store the constraint for later on
			scip_constraints.push_back(cons);
		}

		// set objective

		// determine the coefficient of each variable in the objective function
		const std::unordered_map<const Variable*, double>& obj_coeffs = objective_->coefficients();
		std::unordered_map<const Variable*, double>::const_iterator cur = obj_coeffs.begin();
		for (; cur != obj_coeffs.end(); ++cur) {
			const Variable* var = cur->first;
			double coeff = cur->second;
			SCIP_CALL(SCIPchgVarObj(scip, scip_variables[var->index()], coeff));
		}

		// set the objective sense
		bool minimize = (objective_->sense() == Linear_objective::MINIMIZE);
		SCIP_CALL(SCIPsetObjsense(scip, minimize ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE));

		// Always turn presolve on (it's the SCIP default).
		bool presolve = true;
		if (presolve)
			SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrounds", -1)); // maximal number of presolving rounds (-1: unlimited, 0: off)
		else
			SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrounds", 0));  // disable presolve

		bool status = false;
		// this tells scip to start the solution process
		if (SCIPsolve(scip) == SCIP_OKAY) {
			// get the best found solution from scip
			SCIP_SOL* sol = SCIPgetBestSol(scip);
			if (sol) {
				// If optimal or feasible solution is found.
				result_.resize(variables_.size());
				for (std::size_t i = 0; i < variables_.size(); ++i) {
					double x = SCIPgetSolVal(scip, sol, scip_variables[i]);
					Variable* v = variables_[i];
					v->set_solution_value(x);
					if (v->variable_type() != Variable::CONTINUOUS)
						result_[i] = static_cast<int>(std::round(x));
				}
				status = true;
			}
		}

		// report the status: optimal, infeasible, etc.
		SCIP_STATUS scip_status = SCIPgetStatus(scip);
		switch (scip_status) {
		case SCIP_STATUS_OPTIMAL:
			// provides info only if fails.
			break;
		case SCIP_STATUS_GAPLIMIT:
			// To be consistent with the other solvers.
			// provides info only if fails.
			break;
		case SCIP_STATUS_INFEASIBLE:
			error_message_ = "model was infeasible";
			break;
		case SCIP_STATUS_UNBOUNDED:
			error_message_ = "model was unbounded";
			break;
		case SCIP_STATUS_INFORUNBD:
			error_message_ = "model was either infeasible or unbounded";
			break;
		case SCIP_STATUS_TIMELIMIT:
			error_message_ = "aborted due to time limit";
			break;
		default:
			error_message_ = "aborted with status: " + std::to_string(scip_status);
			break;
		}

		SCIP_CALL(SCIPresetParams(scip));

		// since the SCIPcreateVar captures all variables, we have to release them now
		for (std::size_t i = 0; i < scip_variables.size(); ++i)
			SCIP_CALL(SCIPreleaseVar(scip, &scip_variables[i]));
		scip_variables.clear();

		// the same for the constraints
		for (std::size_t i = 0; i < scip_constraints.size(); ++i)
			SCIP_CALL(SCIPreleaseCons(scip, &scip_constraints[i]));
		scip_constraints.clear();

		// after releasing all vars and cons we can free the scip problem
		// remember this has always to be the last call to scip
		SCIP_CALL(SCIPfree(&scip));

		return status;
	}

} // namespace CGAL

#endif // CGAL_GLPK_MIXED_INTEGER_PROGRAM_TRAITS_H
