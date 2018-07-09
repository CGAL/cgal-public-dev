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

#ifndef CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_MIP_MODEL_H
#define CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_MIP_MODEL_H

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

/*!
\file mip_model.h
*/

namespace CGAL {

	class MIP_model;

	/// The base class of linear model elements, e.g., Variable, Linear_constraint, and Linear_objective 
	class Model_entry
	{
	private:
		/// A model entry (e.g., variable, constraint, objective) cannot belong to multiple models.
		/// "model" owns this entry.
		Model_entry(MIP_model* model, const std::string& name = "", int idx = 0) : model_(model), name_(name), index_(idx) {}

	public:
		const std::string& name() const { return name_; }
		void set_name(const std::string& n) { name_ = n; }

		int  index() const { return index_; }
		void set_index(int idx) { index_ = idx; }

		/// the model that owns this entry
		const MIP_model* model() const { return model_; }
		MIP_model* model() { return model_; }

	private:
		MIP_model * model_; // the model that owns this entry
		std::string	name_;
		int			index_;

		friend class Variable;
		friend class Linear_expression;
		friend class MIP_model;
		friend class MIP_solver;
	};

	
	/// The base class of linear model elements that might have bound constraints.
	class Bound
	{
	private:
		Bound(double lb = -infinity(), double ub = +infinity());

	public:
		void   set_bounds(double lb, double ub);
		void   set_lower_bound(double lb) { lower_bound_ = lb; }
		void   set_upper_bound(double ub) { upper_bound_ = ub; }

		void   get_bounds(double& lb, double& ub) const;
		double lower_bound() const { return lower_bound_; }
		double upper_bound() const { return upper_bound_; }

		static double infinity();

	private:
		double		lower_bound_;
		double		upper_bound_;

		static double infinity_;

		friend class Variable;
		friend class Linear_constraint;
		friend class MIP_model;
		friend class MIP_solver;
	};


	/// The variables of linear models.
	class Variable : public Model_entry, public Bound
	{
	public:
		enum VariableType { CONTINUOUS, INTEGER, BINARY };

	private:
		/// A variable cannot belong to several models.
		/// "model" owns this variable.
		Variable(MIP_model* model, VariableType type = CONTINUOUS, double lb = -infinity(), double ub = +infinity(), const std::string& name = "", int idx = 0);

	public:
		VariableType variable_type() const { return variable_type_; }
		void set_variable_type(VariableType t);

		/// Returns the value of the variable in the current solution.
		/// Note: (1) valid only if the problem was successfully solved.
		///       (2) if the variable is integer and rounded == true, then the 
		///           value will be rounded to the nearest integer.
		double solution_value(bool rounded = false) const;

	private:
		void set_solution_value(double value) { solution_value_ = value; }

	private:
		VariableType variable_type_;
		double		 solution_value_;

		friend class MIP_model;
		friend class MIP_solver;
	};


	/// The base class of Linear_constraint and Linear_objective.
	class Linear_expression : public Model_entry
	{
	private:
		/// An expression cannot belong to several models.
		/// "model" owns this expression.
		Linear_expression(MIP_model* model, const std::string& name = "", int idx = 0);

	public:
		/// Add a coefficient to a variable. 
		void  add_coefficient(const Variable* var, double coeff);

		const std::unordered_map<const Variable*, double>& coefficients() const { return coefficients_; }
		void  set_coefficients(const std::unordered_map<const Variable*, double>& coeffs) { coefficients_ = coeffs; }

		double get_coefficient(const Variable* var) const;

		// the constant term
		void set_offset(double value) { offset_ = value; }
		double offset() const { return offset_; }

		/// evaluates the value of this expression at the solution found.
		/// Note: (1) valid only if the problem was successfully solved.
		///       (2) if a variable is integer and rounded == true, then the 
		///           variable value will be rounded to the nearest integer.
		double solution_value(bool rounded = false) const;

		virtual void clear() { coefficients_.clear(); offset_ = 0.0; }

	private:
		std::unordered_map<const Variable*, double>	coefficients_;
		double								offset_;

		friend class Linear_constraint;
		friend class Linear_objective;
		friend class MIP_model;
		friend class MIP_solver;
	};

	
	/// The linear constraint of linear models.
	class Linear_constraint : public Linear_expression, public Bound
	{
	private:
		/// A constraint cannot belong to several models.
		/// "model" owns this constraint.
		Linear_constraint(MIP_model* model, double lb = -infinity(), double ub = +infinity(), const std::string& name = "", int idx = 0);
		virtual ~Linear_constraint() {}

		friend class MIP_model;
		friend class MIP_solver;
	};


	/// The linear objective of linear models.
	class Linear_objective : public Linear_expression
	{
	public:
		enum Sense { MINIMIZE, MAXIMIZE, UNDEFINED };

	private:
		/// An objective cannot belong to several models.
		/// "model" owns this objective.
		Linear_objective(MIP_model* model, Sense sense);
		virtual ~Linear_objective() {}

	public:
		void  set_sense(Sense sense) { sense_ = sense; }
		Sense sense() const { return sense_; }

		void clear();

	private:
		Sense sense_;

		friend class MIP_model;
		friend class MIP_solver;
	};


	/*!
	\ingroup PkgPolygonalSurfaceReconstruction

	The interface of modeling Mixed Integer Programs and general Linear Programs.

	*/

	class MIP_model
	{
	public:
		MIP_model();
		~MIP_model();

		/// create a single variable, it to the model, and returns the pointer.
		/// Note: if name is empty or not provided, a default name (e.g., x0, x1...) will be given.
		Variable* create_variable(
			Variable::VariableType type = Variable::CONTINUOUS,
			double lb = -Variable::infinity(),
			double ub = +Variable::infinity(),
			const std::string& name = ""
		);

		/// create a set of variables and add them to the model.
		/// Note: variables will be given default names, e.g., x0, x1...
		std::vector<Variable*> create_n_variables(std::size_t n);

		/// create a single linear constraint, add it to the model, and returns the pointer.
		/// Note: if name is empty or not provided, a default name (e.g., c0, c1...) will be given.
		Linear_constraint* create_constraint(
			double lb = -Variable::infinity(),
			double ub = +Variable::infinity(),
			const std::string& name = ""
		);

		/// create a set of linear constraints and add them to the model.	
		/// Note: constraints with be given default names, e.g., c0, c1...
		std::vector<Linear_constraint*> create_n_constraints(std::size_t n);

		/// create the objective function and returns the pointer.
		Linear_objective * create_objective(Linear_objective::Sense sense = Linear_objective::MINIMIZE);


		bool has_variable(const Variable* var) const;
		bool has_constraint(const Linear_constraint* cons) const;


		std::size_t num_variables() const { return variables_.size(); }
		const std::vector<Variable*>& variables() const { return variables_; }
		std::vector<Variable*>& variables() { return variables_; }

		std::size_t num_constraints() const { return constraints_.size(); }
		const std::vector<Linear_constraint*>& constraints() const { return constraints_; }
		std::vector<Linear_constraint*>& constraints() { return constraints_; }

		const Linear_objective * objective() const;
		Linear_objective * objective();

		std::size_t num_continuous_variables() const;
		std::size_t num_integer_variables() const;
		std::size_t num_binary_variables() const;

		bool is_continuous() const;			// returns true if all variables are continuous
		bool is_mix_integer_model() const;	// returns true if mixed inter model
		bool is_integer_model() const;		// returns true if inter model
		bool is_binary_proram() const;		// returns true if binary model

											/// print statistics of the model to the stream
		void print_statistics(std::ostream& output) const;

		/// check the model if there are issues like:
		///  -) variables, constraints, and/or objective not owned by the program
		///  -) duplicated variables or constraints
		///  -) variables have the same name or index
		///  -) constraints have the same name or index
		///  -) variables with infeasible bounds (e.g., lb > ub)
		///  -) constraints with infeasible bounds (e.g., lb > ub)
		bool is_valid(bool verbose = true) const;

		// clear all variables, constraints, and the objective.
		void clear();

	private:
		Linear_objective * objective_;
		std::vector<Variable*>			variables_;
		std::vector<Linear_constraint*>	constraints_;

		friend class MIP_solver;
	};


	//////////////////////////////////////////////////////////////////////////

	// implementation

	//double Bounded::infinity_ = std::numeric_limits<double>::max();
	double Bound::infinity_ = 1e20;		// values larger than this value are considered infinity


	Bound::Bound(double lb /* = -infinity() */, double ub /* = +infinity() */)
		: lower_bound_(lb)
		, upper_bound_(ub)
	{
	}

	double Bound::infinity() {
		return infinity_;
	}

	void Bound::set_bounds(double lb, double ub) {
		lower_bound_ = lb;
		upper_bound_ = ub;
	}

	void Bound::get_bounds(double& lb, double& ub) const {
		lb = lower_bound_;
		ub = upper_bound_;
	}


	Variable::Variable(
		MIP_model* model,
		VariableType type /* = CONTINUOUS */,
		double lb /* = -infinity() */,
		double ub /* = +infinity() */,
		const std::string& name /* = "" */,
		int idx /* = 0*/
	)
		: Model_entry(model, name, idx)
		, Bound(lb, ub)
		, variable_type_(type)
		, solution_value_(0.0)
	{
		if (type == BINARY)
			Bound::set_bounds(0.0, 1.0);
	}


	void Variable::set_variable_type(VariableType type) {
		variable_type_ = type;
		if (type == BINARY)
			Bound::set_bounds(0.0, 1.0);
	}


	double Variable::solution_value(bool rounded /* = false*/) const {
		if (rounded && variable_type_ != CONTINUOUS)
			return std::round(solution_value_);
		else
			return solution_value_;
	}


	//////////////////////////////////////////////////////////////////////////

	Linear_expression::Linear_expression(MIP_model* model, const std::string& name, int idx)
		: Model_entry(model, name, idx)
		, offset_(0.0)
	{
	}


	void Linear_expression::add_coefficient(const Variable* var, double coeff) {
		if (!model()->has_variable(var)) {
			std::cerr << "model does not own variable " << var->name() << " (" << var->index() << ")" << std::endl;
			return;
		}

		if (coefficients_.find(var) == coefficients_.end())
			coefficients_[var] = coeff;
		else
			coefficients_[var] += coeff;
	}


	double Linear_expression::get_coefficient(const Variable* var) const {
		if (!model()->has_variable(var)) {
			std::cerr << "model does not own variable " << var->name() << " (" << var->index() << ")" << std::endl;
			return 0.0;
		}

		std::unordered_map<const Variable*, double>::const_iterator pos = coefficients_.find(var);
		if (pos != coefficients_.end())
			return pos->second;
		else {
			std::cerr << "linear expression does not own variable " << var->name() << " (" << var->index() << ")" << std::endl;
			return 0.0;
		}
	}


	double Linear_expression::solution_value(bool rounded /* = false*/) const {
		double solution = offset_;

		std::unordered_map<const Variable*, double>::const_iterator it = coefficients_.begin();
		for (; it != coefficients_.end(); ++it) {
			const Variable* var = it->first;
			double coeff = it->second;
			solution += var->solution_value(rounded) * coeff;
		}
		return solution;
	}


	void Linear_objective::clear() {
		Linear_expression::clear();
		set_name("");
		set_index(0);
		set_offset(0.0);
	}


	Linear_constraint::Linear_constraint(
		MIP_model* model,
		double lb /* = -infinity() */,
		double ub /* = +infinity() */,
		const std::string& name/* = "" */,
		int idx /* = 0*/
	)
		: Linear_expression(model, name, idx)
		, Bound(lb, ub)
	{
	}


	Linear_objective::Linear_objective(MIP_model* model, Sense sense)
		: Linear_expression(model)
		, sense_(sense)
	{
	}


	MIP_model::MIP_model() {
		// intentionally set the objective to UNDEFINED, so it will allow me to warn
		// the user if he/she forgot to set the objective sense.
		objective_ = new Linear_objective(this, Linear_objective::UNDEFINED);
	}


	MIP_model::~MIP_model() {
		clear();
		delete objective_;
	}


	void MIP_model::clear() {
		for (std::size_t i = 0; i < variables_.size(); ++i)
			delete variables_[i];
		variables_.clear();

		for (std::size_t i = 0; i < constraints_.size(); ++i)
			delete constraints_[i];
		constraints_.clear();

		objective_->clear();
	}

	/// \cond SKIP_IN_MANUAL
	namespace internal {
		/**
		* Converts an integer v to a string of specified 'width' by
		* filling with character 'fill'
		*/
		template <class Int>
		inline std::string from_integer(Int v, int width, char fill) {
			std::ostringstream string_stream;
			string_stream << std::setfill(fill) << std::setw(width) << v;
			return string_stream.str();
		}
	}
	/// \endcond


	Variable* MIP_model::create_variable(
		Variable::VariableType type /* = Variable::CONTINUOUS */,
		double lb /* = -Variable::infinity() */,
		double ub /* = +Variable::infinity() */,
		const std::string& name /* = "" */)
	{
		Variable* v = new Variable(this, type, lb, ub);

		std::size_t idx = variables_.size();
		v->set_index(static_cast<int>(idx));

		const std::string& fixed_name = name.empty() ? "x" + internal::from_integer(idx, 9, '0') : name;
		v->set_name(fixed_name);

		variables_.push_back(v);
		return v;
	}


	std::vector<Variable*> MIP_model::create_n_variables(std::size_t n) {
		std::vector<Variable*> variables;
		for (std::size_t i = 0; i < n; ++i) {
			Variable* v = create_variable();
			variables.push_back(v);
		}
		return variables;
	}


	Linear_constraint* MIP_model::create_constraint(
		double lb /* = -Variable::infinity() */,
		double ub /* = +Variable::infinity() */,
		const std::string& name /* = "" */)
	{
		Linear_constraint* c = new Linear_constraint(this, lb, ub);

		std::size_t idx = constraints_.size();
		c->set_index(static_cast<int>(idx));

		const std::string& fixed_name = name.empty() ? "c" + internal::from_integer(idx, 9, '0') : name;
		c->set_name(fixed_name);

		constraints_.push_back(c);
		return c;
	}


	std::vector<Linear_constraint*> MIP_model::create_n_constraints(std::size_t n) {
		std::vector<Linear_constraint*> constraints;
		for (std::size_t i = 0; i < n; ++i) {
			Linear_constraint* v = create_constraint();
			constraints.push_back(v);
		}
		return constraints;
	}


	Linear_objective * MIP_model::create_objective(Linear_objective::Sense sense /* = Linear_objective ::MINIMIZE*/) {
		if (objective_)
			delete objective_;

		objective_ = new Linear_objective(this, sense);
		return objective_;
	}


	bool MIP_model::has_variable(const Variable* var) const {
		if (var == nullptr)
			return false;

		if (var->index() >= 0 && var->index() < variables_.size()) {
			// Then, verify that the variable with this index has the same address.
			return variables_[var->index()] == var;
		}
		return false;
	}


	bool MIP_model::has_constraint(const Linear_constraint* cons) const {
		if (cons == nullptr)
			return false;

		if (cons->index() >= 0 && cons->index() < constraints_.size()) {
			// Then, verify that the constraint with this index has the same address.
			return constraints_[cons->index()] == cons;
		}
		return false;
	}


	const Linear_objective * MIP_model::objective() const {
		if (!objective_)
			std::cerr << "please call \'create_objective()\' to create an objective first" << std::endl;

		return objective_;
	}


	Linear_objective * MIP_model::objective() {
		if (!objective_)
			std::cerr << "please call \'create_objective()\' to create an objective first" << std::endl;

		return objective_;
	}


	std::size_t MIP_model::num_continuous_variables() const {
		std::size_t num_continuous_var = 0;
		for (std::size_t i = 0; i < variables_.size(); ++i) {
			const Variable* v = variables_[i];
			if (v->variable_type() == Variable::CONTINUOUS)
				++num_continuous_var;
		}
		return num_continuous_var;
	}


	std::size_t MIP_model::num_integer_variables() const {
		std::size_t num_iteger_var = 0;
		for (std::size_t i = 0; i < variables_.size(); ++i) {
			const Variable* v = variables_[i];
			if (v->variable_type() == Variable::INTEGER/* || v->variable_type() == Variable::BINARY*/)
				++num_iteger_var;
		}
		return num_iteger_var;
	}


	std::size_t MIP_model::num_binary_variables() const {
		std::size_t num_binary_var = 0;
		for (std::size_t i = 0; i < variables_.size(); ++i) {
			const Variable* v = variables_[i];
			if (v->variable_type() == Variable::BINARY)
				++num_binary_var;
		}
		return num_binary_var;
	}


	// returns true if all variables are continuous
	bool MIP_model::is_continuous() const {
		std::size_t num = num_continuous_variables();
		return (num > 0) && (num == variables_.size());
	}


	// returns true if mixed inter model
	bool MIP_model::is_mix_integer_model() const {
		std::size_t num = num_continuous_variables();
		return (num > 0) && (num < variables_.size());
	}


	// returns true if inter model
	bool MIP_model::is_integer_model() const {
		std::size_t num = num_integer_variables();
		return (num > 0) && (num == variables_.size()); // or (num_continuous_variables() == 0)
	}


	// returns true if binary model
	bool MIP_model::is_binary_proram() const {
		std::size_t num = num_binary_variables();
		return (num > 0) && (num == variables_.size());
	}


	// print statistics of the model
	void MIP_model::print_statistics(std::ostream& output) const {
		std::cout << "\tProgram type: ";
		if (is_binary_proram())
			output << "Binary" << std::endl;
		else if (is_integer_model())
			output << "Integer" << std::endl;
		else if (is_mix_integer_model())
			output << "Mixed Integer" << std::endl;
		else
			output << "Continuous" << std::endl;
		output << "\t#Variables: " << num_variables() << std::endl;
		output << "\t\tBinary: " << num_binary_variables() << std::endl;
		output << "\t\tInteger: " << num_integer_variables() << std::endl;
		output << "\t\tContinuous: " << num_continuous_variables() << std::endl;
		output << "\t#Constraints: " << num_constraints() << std::endl;
		output << "\t#Objective sense: " << ((objective()->sense() == Linear_objective::MINIMIZE) ? "Minimize" : "Maximize") << std::endl;
		output << "\t#Objective offset: " << objective()->offset() << std::endl;
	}


	bool MIP_model::is_valid(bool verbose /* = true*/) const {
		bool valid = true;

		if (objective()->sense() == Linear_objective::UNDEFINED) {
			valid = false;
			if (verbose)
				std::cerr << "incomplete objective: undefined objective sense." << std::endl;
		}

		if (variables_.empty()) {
			valid = false;
			if (verbose)
				std::cerr << "variable set is empty" << std::endl;
		}

		// check if multiple variables have the same name or index
		std::unordered_map<int, const Variable*> vindices;
		std::unordered_map<std::string, const Variable*> vnames;
		for (std::size_t i = 0; i < variables_.size(); ++i) {
			const Variable* v = variables_[i];

			const std::string& name = v->name();
			int idx = v->index();

			if (v->model() != this) {
				valid = false;
				if (verbose)
					std::cerr << "variable " << v->name() << " (index " << idx << ") is not owned by this program" << std::endl;
			}

			if (idx != i) {
				valid = false;
				if (verbose)
					std::cerr << "variable " << v->name() << " (index " << idx << ") has an inconsistent index" << std::endl;
			}

			if (v->lower_bound() > v->upper_bound()) {
				valid = false;
				if (verbose)
					std::cerr << "variable " << v->name() << " (index " << idx << ") has contradictory bounds (i.e. lb > ub): "
					<< " lb = " << v->lower_bound() << ", ub = " << v->upper_bound() << std::endl;
			}

			std::unordered_map<int, const Variable*>::const_iterator vpos_id = vindices.find(idx);
			if (vpos_id == vindices.end())
				vindices[idx] = v;
			else {
				valid = false;
				if (verbose) {
					if (vpos_id->second == v)
						std::cerr << "duplicated variable " << name << " (index " << idx << ")" << std::endl;
					else
						std::cerr << "variables " << vpos_id->second->name() << " and " << v->name() << " have the same index " << idx << std::endl;
				}
			}

			std::unordered_map<std::string, const Variable*>::const_iterator vpos_name = vnames.find(name);
			if (vpos_name == vnames.end())
				vnames[name] = v;
			else {
				valid = false;
				if (verbose) {
					if (vpos_name->second == v)
						std::cerr << "duplicated variable " << name << " (index " << idx << ")" << std::endl;
					else
						std::cerr << "variables " << vpos_name->second->index() << " and " << v->index() << " have the same name: " << vpos_name->first << std::endl;
				}
			}
		}

		// check if multiple constraints have the same name or index
		std::unordered_map<int, const Linear_constraint*> cindices;
		std::unordered_map<std::string, const Linear_constraint*> cnames;
		for (std::size_t i = 0; i < constraints_.size(); ++i) {
			const Linear_constraint* c = constraints_[i];
			const std::string& name = c->name();
			int idx = c->index();

			if (c->model() != this) {
				valid = false;
				if (verbose)
					std::cerr << "constraint " << c->name() << " (index " << idx << ") is not owned by this program" << std::endl;
			}

			if (idx != i) {
				valid = false;
				if (verbose)
					std::cerr << "constraint indices are not consistent" << std::endl;
			}

			if (c->lower_bound() > c->upper_bound()) {
				valid = false;
				if (verbose)
					std::cerr << "constraint " << c->name() << " (index " << idx << ") has contradictory bounds (i.e. lb > ub): "
					<< " lb = " << c->lower_bound() << ", ub = " << c->upper_bound() << std::endl;
			}

			std::unordered_map<int, const Linear_constraint*>::const_iterator cpos_id = cindices.find(idx);
			if (cpos_id == cindices.end())
				cindices[idx] = c;
			else {
				valid = false;
				if (verbose) {
					if (cpos_id->second == c)
						std::cerr << "duplicated constraint " << name << " (index " << idx << ")" << std::endl;
					else
						std::cerr << "constraints " << cpos_id->second->name() << " and " << name << " have the same index " << idx << std::endl;
				}
			}

			std::unordered_map<std::string, const Linear_constraint*>::const_iterator cpos_name = cnames.find(name);
			if (cpos_name == cnames.end())
				cnames[name] = c;
			else {
				valid = false;
				if (verbose) {
					if (cpos_name->second == c)
						std::cerr << "duplicated constraint " << name << " (index " << idx << ")" << std::endl;
					else
						if (verbose)
							std::cerr << "constraints " << cpos_name->second->index() << " and " << c->index() << " have the same name: " << cpos_name->first << std::endl;
				}
			}
		}

		return valid;
	}


} //namespace CGAL

#endif	// CGAL_POLYGONAL_SURFACE_RECONSTRUCTION_MIP_MODEL_H
