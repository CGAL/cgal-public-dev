/*!
\ingroup PkgSolverConcepts
\cgalConcept

@brief Concept describing the set of requirements for (constrained or unconstrained) 
Mixed Integer Programming (MIP) problems. A model of this concept stores the integer
variables, linear objective, and linear constraints (if any) and provides a method 
to solve the problem.

\cgalHasModel `CGAL::Mixed_integer_program_traits<T>`
\cgalHasModel `CGAL::GLPK_mixed_integer_program_traits<T>`
\cgalHasModel `CGAL::SCIP_mixed_integer_program_traits<T>`
*/

class MixedIntegerProgramTraits
{
public:
	/// \name Types
	/// @{

	/*!

	*/
	typedef unspecified_type Variable;

	/*!

	*/
	typedef unspecified_type Linear_constraint;

	/*!

	*/
	typedef unspecified_type Linear_objective;

	/// @}

	/// \name Creation
	/// @{

	/*!
	Default constructor.
	*/
	MixedIntegerProgramTraits();

	/// @}

	/// \name Operations
	/// @{

	/// Creates a single variable, adds it to the solver, and returns its pointer.
	Variable* create_variable(Variable_type type, FT lb, FT ub, const std::string& name);

	/// Creates a set of variables, adds them to the solver, and returns their pointers.
	/// \note Variables will be given default names, e.g., x0, x1...
	std::vector<Variable*> create_n_variables(std::size_t n);

	/// Creates a single linear constraint, adds it to the solver, and returns the pointer.
	Linear_constraint* create_constraint(FT lb, FT ub, const std::string& name);

	/// Creates a set of linear constraints, adds them to the solver, and returns their pointers.	
	/// \note Constraints will be given default names, e.g., c0, c1...
	std::vector<Linear_constraint*> create_n_constraints(std::size_t n);

	/// Creates the objective function and returns the pointer.
	Linear_objective * create_objective(Sense sense);

	/// Returns the number of variables
	std::size_t num_variables() const;

	/// Returns the variables
	const std::vector<Variable*>& variables() const;
	std::vector<Variable*>& variables();

	/// Returns the number of constraints
	std::size_t num_constraints() const;

	/// Returns the constraints
	const std::vector<Linear_constraint*>& constraints() const;
	std::vector<Linear_constraint*>& constraints();

	/// Is the variable owned by this program?
	bool has_variable(const Variable* var) const;

	/// Is the constraint owned by this program?
	bool has_constraint(const Linear_constraint* cons) const;

	/// Returns the number of continuous variables
	std::size_t num_continuous_variables() const;

	/// Returns the number of integer variables
	std::size_t num_integer_variables() const;

	/// Returns the number of binary variables
	std::size_t num_binary_variables() const;

	/// Returns true if all variables are continuous
	bool is_continuous() const;

	/// Returns true if this is a mixed integer program
	bool is_mixed_integer_program() const;

	/// Returns true if this is an integer program
	bool is_integer_program() const;

	/// Returns true if binary program
	bool is_binary_program() const;

	/// Returns the objective
	const Linear_objective * objective() const;
	Linear_objective * objective();

	/// Solves the program. Returns false if fails.
	bool solve();

	/// Returns the result. 
	/// \note (1) Result is valid only if the solver succeeded.
	///       (2) Each entry in the result corresponds to the variable with the
	///			 same index in the program.
	const std::vector<FT>& solution() const;

	/// Returns the error message.
	/// \note This function should be called after call to solve().
	const std::string& error_message() const { return error_message_; }

	/// Clears all variables, constraints, and the objective.
	void clear();

	/// @}
}; /* end MixedIntegerProgramTraits */



   /*!
   \cgalConcept

   `MixedIntegerProgramTraits::Variable` is a concept of a variable in
   a Mixed Integer Programming (MIP) problem.

   \cgalHasModel `CGAL::Variable<FT>`

   */
class MixedIntegerProgramTraits::Variable
{
public:
	/// \name Types
	/// @{

	/*!
	A variable can be continuous, integer, or binary
	*/
	enum Variable_type { CONTINUOUS, INTEGER, BINARY };

	/*!

	*/
	typedef	MixedIntegerProgramTraits	Solver;


	/// @}

	/// \name Creation 
	/// @{

	/*!
	Create a variable, initialized with the solver it belongs to,
	the variable type, lower bound, upper bound, name, and index.
	*/
	Variable(Solver* solver, Variable_type type, FT lb =, FT ub, const std::string& name, int idx);

	/// \name Operations 
	/// @{

	/// Return the variable type
	Variable_type variable_type() const;

	/// Set/Change the variable type
	void set_variable_type(Variable_type t);

	/*!
	Return the name of the variable.
	*/
	const std::string& name() const;

	/*!
	Set the name of the variable.
	*/
	void set_name(const std::string& n);

	/*!
	Return the index of the variable.
	*/
	int  index() const;

	/*!
	Set the index of the variable.
	*/
	void set_index(int idx);

	/// Return the solver that owns this variable
	const Solver* solver() const;
	Solver* solver();

	/// Set the lower bound
	void set_lower_bound(FT lb);

	/// Set the upper bound
	void set_upper_bound(FT ub);

	/// Set both lower and upper bounds
	void set_bounds(FT lb, FT ub);

	/// Get the lower bound
	FT lower_bound() const;

	/// Get the upper bound
	FT upper_bound() const;

	/// Get both lower and upper bounds
	void get_bounds(FT& lb, FT& ub) const;

	/// Get the infinity threshold (e.g., 1e20). 
	/// Values greater than this value are considered as infinity.
	static FT infinity();

	/// Returns the value of the variable in the current solution.
	/// \note (1) Valid only if the program was successfully solved.
	///       (2) If the variable is integer and rounded == true, then the 
	///           value will be rounded to the nearest integer.
	FT solution_value(bool rounded = false) const;

	/// Set the solution value (should be called in the solve function of
	/// its owner solver).
	void set_solution_value(FT value);

	/// @}

}; /* end Variable */



   /*!

   \cgalConcept

   `MixedIntegerProgramTraits::Linear_constraint` is a concept of a linear
   constraint in a Mixed Integer Programming (MIP) problem.

   \cgalHasModel `CGAL::Linear_constraint<FT>`
   */
class MixedIntegerProgramTraits::Linear_constraint
{
public:
	/// \name Types
	/// @{

	/*!

	*/
	typedef	MixedIntegerProgramTraits			Solver;


	/*!

	*/
	typedef	MixedIntegerProgramTraits::Variable	Variable;


	/// @}

	/// \name Creation 
	/// @{

	/*!
	Create a linear constraint, initialized with the solver it belongs to,
	the lower bound, upper bound, name, and index.
	*/
	Linear_constraint(Solver* solver, FT lb, FT ub, const std::string& name, int idx);

	/// \name Operations 
	/// @{

	/*!
	Return the name of the constraint.
	*/
	const std::string& name() const;

	/*!
	Set the name of the constraint.
	*/
	void set_name(const std::string& n);

	/*!
	Return the index of the constraint.
	*/
	int  index() const;

	/*!
	Set the index of the constraint.
	*/
	void set_index(int idx);

	/// Return the solver that owns this constraint
	const Solver* solver() const;
	Solver* solver();

	/// Set the lower bound
	void set_lower_bound(FT lb);

	/// Set the upper bound
	void set_upper_bound(FT ub);

	/// Set both lower and upper bounds
	void set_bounds(FT lb, FT ub);

	/// Get the lower bound
	FT lower_bound() const;

	/// Get the upper bound
	FT upper_bound() const;

	/// Get both lower and upper bounds
	void get_bounds(FT& lb, FT& ub) const;

	/// Get the infinity threshold (e.g., 1e20). 
	/// Values greater than this value are considered as infinity.
	static FT infinity();

	/// Set the coefficients of the constraint. 
	void  set_coefficients(const std::unordered_map<const Variable*, FT>& coeffs);

	/// Add a coefficient to a variable of the constraint. 
	void  add_coefficient(const Variable* var, FT coeff);

	/// Return the coefficients of the constraint. 
	const std::unordered_map<const Variable*, FT>& coefficients() const;

	/// Get the coefficient of the variable in this constraint. 
	FT get_coefficient(const Variable* var) const;

	/// Set the constant term.
	void set_offset(FT value);

	/// Get the constant term.
	FT offset() const;

	/// Clear all variables, set the constant term to zero.
	/// Useful to reuse the object to define a new linear constraint.
	void clear();

	/// @}

}; /* end Linear_constraint */



   /*!

   \cgalConcept

   `MixedIntegerProgramTraits::Linear_objective` is a concept of the linear
   objective function in a Mixed Integer Programming (MIP) problem.

   \cgalHasModel `CGAL::Linear_objective<FT>`
   */
class MixedIntegerProgramTraits::Linear_objective
{
public:
	/// \name Types
	/// @{

	/// The objective sense (i.e., optimization direction)
	enum Sense { MINIMIZE, MAXIMIZE, UNDEFINED };

	/*!

	*/
	typedef	MixedIntegerProgramTraits			Solver;


	/*!

	*/
	typedef	MixedIntegerProgramTraits::Variable	Variable;


	/// @}

	/// \name Creation 
	/// @{

	/*!
	Create a linear objective, initialized with the solver it belongs to
	and the objective sense.
	*/
	Linear_objective(Solver* solver, Sense sense);

	/// \name Operations 
	/// @{

	/// Set the objective sense.
	void  set_sense(Sense sense);

	/// Get the objective sense.
	Sense sense() const;

	/// Set the coefficients of the constraint. 
	void  set_coefficients(const std::unordered_map<const Variable*, FT>& coeffs);

	/// Add a coefficient to a variable of the constraint. 
	void  add_coefficient(const Variable* var, FT coeff);

	/// Return the coefficients of the constraint. 
	const std::unordered_map<const Variable*, FT>& coefficients() const;

	/// Get the coefficient of the variable in this constraint. 
	FT get_coefficient(const Variable* var) const;

	/// Set the constant term.
	void set_offset(FT value);

	/// Get the constant term.
	FT offset() const;

	/// Clear the objective (i.e., remove all variables, reset the 
	/// objective sense to UNDEFINED). Useful to reuse the object 
	/// to define a new linear objective.
	void clear();

	/// @}

}; /* end Linear_objective */
