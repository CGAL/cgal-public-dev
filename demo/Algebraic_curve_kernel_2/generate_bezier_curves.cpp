#include <CGAL/basic.h>
#include <sstream>
#include <CGAL/Polynomial.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Random.h>
#include <CGAL/ipower.h>

#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <cassert>



void print_help(char* execname) {
  std::cout  
    << "This program creates a set of semi random bezier curves \n"
    << "It requires the following input parameters: \n" 
    << "[number of curves] \n"
    << "[degree of the curves]\n "
    << "[number of interpolation points] \n"
    // << "[percentage to use a tangent at a point] in [0..100] \n"
    << "[range for random integer points / tangent] \n" << std::endl;
}

int main(int argc, char** argv) {
    if(argc<3) {
        print_help(argv[0]);
        std::exit(-1);
    }


    int curr_arg=1;
    int no_curves = atoi(argv[curr_arg++]);
    int degree = atoi(argv[curr_arg++]);
    int no_points = atoi(argv[curr_arg++]);
    // int tangent_prob = atoi(argv[curr_arg++]);
    int range = atoi(argv[curr_arg++]);

    std::cerr
      << "create bezier curves for the following parameters: \n" 
      << "[number of curves]                                   " << no_curves << "\n"
      << "[degree of the curves]                               " << degree << "\n"
      << "[number of interpolation points]                     " << no_points << "\n"
      // << "[percentage to use a tangent at a point] in [0..100] " << tangent_prob << "  \n"
      << "[range for random integer points / tangent]          " << range << " \n" 
      << std::endl;

    
    typedef CGAL::Arithmetic_kernel::Integer Integer;
    typedef CGAL::Arithmetic_kernel::Rational Rational;
    typedef CGAL::Polynomial<Rational> Poly_rat_1;
    std::vector<Rational> xpoints,ypoints;
    std::vector<Rational> xtangents,ytangents; 

    CGAL::Random random(31); // that is my current age .-) 
    for(int i = 0; i < no_points;i++){
      xpoints.push_back(Rational(random.get_int(-range,range)));
      xtangents.push_back(Rational(random.get_int(-range,range)));
      ypoints.push_back(Rational(random.get_int(-range,range)));
      ytangents.push_back(Rational(random.get_int(-range,range)));
    }
    std::vector<Poly_rat_1> Bersteins;   
    std::vector<Poly_rat_1> DBersteins; 
    {
      Poly_rat_1 t = CGAL::shift(Poly_rat_1(1),1,0);
      for(int i = 0; i <= degree; i++){
        Bersteins.push_back(CGAL::ipower(t-1,degree-i)*CGAL::ipower(t,i));
      }
      for(int i = 0; i <= degree; i++){
        DBersteins.push_back(CGAL::differentiate(Bersteins[i]));
      }
    }
    

    std::cout << no_curves << std::endl; 
    for(int nc = 0; nc <no_curves; nc++){
      // lets try to procude one curve:
      // program and solution types
      typedef CGAL::Quadratic_program<Rational> Program;
      typedef CGAL::Quadratic_program_solution<Rational> Solution;

      Program xlp (CGAL::EQUAL, false, 0, false, 0);  
      Program ylp (CGAL::EQUAL, false, 0, false, 0);  
      
      CGAL::Quadratic_program_options options;
      options.set_verbosity (0) ;

      for (int i = 0; i<=degree; i++){ // generate degree+1 equations 
        Rational t0 = Rational(1)/(degree+2);
        Rational tj = 0; 
        std::vector<Rational> equation; 
        for (int j = 0; j<=degree; j++){ // generate degree+1 entries 
          tj += t0;
          xlp.set_a(j,i,Bersteins[i].evaluate(tj));
          ylp.set_a(j,i,Bersteins[i].evaluate(tj));
          // equation.push_back(Bersteins[i].evaluate(tj));
        }
        int index = random.get_int(0,no_points-1);
        xlp.set_r(i,CGAL::EQUAL); 
        xlp.set_b(i,xpoints[index]);
        ylp.set_b(i,ypoints[index]);   
      }
      
      Solution xs = CGAL::solve_linear_program(xlp, Rational(),options);
      Solution ys = CGAL::solve_linear_program(ylp, Rational(),options);
      
      std::vector<Rational> xcoord, ycoord; 
      
      for(Solution::Variable_value_iterator it = xs.variable_values_begin (); 
          it != xs.variable_values_end(); it++){ 
        xcoord.push_back(it->numerator() / it->denominator() );
      }
      for(Solution::Variable_value_iterator it = ys.variable_values_begin (); 
          it != ys.variable_values_end(); it++){ 
        ycoord.push_back(it->numerator() / it->denominator() );
      }

      // OUTPUT CURVE: 
      std::cout << xcoord.size() <<" "; // number of control points  
      for(int i = 0; i < xcoord.size();i++){
        CGAL::simplify(xcoord[i]);
        CGAL::simplify(ycoord[i]);
        std::cout << xcoord[i] << " " << ycoord[i] << " "; // points  
      }
      std::cout << std::endl; 
    }
}
