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
    << "[degree of the curves]\n"
    << "[number of interpolation points] \n"
    << "[percentage to use a tangent at a point] in [0..100] \n"
    << "[range for random integer points / tangent] \n" 
    << "[number of points equal to each curve ]\n" << std::endl;
}

template<class NT> NT fac(const NT& n){
  if(n<=1) return 1;
  return fac(n-1)*n;
}
template<class NT> NT binom(const NT& n,const NT& i){
  return fac(n)/(fac(i)*fac(n-i));
}


int main(int argc, char** argv) {
    if(argc<7) {
        print_help(argv[0]);
        std::exit(-1);
    }


    int curr_arg=1;
    int no_curves = atoi(argv[curr_arg++]);
    int degree = atoi(argv[curr_arg++]);
    int no_points = atoi(argv[curr_arg++]);
    int tangent_prob = atoi(argv[curr_arg++]);
    int range = atoi(argv[curr_arg++]);
    int no_equal_p = atoi(argv[curr_arg++]);

    std::cerr
      << "create bezier curves for the following parameters: \n" 
      << "[number of curves]                                   " << no_curves << "\n"
      << "[degree of the curves]                               " << degree << "\n"
      << "[number of interpolation points]                     " << no_points << "\n"
      << "[percentage to use a tangent at a point] in [0..100] " << tangent_prob << "\n"
      << "[range for random integer points / tangent]          " << range << "\n" 
      << std::endl;

    
    typedef CGAL::Arithmetic_kernel::Integer Integer;
    typedef CGAL::Arithmetic_kernel::Rational Rational;
    CGAL::Fraction_traits<Rational>::Decompose decompose;
      
    typedef CGAL::Polynomial<Rational> Poly_rat_1;
    std::vector<Integer> xpoints,ypoints;
    std::vector<Integer> xtangents,ytangents; 

    CGAL::Random random; // that is my current age .-) 
    for(int i = 0; i < no_points;i++){
      xpoints.push_back(Integer(random.get_int(-range,range)));
      xtangents.push_back(Integer(random.get_int(-range,range)));
      ypoints.push_back(Integer(random.get_int(-range,range)));
      ytangents.push_back(Integer(random.get_int(-range,range)));
    }
    std::vector<Poly_rat_1> Bersteins;   
    std::vector<Poly_rat_1> DBersteins; 
    
    CGAL::set_pretty_mode(std::cout);
    {
      Poly_rat_1 t = CGAL::shift(Poly_rat_1(1),1,0);
      for(int i = 0; i <= degree; i++){
        Bersteins.push_back(binom(degree,i)*CGAL::ipower(t-1,degree-i)*CGAL::ipower(t,i));
      }
      for(int i = 0; i <= degree; i++){
        DBersteins.push_back(CGAL::differentiate(Bersteins[i]));
      }
    }
    Rational step_size = Rational(6)/Rational(21);
    
    
    Integer lcm(1);
    {
      Rational ti = 0; 
      Integer num,den;
      for(int i = 0; i<=degree;i++){
        ti+= step_size;
        for(int j = 0; j<=degree;j++){
          decompose(Bersteins[j].evaluate(ti),num,den);
          lcm *= (den/CGAL::gcd(den,lcm));
          decompose(DBersteins[j].evaluate(ti),num,den);
          lcm *= (den/CGAL::gcd(den,lcm));
        }
      }
    }
  

    std::cout << no_curves << std::endl; 
    for(int nc = 0; nc <no_curves; nc++){
      // lets try to procude one curve:
      // program and solution types
      typedef CGAL::Quadratic_program<Integer> Program;
      typedef CGAL::Quadratic_program_solution<Rational> Solution;

      Program xlp (CGAL::EQUAL, false, 0, false, 0);  
      Program ylp (CGAL::EQUAL, false, 0, false, 0);  
      
      CGAL::Quadratic_program_options options;
      options.set_verbosity (0) ;
      
      // one equation for each interpolation point or tangent
      Rational ti = 0;      // parameter value for that point
      for (int i = 0; i<=degree; i++){ // generate degree+1 equations 
        bool use_tangent = 
          random.get_int(0,100)<tangent_prob && i < degree; 
        ti += step_size; // increase parameter value for next point; 

        std::vector<Rational> equation; 
        for (int j = 0; j<=degree; j++){ // generate degree+1 entries
          Integer num,den;
          decompose(Bersteins[j].evaluate(ti)*lcm,num,den);
          // std::cout << CGAL::integral_division(num,den) << " " ; 
          xlp.set_a(j,i,CGAL::integral_division(num,den));
          ylp.set_a(j,i,CGAL::integral_division(num,den));
          if(use_tangent){
            decompose(DBersteins[j].evaluate(ti)*lcm,num,den);
            xlp.set_a(j,i+1,CGAL::integral_division(num,den));
            ylp.set_a(j,i+1,CGAL::integral_division(num,den));
          }
          // equation.push_back(Bersteins[i].evaluate(tj));
        }
        int index = random.get_int(0,no_points-1);
        xlp.set_r(i,CGAL::EQUAL); 
        // std::cout << xpoints[index] << " " ; 
        if(i<no_equal_p){
          xlp.set_b(i,xpoints[i]*lcm);
          ylp.set_b(i,ypoints[i]*lcm);
        }else{
          xlp.set_b(i,xpoints[index]*lcm);
          ylp.set_b(i,ypoints[index]*lcm);
        }
        if(use_tangent){
          xlp.set_b(i+1,xtangents[index]*lcm);
          ylp.set_b(i+1,ytangents[index]*lcm);
          i++;
        }
        //std::cout << std::endl; 
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

      {
        Integer lcm = 1;
        // OUTPUT CURVE: 
        std::cout << xcoord.size() <<" "; // number of control points  
        std::cerr << xcoord.size() <<" "; // number of control points  
        for(int i = 0; i < xcoord.size();i++){
          CGAL::simplify(xcoord[i]);
          CGAL::simplify(ycoord[i]);
          Integer num, den, g; 
          decompose(xcoord[i],num,den);
          lcm *= den / CGAL::gcd(lcm,den);
          decompose(ycoord[i],num,den);
          lcm *= den / CGAL::gcd(lcm,den);
        }   
      
        
        for(int i = 0; i < xcoord.size();i++){
          CGAL::simplify(xcoord[i]*=lcm);
          CGAL::simplify(ycoord[i]*=lcm);
          Integer num, den;
          decompose(xcoord[i],num,den);
          std::cout << CGAL::integral_division(num,den) << " "; 
          std::cerr << CGAL::integral_division(num,den) << " ";
          decompose(ycoord[i],num,den);
          std::cerr << CGAL::integral_division(num,den)<< " "; 
          std::cout << CGAL::integral_division(num,den)<< " "; 
          
        }
        std::cout << std::endl; 
        std::cerr << std::endl; 
      }
    }
}
