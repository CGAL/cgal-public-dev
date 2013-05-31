//#define CGAL_ACK_DEBUG_FLAG 1

#include <CGAL/basic.h>

#ifndef CGAL_USE_GMP
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs GMP ..." << std::endl; 
  return 0;
}
#else

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Timer.h>

#include <string>
#include <fstream>

typedef CGAL::GMP_arithmetic_kernel::Integer Integer;
typedef CGAL::Arr_algebraic_segment_traits_2<Integer> Arr_traits_2;
typedef Arr_traits_2::CKvA_2 Curved_kernel_via_analysis_2;
typedef Arr_traits_2::Algebraic_kernel_d_2 Algebraic_kernel_d_2;
typedef Arr_traits_2::Algebraic_kernel_d_1 Algebraic_kernel_d_1;
typedef CGAL::Arrangement_2<Arr_traits_2> Arrangement_2;
typedef Arr_traits_2::Curve_2 Curve_2;
typedef Arr_traits_2::Polynomial_2 Polynomial_2;
typedef Arr_traits_2::X_monotone_curve_2 Arc_2;
typedef Arr_traits_2::Point_2 Point_2;

typedef Algebraic_kernel_d_2::Bound Rational;
typedef Algebraic_kernel_d_2::Algebraic_real_1 Algebraic_real_1;
typedef Algebraic_kernel_d_2::Algebraic_real_2 Algebraic_real_2;


// Traits class for Polynomial type
typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;
typedef Polynomial_traits_2::Coefficient_type Polynomial_1;

typedef CGAL::Polynomial_type_generator<Integer,3>::Type Polynomial_3;



Polynomial_1 t = CGAL::shift(Polynomial_1(1),1);
Polynomial_3 x = CGAL::shift(Polynomial_3(1),1,0);
Polynomial_3 y = CGAL::shift(Polynomial_3(1),1,1);

template<typename OutputIterator>
int is_in_zero_one_param(Polynomial_1 X,
			 Polynomial_1 Y,
			 Curve_2 C,
			 Rational x,
			 OutputIterator out) {
  typedef CGAL::Fraction_traits<Rational> FT;
  FT::Decompose decompose;
  Integer xnum, xdenom;
  decompose(x,xnum,xdenom);
  std::vector<Algebraic_real_1> t_vec;
  Algebraic_kernel_d_1::Solve_1()(X*xdenom-xnum,false,Rational(0),Rational(1),
				std::back_inserter(t_vec));
  Curve_2::Status_line_1 status_line=C.status_line_at_exact_x(x);
  int no_arcs=status_line.number_of_events();
  
  for(int i=0; i<(int)t_vec.size();i++) {
    Algebraic_real_1 t = t_vec[i];
    CGAL::internal::Interval_evaluate_1<Polynomial_1,Rational> 
      interval_eval;
    int stack_index=0;
    int prec=4;
    
    while(true) {
      std::pair<Rational,Rational> curr_approx = interval_eval
	(Y, 
	 Algebraic_kernel_d_1::Approximate_absolute_1()(t,prec));
      while(curr_approx.first > status_line.upper_bound(stack_index)) {
	stack_index++;
      }
      if(curr_approx.second < status_line.upper_bound(stack_index)) {
	break;
      }
      prec*=2;
    }
    *out++=stack_index;
  }

  return (int)t_vec.size();

}

Integer binom(int n, int k) {
  if(n/2<k) {
    return binom(n,n-k);
  }
  Integer res(1);
  for(int j=n-k+1;j<=n;j++) {
    res*=Integer(j);
  }
  for(int j=2;j<=k;j++) {
    res/=Integer(j);
  }
  return res;
}

template<typename OutputStream>
std::pair<Polynomial_1,Polynomial_1>
get_param(OutputStream& str, Integer& start_x, Integer& end_x,
	  Integer& start_y, Integer& end_y) {

  int n;
  str >> n;
  n--;
  Polynomial_1 omt=-t+Polynomial_1(1);
  Polynomial_1 res_x=Polynomial_1(0),res_y=Polynomial_1(0);
  for(int k=0;k<=n;k++) {
    Integer pkx, pky;
    str >> pkx >> pky;
    if(k==0) {
      start_x=pkx;
      start_y=pky;
    } 
    if(k==n) {
      end_x=pkx;
      end_y=pky;
    }
    Integer curr_binom=binom(n,k);
    res_x+=CGAL::ipower(t,k)*CGAL::ipower(omt,n-k)*pkx*curr_binom;
    res_y+=CGAL::ipower(t,k)*CGAL::ipower(omt,n-k)*pky*curr_binom;
  }
  /*
  std::cout << "HERE: " << std::endl;
  std::cout << res_x << std::endl;
  std::cout << res_y << std::endl;
  */
  return std::make_pair(res_x,res_y);
}

int main(int argc, char** argv) {

  if(argc<2) {
    std::cerr << "Need filename" << std::endl;
    std::exit(0);
  }

  // For nice printouts
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);

  CGAL::Timer impl_timer, ca_timer,select_timer, arr_timer;
  
  Algebraic_kernel_d_2 ak_2 =Curved_kernel_via_analysis_2::instance().kernel();
  Arr_traits_2 arr_traits;

  Arr_traits_2::Construct_curve_2 construct_curve
    = arr_traits.construct_curve_2_object();

  Arrangement_2 arr(&arr_traits);

  std::ifstream ifstr(argv[1]);
  int no_curves;

  ifstr >> no_curves;

  std::cout << "#Curves: " << no_curves << std::endl;

  for(int i=0;i<no_curves;i++) {

    impl_timer.start();

    Integer start_x,end_x,start_y,end_y;

    std::pair<Polynomial_1,Polynomial_1> poly_pair 
      = get_param(ifstr,start_x,end_x,start_y,end_y);
    Polynomial_3 first_pol=Polynomial_3(Polynomial_2(poly_pair.first));
    first_pol=x-CGAL::Polynomial_traits_d<Polynomial_3>::Move()(first_pol,0,2);
    Polynomial_3 second_pol=Polynomial_3(Polynomial_2(poly_pair.second));
    second_pol=y-CGAL::Polynomial_traits_d<Polynomial_3>::Move()
      (second_pol,0,2);
    Polynomial_2 curr_pol 
      = CGAL::make_square_free(CGAL::resultant(first_pol,second_pol));
    std::cerr << curr_pol << std::endl;
    
    impl_timer.stop();

    ca_timer.start();

    Curve_2 cv = construct_curve(curr_pol);
    int n = cv.number_of_status_lines_with_event();
    
    std::vector<Algebraic_real_1> crit_values;

    for(int i=0;i<n;i++) {
      crit_values.push_back(cv.status_line_at_event(i).x());
    }
    ca_timer.stop();
    select_timer.start();
    crit_values.push_back(Algebraic_real_1(start_x));
    crit_values.push_back(Algebraic_real_1(end_x));
    std::sort(crit_values.begin(),crit_values.end());
    std::vector<Algebraic_real_1>::iterator end_uniq 
      = std::unique(crit_values.begin(),crit_values.end());
    std::vector<Rational> interm_vals;
    Algebraic_kernel_d_1 dummy;
    CGAL::internal::find_intermediate_values
      (&dummy,crit_values.begin(),
       end_uniq,std::back_inserter(interm_vals));

    std::vector<std::pair<Rational,int> > arc_info;

    for(int i=0; i<(int)interm_vals.size(); i++) {
      Rational x = interm_vals[i];
      std::cerr << "x=" << CGAL::to_double(x) << std::endl;
      std::vector<int> arcnos;
      int no = is_in_zero_one_param(poly_pair.first,poly_pair.second,cv,x,
				    std::back_inserter(arcnos));
      for(int j=0;j<no;j++) {
	arc_info.push_back(std::make_pair(x,arcnos[j]));
      }

      std::cerr << "found " << no << " arcs on the bezier curve" << std::endl;
    }

    std::vector<CGAL::Object> all_objects;
    std::vector<Arc_2> all_arcs;
    std::vector<Arc_2> sel_arcs;
    arr_traits.make_x_monotone_2_object()(cv,std::back_inserter(all_objects));

    for(int i=0;i<(int)all_objects.size();i++) {
      Arc_2 arc;
      if(CGAL::assign(arc,all_objects[i])) {
	all_arcs.push_back(arc);
      } else {
	//std::cout << "Ignore isolating point" << std::endl;
      }
    }


    Algebraic_kernel_d_2::Construct_algebraic_real_2 construct_alg_real_2
      = ak_2.construct_algebraic_real_2_object();
    Algebraic_real_2 start_point_ar2=construct_alg_real_2(start_x,start_y);
    Point_2 start_point(start_point_ar2.x(),start_point_ar2.curve(),
			start_point_ar2.arcno());
    Algebraic_real_2 end_point_ar2=construct_alg_real_2(end_x,end_y);
    Point_2 end_point(end_point_ar2.x(),end_point_ar2.curve(),
			end_point_ar2.arcno());

    for(std::vector<Arc_2>::iterator it=all_arcs.begin();
	it!=all_arcs.end();it++) {
      if(arr_traits.is_on_2_object()(start_point,*it)) {
	if(! it->is_in_x_range_interior(Algebraic_real_1(start_x))) {
	  break;
	}
	//std::cout << "Found start" << std::endl;
	Arc_2 left, right;
	arr_traits.split_2_object()(*it,start_point,left,right);
	it=all_arcs.erase(it);
	all_arcs.insert(it,left);
	all_arcs.insert(it,right);
	break;
      }
    }
    for(std::vector<Arc_2>::iterator it=all_arcs.begin();
	it!=all_arcs.end();it++) {
      if(arr_traits.is_on_2_object()(end_point,*it)) {
	if(! it->is_in_x_range_interior(Algebraic_real_1(end_x))) {
	  break;
	}
	//std::cout << "Found end" << std::endl;
	Arc_2 left, right;
	arr_traits.split_2_object()(*it,end_point,left,right);
	it=all_arcs.erase(it);
	all_arcs.insert(it,left);
	all_arcs.insert(it,right);
	break;
      }
    }

    std::cerr << "select arcs" << std::endl;
    for(int i=0;i<(int)all_arcs.size();i++) {
      Arc_2 arc = all_arcs[i];
      if((! arc.is_finite(CGAL::ARR_MIN_END)) || 
	 (! arc.is_finite(CGAL::ARR_MAX_END))) {
	continue;
      }
      std::cerr << "curr arc: " << arc << std::endl;
      for(int k=0; k<(int)arc_info.size();k++) {
	Rational x = arc_info[k].first;
	int arcno = arc_info[k].second;
	if(arc.is_in_x_range_interior(Algebraic_real_1(x)) && 
	   arc.arcno()==arcno) {
	  std::cerr << "Select arc" << std::endl;
	  sel_arcs.push_back(arc);
	  k=(int)arc_info.size();
	}
      }
    }

    select_timer.stop();

    arr_timer.start();

    //CGAL::insert(arr,construct_curve(curr_pol));
    CGAL::insert(arr,sel_arcs.begin(),sel_arcs.end());
    std::cerr << "done" << std::endl;
    arr_timer.stop();
  }

  // Print the arrangement size.
    std::cout << "The arrangement size:" << std::endl
              << "   V = " << arr.number_of_vertices()
              << ",  E = " << arr.number_of_edges() 
              << ",  F = " << arr.number_of_faces() 
	      << ", " << arr.number_of_isolated_vertices() << std::endl;
    
    std::cout << "Implicitization: " << impl_timer.time() << "\n"
	      << "Curve analysis:  " << ca_timer.time() << "\n"
	      << "Selection:       " << select_timer.time() << "\n"
	      << "Arrangement:     " << arr_timer.time() << std::endl;


  return 0;
}

#endif
