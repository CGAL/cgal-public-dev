// Copyright (c) 2006-2010 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://algerbya@scm.gforge.inria.fr/svn/cgal/branches/features/Algebraic_kernel_d-RS_bivariate-nancy/Algebraic_kernel_d/include/CGAL/RS/Algebraic_kernel_rs_1.h $
// $Id: Algebraic_kernel_rs_1.h 63060 2011-04-20 12:13:38Z penarand $
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>
// Author: Yacine Bouzidi <Bouzidi.yacine@gmail.com>

#ifndef CGAL_RS_ALGEBRAIC_KERNEL_RS_1
#define CGAL_RS_ALGEBRAIC_KERNEL_RS_1



#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/RS/New_algebraic_1.h> 
#include <CGAL/RS/New_refine_rs_1.h> 
#include <CGAL/RS/rs_calls_1.h>
#include <CGAL/RS/New_sign_1.h>

#ifdef IEEE_DBL_MANT_DIG
#  define CGAL_RS_FUNCTORS_DBL_PREC IEEE_DBL_MANT_DIG
#else
#  define CGAL_RS_FUNCTORS_DBL_PREC 53
#endif

namespace CGAL {

template <class _C>
struct New_Algebraic_kernel_rs_1{

  
  typedef _C                                      Coefficient;
  typedef CGAL::Polynomial<Coefficient>           Polynomial_1;
  typedef CGAL::NewAlg::New_Algebraic_1<Coefficient>                   Algebraic_real_1;
  typedef CGAL::Gmpfr                       Bound;
  typedef int                Multiplicity;
  typedef CGAL::Polynomial_traits_d< Polynomial_1 >   PT_1;
  
        // constructor: we must initialize RS just a time, so this is a
        // good time to do it
  New_Algebraic_kernel_rs_1(){CGAL::init_solver();};
  ~New_Algebraic_kernel_rs_1(){CGAL::reset_solver();};
  
  
  
  
  
  
  struct Compute_polynomial_1 :
    public std::unary_function< Algebraic_real_1, Polynomial_1 > {
    Polynomial_1 operator()(const Algebraic_real_1 &a)const{
    return a.pol();
  }
};      // Compute_polynomial_1



struct Is_square_free_1 
  : public std::unary_function< Polynomial_1, bool > {
  bool operator()( const Polynomial_1& p ) const {
    typename CGAL::Polynomial_traits_d< Polynomial_1 >::Is_square_free isf;
    return isf(p);
  }
}; // Is_square_free


struct Make_square_free_1
  : public std::unary_function< Polynomial_1, Polynomial_1 > {
  Polynomial_1 operator()( const Polynomial_1& p ) const {
    return typename CGAL::Polynomial_traits_d< Polynomial_1 >::Make_square_free()( p );
  }
}; //Make_square_free



struct Square_free_factorize_1 {
    template< class OutputIterator>
    OutputIterator operator()( const Polynomial_1& p, OutputIterator it) const {
      typename PT_1::Square_free_factorize_up_to_constant_factor sqff;
      return sqff(p,it);
    } 
}; //Square_free_factorize_1



struct Is_coprime_1
    : public std::binary_function< Polynomial_1, Polynomial_1, bool > {
    bool operator()( const Polynomial_1& p1, const Polynomial_1& p2 ) const {
      typename CGAL::Polynomial_traits_d< Polynomial_1 >::Total_degree total_degree;
                        
      // TODO: Is GCD already filtered? 
      return( total_degree( gcd_utcf( p1, p2 ) ) == 0 );                        
    } 
}; //Is_coprime_1




struct Make_coprime_1 {
    typedef bool         result_type;
    typedef Polynomial_1 first_argument_type;
    typedef Polynomial_1 second_argument_type;
    typedef Polynomial_1 third_argument_type;
    typedef Polynomial_1 fourth_argument_type;
    typedef Polynomial_1 fifth_argument_type;
                
    bool operator()( const Polynomial_1& p1,
		     const Polynomial_1& p2,
        Polynomial_1& g, // ggT utcf 
        Polynomial_1& q1, // Rest utcf
        Polynomial_1& q2 ) const {
      g = typename CGAL::Polynomial_traits_d< Polynomial_1 >::Gcd_up_to_constant_factor()( p1, p2 );
      q1 = p1 / g;
      q2 = p2 / g;
      return CGAL::is_one(g);
    }                 
}; // Make_coprime_1





  
struct Solve_1{
  
  
  
  typedef typename CGAL::Polynomial_traits_d< Polynomial_1 >::Make_square_free Make_square_free;
  typedef std::vector< std::pair<Polynomial_1,int> >                      sqfrvec;
  typedef typename PT_1::Square_free_factorize_up_to_constant_factor sqff;
  
  template <class OutputIterator>
  OutputIterator operator()(const Polynomial_1 &p,OutputIterator res)const{
    
    int nr,*m;
    mpfi_ptr *x;
    sqfrvec sfv;
    sqff()(p,back_inserter(sfv));
    Polynomial_1 sfp = Make_square_free()( p );
    x=(mpfi_ptr*)malloc(sfp.degree()*sizeof(mpfi_ptr));
    m=(int*)calloc(sfp.degree(),sizeof(int));
    nr=solve_1(x,sfp);
      CGAL_assertion_msg(nr>=0,"error in resolution");
      for(int i=0;i<sfv.size();++i){
            int k=sfv[i].first.degree();
            for(int j=0;k&&j<nr;++j){
                if(!m[j]){
                    CGAL::Sign sg_l=
		      CGAL::RSSign::signat(sfv[i].first,&(x[j]->left));
                    CGAL::Sign sg_r=
		      CGAL::RSSign::signat(sfv[i].first,&(x[j]->right));
                    if((sg_l!=sg_r)||((sg_l==CGAL::ZERO)&&(sg_r==CGAL::ZERO))){
                        m[j]=sfv[i].second;
                        --k;
                    }
                }
            }
        }
        for(int i=0;i<nr;++i){
	  *res++=std::make_pair(*new Algebraic_real_1(x[i],p,i,m[i],
                                                 i?x[i-1]:NULL,
                                                 i==nr-1?NULL:x[i+1]),
                                  m[i]);
        }
        free(m);
        free(x);
        return res;
    }

    template <class OutputIterator>
    OutputIterator operator()(const Polynomial_1 &p,
                              bool known_to_be_square_free,
                              OutputIterator res)const{
        int nr,m;
        mpfi_ptr *x;
        if(known_to_be_square_free){
	  x=(mpfi_ptr*)malloc(p.degree()*sizeof(mpfi_ptr));
	  nr=solve_1(x,p);
	  CGAL_assertion_msg(nr>=0,"error in resolution");
	  m=1;    // we know in this case that multiplicity is 1
        }else{
	  x=(mpfi_ptr*)malloc(Make_square_free()(p).degree()*sizeof(mpfi_ptr));
	  nr=solve_1(x,Make_square_free()(p));
            CGAL_assertion_msg(nr>=0,"error in resolution");
            m=0;    // we won't calculate multiplicities
        }
        for(int i=0;i<nr;++i)
	  *res++=*new Algebraic_real_1(x[i],p,i,m,
                            i?x[i-1]:NULL,
                            i==nr-1?NULL:x[i+1]);
        free(x);
        return res;
    }
};  // Solve_1

/* template <class _P,class _Gcd_policy> */
/* struct Solve_1{ */
/*     typedef _P                  P; */
/*     typedef CGAL::to_rs_poly<P> convert; */
/*     typedef _Gcd_policy         Gcd; */
/*     typedef Solve_RS_1<Gcd>     Solve_RS; */

/*     template <class OutputIterator> */
/*     OutputIterator operator()(const P &p,OutputIterator res)const{ */
/*         return Solve_RS()(convert()(p),res); */
/*     } */

/*     template <class OutputIterator> */
/*     OutputIterator operator()(const P &p, */
/*                               bool known_to_be_square_free, */
/*                               OutputIterator res)const{ */
/*         return Solve_RS()(convert()(p),known_to_be_square_free,res); */
/*     } */

/*   template <class OutputIterator> */
/*   OutputIterator operator()(const P &p, */
/*                             const Bound& lower, */
/*                             const Bound& upper, */
/*                             OutputIterator res)const{ */
/*     typedef std::vector<std::pair<Algebraic,Multiplicity> > RMV; */
/*     RMV roots; */
/*     this->operator()(p,std::back_inserter(roots)); */
/*     for(typename RMV::iterator it = roots.begin(); it != roots.end();it++){ */
/*       if(lower <= it->first && it->first <= upper) */
/*         *res++=*it; */
/*     } */
/*     return res; */
/*   } */

/*   template <class OutputIterator> */
/*   OutputIterator operator()(const P &p, */
/*                             bool known_to_be_square_free, */
/*                             const Bound& lower, */
/*                             const Bound& upper, */
/*                             OutputIterator res)const{ */
/*     typedef std::vector< Algebraic > RV; */
/*     RV roots; */
/*     this->operator()(p,known_to_be_square_free,std::back_inserter(roots)); */
/*     for(typename RV::iterator it = roots.begin(); it != roots.end();it++){ */
/*       if(lower <= *it && *it <= upper) */
/*         *res++=*it; */
/*     } */
/*     return res; */
/*   } */
/* };  // Solve_1 */




 struct Construct_algebraic_real_1{
   

   Algebraic_real_1 operator()(int a) const {
     return Algebraic_real_1(a);
    }

   Algebraic_real_1 operator()(const Bound a) const {
     return Algebraic_real_1(a);
    }

   Algebraic_real_1 operator()(const Coefficient a) const {
     return Algebraic_real_1(a);
    }

   Algebraic_real_1 operator()(const Polynomial_1 &p,int i) const {
      CGAL_precondition(CGAL::is_square_free(p));
      std::vector<Algebraic_real_1> roots;
      std::back_insert_iterator<std::vector<Algebraic_real_1> > rootsit(roots);
      Solve_1()(p,true,rootsit);
      return roots[i];
    }
   
   Algebraic_real_1 operator()(const Polynomial_1 &p,Bound l,Bound u) const {
        mpfi_t i;
        mpfi_init(i);
        mpfr_set(&i->left,l.fr(),GMP_RNDD);
        mpfr_set(&i->right,u.fr(),GMP_RNDU);
        return Algebraic_real_1(i,p,0,0,NULL,NULL);
    }
};  // Construct_algebraic_real_1


  

struct Number_of_solutions_1:
    public std::unary_function<Polynomial_1,int>{
  
  int operator()(const Polynomial_1 &p)const{
    int nr;
    mpfi_ptr *x;
    //CGAL::RS_polynomial_1 rspoly=convert()(p);
    x=(mpfi_ptr*)malloc(p.degree()*sizeof(mpfi_ptr));
    nr=solve_1(x,p);
    CGAL_assertion_msg(nr>=0,"error in resolution");
    free(x);
    return nr;
  }
};  // Number_of_solutions_1

  
struct Sign_at_1:
    public std::binary_function<Polynomial_1,Algebraic_real_1,CGAL::Sign>{
  
  CGAL::Sign operator()(const Polynomial_1 &p,const Algebraic_real_1 &a)const{
    return RS3::sign_1(p,a);
    }
};  // Sign_at_1

  
struct Is_zero_at_1:
    public std::binary_function<Polynomial_1,Algebraic_real_1,bool>{
  
  bool operator()(const Polynomial_1 &p,const Algebraic_real_1 &a)const{
    return (Sign_at_1()(p,a)==CGAL::ZERO);
  }
};  // Is_zero_at_1
  

struct Compare_1:
  public std::binary_function<Algebraic_real_1,Algebraic_real_1,CGAL::Comparison_result>{
  
  typedef CGAL::Comparison_result       Comparison_result;
  typedef CGAL::Gmpz                    Gmpz;
  typedef CGAL::Gmpq                    Gmpq;

  Comparison_result operator()(const Algebraic_real_1 &r1,const Algebraic_real_1 &r2)const{
    return CGAL::NewAlg::RS_COMPARE::compare_1<Coefficient>(r1,r2);
  }
  
  Comparison_result operator()(const int &r1,  const Algebraic_real_1 &r2)const{
    return this->operator()(Algebraic_real_1(r1),r2);}
  Comparison_result operator()(const Bound &r1,const Algebraic_real_1 &r2)const{
    return this->operator()(Algebraic_real_1(r1),r2);}
  Comparison_result operator()(const Gmpz &r1, const Algebraic_real_1 &r2)const{
    return this->operator()(Algebraic_real_1(r1),r2);}
  Comparison_result operator()(const Gmpq &r1, const Algebraic_real_1 &r2)const{
    return this->operator()(Algebraic_real_1(r1),r2);}
  Comparison_result operator()(const Algebraic_real_1 &r1,const int   &r2)const{
    return this->operator()(r1,Algebraic_real_1(r2));}
  Comparison_result operator()(const Algebraic_real_1 &r1,const Bound &r2)const{
    return this->operator()(r1,Algebraic_real_1(r2));}
  Comparison_result operator()(const Algebraic_real_1 &r1,const Gmpz  &r2)const{
    return this->operator()(r1,Algebraic_real_1(r2));}
  Comparison_result operator()(const Algebraic_real_1 &r1,const Gmpq  &r2)const{
    return this->operator()(r1,Algebraic_real_1(r2));}
};  // Compare_1

  
struct Isolate_1:
    public std::binary_function<Algebraic_real_1,Polynomial_1,std::pair<Bound,Bound> >{
  
  std::pair<Bound,Bound> operator()(const Algebraic_real_1 &a,const Polynomial_1 &p)const{
    std::vector<Algebraic_real_1> roots;
    std::back_insert_iterator<std::vector<Algebraic_real_1> > rootsit(roots);
    Solve_1()(p,true,rootsit);
    for(/*std::vector<Algebraic_real_1>::size_type*/ int i=0;i<roots.size();++i)
      Compare_1()(a,roots[i]);
    return std::make_pair(Bound(a.left()),Bound(a.right()));
  }
};  // Isolate_1


struct Bound_between_1:
  public std::binary_function<Algebraic_real_1,Algebraic_real_1,Bound>{
  
  Bound operator()(const Algebraic_real_1 &x1,const Algebraic_real_1 &x2)const{
            double l,r,m;
            switch(CGAL::NewAlg::RS_COMPARE::compare_1<Coefficient>(x1,x2)){
                case CGAL::LARGER:
                    CGAL_assertion(x2.sup()<x1.inf());
                    l=x2.sup().to_double(std::round_toward_infinity);
                    r=x1.inf().to_double(std::round_toward_neg_infinity);
                    m=(l+r)/2;
                    if(l<m&&m<r){
		      return Bound(m,CGAL_RS_FUNCTORS_DBL_PREC);
                    }
                    return Bound::add(x2.sup(),
                                      x1.inf(),
                                      (x2.sup().get_precision()>
                                                x1.inf().get_precision()?
                                       1+x2.sup().get_precision():
                                       1+x1.inf().get_precision()))/2;
                    break;
                case CGAL::SMALLER:
                    CGAL_assertion(x1.sup()<x2.inf());
                    l=x1.sup().to_double(std::round_toward_infinity);
                    r=x2.inf().to_double(std::round_toward_neg_infinity);
                    m=(l+r)/2;
                    if(l<m&&m<r){
                        return Bound(m,CGAL_RS_FUNCTORS_DBL_PREC);
                    }
                    return Bound::add(x1.sup(),
                                      x2.inf(),
                                      (x1.sup().get_precision()>
                                                x2.inf().get_precision()?
                                       1+x1.sup().get_precision():
                                       1+x2.inf().get_precision()))/2;
                    break;
                default:
                    CGAL_error_msg("bound between two equal numbers");
            }
        }
    };  // Bound_between_1

struct Approximate_absolute_1:
    public std::binary_function<Algebraic_real_1,int,std::pair<Bound,Bound> >{
  
  std::pair<Bound,Bound>
    operator()(const Algebraic_real_1& x, int prec) const {
    //--------------------------------------------------
    //     Bound error = CGAL::ipower(Bound(2),CGAL::abs(prec));
    //     while(prec>0?
//           (x.sup()-x.inf())*error>Bound(1):
//           (x.sup()-x.inf())>error){
//       RS3::refine_1(x,CGAL::abs(prec));
//     }
//--------------------------------------------------
    CGAL::NewAlg::refine_1<Coefficient>(x,std::max<unsigned>(CGAL::abs(prec),
                                       mpfi_get_prec(x.mpfi())));
    //if(mpfi_get_prec(x.mpfi())<CGAL::abs(prec)){
    //        RS3::refine_1(x,CGAL::abs(prec));
    //}
    CGAL_assertion(prec>0?
                   (x.sup()-x.inf())*CGAL::ipower(Bound(2),prec)<=Bound(1):
                   (x.sup()-x.inf())<=CGAL::ipower(Bound(2),-prec));
    return std::make_pair(x.inf(),x.sup());
  }
};

struct Approximate_relative_1
  :public std::binary_function<Algebraic_real_1,int,std::pair<Bound,Bound> >{

  std::pair<Bound,Bound> operator()(const Algebraic_real_1 &x, int prec) const {
    if(CGAL::is_zero(x))
      return std::make_pair(Bound(0),Bound(0));

    Bound error = CGAL::ipower(Bound(2),CGAL::abs(prec));
    Bound max_b = (CGAL::max)(CGAL::abs(x.sup()),CGAL::abs(x.inf()));
    while(prec>0?
          (x.sup()-x.inf())*error>max_b:
          (x.sup()-x.inf())>error*max_b){
      CGAL::NewAlg::refine_1<Coefficient>(x,std::max<unsigned>(CGAL::abs(prec),
                                         mpfi_get_prec(x.mpfi())));
      max_b = (CGAL::max)(CGAL::abs(x.sup()),CGAL::abs(x.inf()));
    }
    CGAL_assertion(prec>0?
                   (x.sup()-x.inf())*CGAL::ipower(Bound(2),prec)<=max_b:
                   (x.sup()-x.inf())<=CGAL::ipower(Bound(2),-prec)*max_b);
    return std::make_pair(x.inf(),x.sup());
  }
};


#define CGALRS_CREATE_FUNCTION_OBJECT(T,N)	\
  T N##_object()const{return T();} 
  CGALRS_CREATE_FUNCTION_OBJECT(Construct_algebraic_real_1, 
                                         construct_algebraic_real_1) 
  CGALRS_CREATE_FUNCTION_OBJECT(Compute_polynomial_1,compute_polynomial_1) 
  CGALRS_CREATE_FUNCTION_OBJECT(Isolate_1,isolate_1) 
	   CGALRS_CREATE_FUNCTION_OBJECT(Is_square_free_1,is_square_free_1) 
	   CGALRS_CREATE_FUNCTION_OBJECT(Make_square_free_1,make_square_free_1)
	   CGALRS_CREATE_FUNCTION_OBJECT(Square_free_factorize_1, 
				square_free_factorize_1) 
	   CGALRS_CREATE_FUNCTION_OBJECT(Is_coprime_1,is_coprime_1) 
	   CGALRS_CREATE_FUNCTION_OBJECT(Make_coprime_1,make_coprime_1) 
	   CGALRS_CREATE_FUNCTION_OBJECT(Solve_1,solve_1)
	   CGALRS_CREATE_FUNCTION_OBJECT(Number_of_solutions_1, 
					 number_of_solutions_1) 
	   CGALRS_CREATE_FUNCTION_OBJECT(Sign_at_1,sign_at_1) 
           CGALRS_CREATE_FUNCTION_OBJECT(Is_zero_at_1,is_zero_at_1) 
	   CGALRS_CREATE_FUNCTION_OBJECT(Compare_1,compare_1)
	   CGALRS_CREATE_FUNCTION_OBJECT(Bound_between_1,bound_between_1) 
	   CGALRS_CREATE_FUNCTION_OBJECT(Approximate_absolute_1, 
					 approximate_absolute_1) 
	   CGALRS_CREATE_FUNCTION_OBJECT(Approximate_relative_1, 
					 approximate_relative_1) 
#undef CGALRS_CREATE_FUNCTION_OBJECT 
};  // Algebraic_kernel_d_1_RS

} // namespace CGAL
#endif  // CGAL_RS_ALGEBRAIC_KERNEL_RS_1
