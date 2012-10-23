
#ifndef CGAL_RESIDUE_TYPE_H
#define CGAL_RESIDUE_TYPE_H

#include <vector>
#include <cfloat>

class Residue;
    
Residue operator + (const Residue&);
Residue operator - (const Residue&);
Residue mod_inverse(const Residue&);

std::ostream& operator << (std::ostream& os, const Residue& p);
std::istream& operator >> (std::istream& is, Residue& p);

/*! \ingroup CGAL_Modular_traits
 * \brief This class represents the Field Z mod p. 
 *  
 * This class uses the type double for representation. 
 * Therefore the value of p is restricted to primes less than 2^26.
 * By default p is set to 67111067.
 *
 * It provides the standard operators +,-,*,/ as well as in&output.
 * 
 * \see Modular_traits
 */


class Residue {
    
public:
    typedef Residue Self;
    typedef Residue NT;

//private:
    static const double  CST_CUT; 

  static int modulus;
  static double prime;
  static double prime_inv;
  static int get_prime_int(){ return modulus;}
  static double get_prime()    { return prime;}
  static double get_prime_inv(){ return prime_inv;}  

    /* Quick integer rounding, valid if a<2^51. for double */ 
    static inline 
    double RES_round (double a){
      // call CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST)
      // before using modular arithmetic 
      return ( (a + CST_CUT)  - CST_CUT);      
    }

    /* Big modular reduction (e.g. after multiplication) */
    static inline 
    double RES_reduce (double a){
      double result = a - get_prime() * RES_round(a * get_prime_inv());
      return result;
    }

    /* Little modular reduction (e.g. after a simple addition). */
    static inline 
    double RES_soft_reduce (double a){
      double p = get_prime();
        double b = 2*a;
        return (b>p) ? a-p :
            ((b<-p) ? a+p : a);
    }

    
    /* -a */
    static inline 
    double RES_negate(double a){
        return RES_soft_reduce(-a);
    }


    /* a*b */
    static inline 
    double RES_mul (double a, double b){
        double c = a*b;

        double r = RES_reduce(c);

//         printf("mul a = %f b = %f res = %f\n", a, b, r);

        return RES_reduce(c);
    }


    /* a+b */
    static inline 
    double RES_add (double a, double b){
        double c = a+b;
        return RES_soft_reduce(c);
    }

    
    static inline 
    double RES_pow(double a, unsigned e) {
        if ( 0==e )  return 1.0;
        double z = a;
        double y = 1;
        while ( 1 )
        {
            if ( e & 1 )  y = RES_mul(y, z);  // y *= z;
            e >>= 1;
            if ( 0==e )  break;
            z = RES_mul(z, z);  // z *= z;
        }
        return  y;
    }

    /* a^-1, using Bezout (extended Euclidian algorithm). */
    static inline 
    double RES_inv (double ri1){
        double bi = 0.0;
        double bi1 = 1.0;
        double ri = get_prime();
        double p, tmp, tmp2;

        int niters = 50;    

//         printf("CST_CUT = %g\n", CST_CUT);
        while (fabs(ri1) != 1.0)
        {
            p = RES_round(ri/ri1);
            tmp = bi - p * bi1;
            tmp2 = ri - p * ri1;
            bi = bi1;
            ri = ri1;
            bi1 = tmp;
            ri1 = tmp2;
//             printf("p = %g ri = %g ri1 = %g tmp = %g tmp2 = %g\n", p,
//                 ri, ri1, tmp, tmp2);
            if(niters-- == 0)break;
        };

        return ri1 * RES_soft_reduce(bi1);	/* Quicker !!!! */
    }
    
    /* a/b */
    static inline 
    double RES_div (double a, double b){
        return RES_mul(a, RES_inv(b));
    }    

public:
    /*! \brief sets the current prime. 
     *  
     *  Note that you are going to change a static member!
     *  \pre p is prime, but we abstained from such a test.
     *  \pre 0 < p < 2^26
     *  
     */
    static int 
    set_current_prime(int p){   
      int old_prime = get_prime_int();  

      modulus = p;
      prime = double(p);
      prime_inv = 1.0 / prime;
      return old_prime; 
    }
 
  /*! \brief return the current prime.  */
    static int get_current_prime(){
      return get_prime_int();
    }
  
  int  get_value() const{
    return int(x_);
  }
    
private:
    double x_;

public: 

    //! constructor of Residue, from int 
    Residue(int n = 0){
        x_= RES_reduce(n);
    }

    //! constructor of Residue, from long 
    Residue(long n){
        x_= RES_reduce(n);
    }

    Residue(unsigned int n) {
        x_ = RES_reduce(n);
    }
   
    Residue(unsigned long n) {
        x_ = RES_reduce(n);
    }

    //! Access operator for x, \c const 
    const double& x() const { return x_; }
    //! Access operator for x
    double&       x()       { return x_; }                     

    Self& operator += (const Self& p2) { 
        x() = RES_add(x(),p2.x()); 
        return (*this); 
    }
    Self& operator -= (const Self& p2){ 
        x() = RES_add(x(),RES_negate(p2.x())); 
        return (*this); 
    }
    Self& operator *= (const Self& p2){ 
        x() = RES_mul(x(),p2.x()); 
        return (*this); 
    }
    Self& operator /= (const Self& p2) { 
        x() = RES_div(x(),p2.x()); 
        return (*this); 
    }

        // 
    Self& operator += (int p2) { 
        x() = RES_add(x(),Residue(p2).x()); 
        return (*this); 
    }
    Self& operator -= (int p2){ 
        x() = RES_add(x(),Residue(-p2).x()); 
        return (*this); 
    }

    Self& operator *= (int p2){ 
        x() = RES_mul(x(),Residue(p2).x()); 
        return (*this); 
    }

    Self& operator /= (int p2) { 
        x() = RES_div(x(),Residue(p2).x()); 
        return (*this); 
    }

    Self pow(unsigned e) {
        Self y;
        y.x() = RES_pow(x(), e);
        return y;
    }

    friend Self mod_inverse(const Self&);
    friend Self operator + (const Self&);
    friend Self operator - (const Self&);

/*    friend Self operator + (const Self&, const Self&);
    friend Self operator - (const Self&, const Self&);
    friend Self operator * (const Self&, const Self&);
    friend Self operator / (const Self&, const Self&);  */                 
};

inline Residue mod_inverse(const Residue& p1) { 
    Residue r;
    r.x() = Residue::RES_inv(p1.x());
    return r;
}

inline Residue operator + (const Residue& p1)
{ return p1; }

inline Residue operator - (const Residue& p1){ 
    typedef Residue RES;
    Residue r; 
    r.x() = RES::RES_negate(p1.x());
    return r; 
}

inline Residue operator + (const Residue& p1, const Residue& p2) {

    Residue p(p1);
    p += p2;
    return p;
}

inline Residue operator - (const Residue& p1, const Residue& p2) {
    Residue p(p1);
    p -= p2;
    return p;
}

inline Residue operator * (const Residue& p1, const Residue& p2) {
   Residue p(p1);
   p *= p2;
   return p;
}
    
inline Residue operator / (const Residue& p1, const Residue& p2) {
   Residue p(p1);
   p /= p2;
   return p;
}            

inline bool operator == (const Residue& p1, const Residue& p2)
{ return ( p1.x()==p2.x() ); }   
inline bool operator == (const Residue& p1, int p2)
{ return ( p1 == Residue(p2) ); }   

inline bool operator != (const Residue& p1, const Residue& p2)
{ return (!( p1.x()==p2.x() )); }   

inline bool operator < (const Residue& p1, const Residue& p2)
{ return ( p1.x() < p2.x() ); }   
inline bool operator < (const Residue& p1, int p2)
{ return ( p1.x() < Residue(p2).x() ); }   

inline bool operator > (const Residue& p1, const Residue& p2)
{ return ( p1.x() > p2.x() ); }   

// I/O 
inline std::ostream& operator << (std::ostream& os, const Residue& p) {   
    typedef Residue RES;
    os << int(p.x()) ;
    return os;
}

int Residue::modulus = 67111067;
double Residue::prime = 67111067.0;
double Residue::prime_inv =1/67111067.0;
const double Residue::CST_CUT = ldexp(3., 51);

#endif // CGAL_RESIDUE_TYPE_H
