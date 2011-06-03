// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_MPQ_CLASS_FIELD_H
#define CGAL_LINBOX_MPQ_CLASS_FIELD_H

#include <CGAL/config.h>

#if !defined CGAL_USE_GMPXX && !defined CGAL_USE_TOPCOM
#error This needs GMPXX or TOPCOM
#else

#ifdef CGAL_USE_GMPXX
#include <gmpxx.h>
#else
#include <Rational.h> // include the TOPCOM wrapper for mpq_class
#endif

#include <CGAL/LinBox/rational_field.h>

namespace CGAL{

        template<>
        class Linbox_rational_field<mpq_class>:
        public Linbox_generic_structure<mpq_class>{
                private:
                typedef mpq_class                       FT;
                typedef Linbox_generic_structure<FT>    LG;
                typedef Linbox_rational_field<FT>       LF;
                typedef LG::LI                 LI;

                public:
                typedef LG::Element            Element;

                Linbox_rational_field(int p=0,int exp=1):LG(p,exp){};

                Linbox_rational_field(const LF &f):LG(){};

                LF& operator=(const LF&){return *this;};

                ~Linbox_rational_field(){};

                template<class T>
                Element& init(Element &x,const T &y=T(0))const{
                        x=Element(y);
                        return x;
                }

                Element& init(Element &x,const LI &num,const LI &den)const{
                        mpq_t q;
                        mpq_init(q);
                        mpq_set_num(q,LinBox::SpyInteger::get_mpz(num));
                        mpq_set_den(q,LinBox::SpyInteger::get_mpz(den));
                        mpq_clear(q);
                        // TODO: use mpq_numref and mpq_denref to avoid
                        // copies; then remove mpq_clear
                        return x=Element(q);
                }

                // this function returns floor(num/den)
                LI& convert(LI &x,const Element &y=0)const{
                        mpz_tdiv_q(x.get_mpz(),
                                   y.get_num_mpz_t(),
                                   y.get_den_mpz_t());
                        return x;
                }

                double& convert(double &x,const Element &y=0)const{
                        x=mpq_get_d(y.get_mpq_t());
                        CGAL_exactness_warning_code(mpq_t z;)
                        CGAL_exactness_warning_code(mpq_init(z);)
                        CGAL_exactness_warning_code(mpq_set_d(z,x);)
                        CGAL_exactness_warning_msg(mpq_cmp(z,y.mpz())==0,
                                                   "conversion was not exact");
                        CGAL_exactness_warning_code(mpq_clear(z);)
                        return x;
                }

                float& convert(float &x,const Element &y=0)const{
                        x=(float)mpq_get_d(y.get_mpq_t());
                        CGAL_exactness_warning_code(mpq_t z;)
                        CGAL_exactness_warning_code(mpq_init(z);)
                        CGAL_exactness_warning_code(mpq_set_f(z,x);)
                        CGAL_exactness_warning_msg(mpq_cmp(z,y.mpz())==0,
                                                   "conversion was not exact");
                        CGAL_exactness_warning_code(mpq_clear(z);)
                        return x;
                }

                LI& get_num(LI &x,const Element &y)const{
                        mpz_set(x.get_mpz(),y.get_num_mpz_t());
                        return x;
                }

                LI& get_den(LI &x,const Element &y)const{
                        mpz_set(x.get_mpz(),y.get_den_mpz_t());
                        return x;
                }

                std::ostream& write(std::ostream &os)const{
                        os<<"Linbox field constructed from mpq_class";
                        return os;
                }

                std::istream& read(std::istream &is){return is;}

                std::ostream& write(std::ostream &os,const Element &x)const{
                        return (os<<x);
                }

                std::istream& read(std::istream &is,Element &x)const{
                        return (is>>x);
                }

        };

} // namespace CGAL

#endif // !defined CGAL_USE_GMPXX && !defined CGAL_USE_TOPCOM
#endif // CGAL_LINBOX_MPQ_CLASS_FIELD_H
