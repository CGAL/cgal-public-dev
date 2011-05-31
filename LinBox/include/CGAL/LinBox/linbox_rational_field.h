// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_LINBOX_RATIONAL_FIELD_H
#define CGAL_LINBOX_LINBOX_RATIONAL_FIELD_H

#include <CGAL/LinBox/linbox_generic_structure.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Fraction_traits.h>
#include <linbox/integer.h>
#include <linbox/field/field-traits.h>

namespace CGAL{

        // This class is a wrapper for CGAL classes: it takes as template
        // parameter a CGAL rational field type, which must be
        // constructible from Gmpz. This way, any CGAL rational field can
        // bu used as a number type in LinBox (for instance, by
        // parameterizing the LinBox matrices).
        template<class _FT>
        class Linbox_rational_field:
        public Linbox_generic_structure<_FT>{
                private:
                typedef _FT                             FT;
                typedef Linbox_generic_structure<FT>    LG;
                typedef Linbox_rational_field<FT>       LF;
                typedef typename LG::LI                 LI;
                /*typedef Fraction_traits<FT>             Ftraits;
                typedef typename Ftraits::Compose       Compose;
                typedef typename Ftraits::Decompose     Decompose;*/

                public:
                typedef typename LG::Element            Element;

                Linbox_rational_field(int p=0,int exp=1):LG(p,exp){};

                Linbox_rational_field(const LF &f):LG(){};

                LF& operator=(const LF&){return *this;};

                ~Linbox_rational_field(){};

                template<class T>
                Element& init(Element &x,const T &y=T(0))const{
                        return (x=y);
                }

                Element& init(Element &x,const LI &num,const LI &den)const{
                        /*Gmpz z1(0),z2(0);
                        mpz_swap(z1.mpz(),num.get_mpz());
                        mpz_swap(z2.mpz(),den.get_mpz());
                        x=Compose()(z1,z2);
                        mpz_swap(z1.mpz(),num.get_mpz());
                        mpz_swap(z2.mpz(),den.get_mpz());
                        return x;*/
                        x=Gmpz(LinBox::SpyInteger::get_mpz(num));
                        x/=Gmpz(LinBox::SpyInteger::get_mpz(den));
                        return x;
                }

                // this function returns floor(num/den)
                LI& convert(LI &x,const Element &y=0)const{
                        /*Gmpz num,den;
                        Decompose()(y,num,den);
                        num/=den;
                        mpz_set(x.get_mpz(),num.mpz());
                        CGAL_assertion(mpz_cmp(x.get_mpz(),num.mpz())==0);
                        return x;
                        */
                        mpz_tdiv_q(x.get_mpz(),
                                   mpq_numref(y.mpq()),
                                   mpq_denref(y.mpq()));
                        return x;
                }

                double& convert(double &x,const Element &y=0)const{
                        x=mpq_get_d(y.mpq());
                        CGAL_exactness_warning_code(mpq_t z;)
                        CGAL_exactness_warning_code(mpq_init(z);)
                        CGAL_exactness_warning_code(mpq_set_d(z,x);)
                        CGAL_exactness_warning_msg(mpq_cmp(z,y.mpz())==0,
                                                   "conversion was not exact");
                        CGAL_exactness_warning_code(mpq_clear(z);)
                        return x;
                }

                float& convert(float &x,const Element &y=0)const{
                        x=(float)mpq_get_d(y.mpq());
                        CGAL_exactness_warning_code(mpq_t z;)
                        CGAL_exactness_warning_code(mpq_init(z);)
                        CGAL_exactness_warning_code(mpq_set_f(z,x);)
                        CGAL_exactness_warning_msg(mpq_cmp(z,y.mpz())==0,
                                                   "conversion was not exact");
                        CGAL_exactness_warning_code(mpq_clear(z);)
                        return x;
                }

                LI& get_num(LI &x,const Element &y)const{
                        /*Gmpz num,den;
                        Decompose()(y,num,den);
                        mpz_set(x.get_mpz(),num.mpz());
                        return x;
                        */
                        mpz_set(x.get_mpz(),mpq_numref(y.mpq()));
                        return x;
                }

                LI& get_den(LI &x,const Element &y)const{
                        /*Gmpz num,den;
                        Decompose()(y,num,den);
                        mpz_set(x.get_mpz(),den.mpz());
                        return x;*/
                        mpz_set(x.get_mpz(),mpq_denref(y.mpq()));
                        return x;
                }

                std::ostream& write(std::ostream &os)const{
                        os<<"Linbox rational field constructed from CGAL type";
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

namespace LinBox{

        template <class FT_>
        struct ClassifyRing<CGAL::Linbox_rational_field<FT_> >{
                typedef RingCategories::RationalTag             categoryTag;
        };

} // namespace LinBox

#endif // CGAL_LINBOX_LINBOX_RATIONAL_FIELD_H
