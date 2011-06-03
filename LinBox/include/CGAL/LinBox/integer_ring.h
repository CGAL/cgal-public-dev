// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_INTEGER_RING_H
#define CGAL_LINBOX_INTEGER_RING_H

#include <CGAL/LinBox/generic_structure.h>
#include <CGAL/Gmpz.h>
#include <linbox/field/field-traits.h>
#include <CGAL/LinBox/integer_randiter.h>

namespace CGAL{

        // This class is a wrapper for CGAL classes: it takes as template
        // parameter a CGAL integer type, which must be constructible from
        // Gmpz. This way, any CGAL integer type can be used as a number
        // type in LinBox (for instance, by parameterizing the LinBox
        // matrices).
        template<class _IT>
        class Linbox_integer_ring:
        public Linbox_generic_structure<_IT>{
                private:
                typedef _IT                             IT;
                typedef Linbox_generic_structure<_IT>   LG;
                typedef Linbox_integer_ring<IT>         LR;
                typedef typename LG::LI                 LI;

                public:
                typedef IntegerRingRandIter<LR>         RandIter;
                typedef typename LG::Element            Element;

                template<class T>
                Element& init(Element &x,const T &y=T(0))const{
                        return (x=y);
                }

                Element& init(Element &x,const LI &y)const{
                        return x=Gmpz(LinBox::SpyInteger::get_mpz(y));
                }

                LI& convert(LI &x,const Element &y)const{
                        mpz_set(x.get_mpz(),y.mpz());
                        return x;
                }

                double& convert(double &x,const Element &y)const{
                        x=mpz_get_d(y.mpz());
                        CGAL_exactness_warning_msg(mpz_cmp_d(y.mpz(),x)==0,
                                                   "conversion was not exact");
                        return x;
                }

                float& convert(float &x,const Element &y)const{
                        x=(float)mpz_get_d(y.mpz());
                        CGAL_exactness_warning_msg(
                                        mpz_cmp_d(y.mpz(),(double)x)==0,
                                        "conversion was not exact");
                        return x;
                }

                Element& gcd(Element &g,
                             const Element &a,
                             const Element &b)const{
                        mpz_gcd(g.mpz(),a.mpz(),b.mpz());
                        return g;
                }

                Element& gcdin(Element &g,const Element &b)const{
                        mpz_gcd(g.mpz(),g.mpz(),b.mpz());
                        return g;
                }

                Element& xgcd(Element &g,
                              Element &s,
                              Element &t,
                              const Element &a,
                              const Element &b)const{
                        mpz_gcdext(g.mpz(),s.mpz(),t.mpz(),a.mpz(),b.mpz());
                        return g;
                }

                Element& lcm(Element &c,
                             const Element &a,
                             const Element &b)const{
                        mpz_gcd(c.mpz(),a.mpz(),b.mpz());
                        return c;
                }

                Element& lcmin(Element &c,const Element &b)const{
                        mpz_gcd(c.mpz(),c.mpz(),b.mpz());
                        return c;
                }

                Element& sqrt(Element &x,const Element &y)const{
                        mpz_sqrt(x.mpz(),y.mpz());
                        return x;
                }

                Element& quo(Element &q,
                             const Element &a,
                             const Element &b)const{
                        return (q=a/b);
                }

                Element& rem(Element &r,
                             const Element &a,
                             const Element &b)const{
                        mpz_tdiv_r(r.mpz(),a.mpz(),b.mpz());
                        return r;
                }

                Element& quoin(Element &a,const Element &b)const{
                        return (a=a/b);
                }

                Element& remin(Element &r,const Element &b)const{
                        mpz_tdiv_r(r.mpz(),r.mpz(),b.mpz());
                        return r;
                }

                bool isDivisor(const Element &a,const Element &b)const{
                        return mpz_divisible_p(b.mpz(),a.mpz());
                }

                std::ostream& write(std::ostream &os)const{
                        os<<"Linbox integer ring constructed from a CGAL type";
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

        template <class RT_>
        struct ClassifyRing<CGAL::Linbox_integer_ring<RT_> >{
                typedef RingCategories::IntegerTag              categoryTag;
        };

} // namespace LinBox

#endif // CGAL_LINBOX_INTEGER_RING_H
