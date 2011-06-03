// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_GENERIC_STRUCTURE_H
#define CGAL_LINBOX_GENERIC_STRUCTURE_H

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <linbox/integer.h>

namespace CGAL{

        // This class provides the members of a generic algebraic structure
        // needed by LinBox. Fields and rings can inherit from this class.
        template<class _NT>
        class Linbox_generic_structure:
        public _NT{
                protected:
                typedef _NT                             NT;
                typedef LinBox::integer                 LI;
                typedef Linbox_generic_structure<NT>    LG;
                public:
                typedef NT                              Element;
                //typedef SOMETHING                       RandIter;

                public:

                Linbox_generic_structure(int p=0,int exp=1){
                        CGAL_assertion_msg(p==0,"modulus must be 0");
                        CGAL_assertion_msg(exp==1,"exponent must be 1");
                };
                Linbox_generic_structure(const LG&){};
                Linbox_generic_structure& operator=(const LG&){return *this;};
                ~Linbox_generic_structure(){};

                Element& assign(Element &x,const Element &y)const{
                        return (x=y);
                }

                LI& cardinality(LI &c)const{
                        // TODO: return a non-negative value for finite
                        // domains
                        return (c=-1);
                }

                // LinBox < 1.2 requires an LI characteristic
                LI& characteristic(LI &c)const{
                        // TODO: return correct characteristic for finite
                        // domains
                        return (c=0);
                }

                // LinBox 1.2 requires a long unsigned characteristic
                long unsigned& characteristic(long unsigned &c)const{
                        // TODO: return correct characteristic for finite
                        // domains
                        return (c=0);
                }

                bool areEqual(const Element &x,const Element &y)const{
                        return (x==y);
                }

#define CGAL_LINBOX_BINARY_OPERATION(_name,_op) \
                Element& _name(Element &x, \
                               const Element &y, \
                               const Element &z)const{ \
                        return (x=y _op z); \
                }
                CGAL_LINBOX_BINARY_OPERATION(add,+)
                CGAL_LINBOX_BINARY_OPERATION(sub,-)
                CGAL_LINBOX_BINARY_OPERATION(mul,*)
#undef CGAL_LINBOX_BINARY_OPERATION

                Element& div(Element &x,
                             const Element &y,
                             const Element &z)const{
                        x=y/z;
                        CGAL_warning_msg(x*z==y,"division was not exact");
                        return x;
                }

                Element& neg(Element &x,const Element &y)const{
                        x=-y;
                        return x;
                }

                Element& inv(Element &x,const Element &y)const{
                        x=1/y;
                        CGAL_warning_msg(x*y==1,"inverse was not exact");
                        return x;
                }

                Element& axpy(Element &r,
                              const Element &a,
                              const Element &x,
                              const Element &y)const{return (r=a*x+y);}

                bool isZero(const Element &x)const{return (x==Element(0));}

                bool isOne(const Element &x)const{return (x==Element(1));}

#define CGAL_LINBOX_INPLACE_OPERATION(_name,_op) \
                Element& _name(Element &x,const Element &y)const{ \
                        return (x _op y); \
                }
                CGAL_LINBOX_INPLACE_OPERATION(addin,+=)
                CGAL_LINBOX_INPLACE_OPERATION(subin,-=)
                CGAL_LINBOX_INPLACE_OPERATION(mulin,*=)
#undef CGAL_LINBOX_INPLACE_OPERATION
                Element& divin(Element &x,const Element &y)const{
                        CGAL_warning_msg(x=(x/y)*y,"division was not exact");
                        return (x/=y);
                }

                Element& negin(Element &x)const{
                        x=-x;
                        return x;
                }

                Element& invin(Element &x)const{
                        CGAL_warning_msg(x=1/(1/x),"inverse was not exact");
                        return (x=1/x);
                }

                Element& axpyin(Element &r,
                                const Element &a,
                                const Element &x)const{return (r+=a*x);}

                long compare(const Element &a,const Element &b)const{
                        return (a==b?0:(a>b?1:-1));
                }

                int sign(const Element &x)const{
                        switch(x.sign()){
                                case CGAL::POSITIVE: return 1; break;
                                case CGAL::NEGATIVE: return -1; break;
                                default: return 0;
                        }
                }

                Element& abs(Element &x,const Element &a)const{
                        if(a.sign()==CGAL::NEGATIVE)
                                x=-a;
                        else
                                x=a;
                        CGAL_assertion(x>=0);
                        return x;
                }

                static int getMaxModulus(){
                        // TODO: return a positive modulus when the CGAL
                        // number type is finite
                        return -1;
                }
        };

} // namespace CGAL

#endif // CGAL_LINBOX_GENERIC_STRUCTURE_H
