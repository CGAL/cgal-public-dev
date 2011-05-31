// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_LINBOX_INTEGER_RANDITER_H
#define CGAL_LINBOX_LINBOX_INTEGER_RANDITER_H

#include <linbox/integer.h>
#include <linbox/element/abstract.h>
#include <CGAL/Gmpz.h>

namespace CGAL{

// The template parameter _IR is an integer ring, not an integer type.
template <class _IR>
class IntegerRingRandIter{
        public:
        typedef _IR                                     IR;
        typedef IntegerRingRandIter<IR>                 IRRI;
        typedef typename IR::Element                    Element;
        typedef LinBox::integer                         integer;
        typedef LinBox::ElementAbstract                 ElementAbstract;

        IntegerRingRandIter(const IR&,const integer&,const integer&);
        IntegerRingRandIter(const IRRI&);
        ~IntegerRingRandIter();
        IntegerRingRandIter& operator=(const IRRI&);
        Element &random(Element&)const;
        ElementAbstract& random(ElementAbstract&)const;
        private:
        IR _F;
        integer _size;
        integer _seed;
}; // class IntegerRingRandIter

template <class IR_>
IntegerRingRandIter<IR_>::
IntegerRingRandIter(const IR &F,const integer &size=0,const integer &seed=0):
_F(F),_size(size),_seed(seed){
        if(seed==0)
                _seed=time(NULL);
};

template <class IR_>
IntegerRingRandIter<IR_>::
IntegerRingRandIter(const IntegerRingRandIter<IR_> &R):
_F(R._F),_size(R._size),_seed(R._seed){};

template <class IR_>
IntegerRingRandIter<IR_>::
~IntegerRingRandIter(){};

template <class IR_>
IntegerRingRandIter<IR_>&
IntegerRingRandIter<IR_>::operator=(const IntegerRingRandIter<IR_> &R){
        if(this!=&R){
                _F=R._F;
                _seed=R._seed;
                _size=R._size;
        }
        return *this;
};

template <class IR_>
typename IR_::Element&
IntegerRingRandIter<IR_>::random(typename IR_::Element &a)const{
        CGAL_error_msg("not implemented");
        return a;
};

template <class IR_>
LinBox::ElementAbstract&
IntegerRingRandIter<IR_>::random(LinBox::ElementAbstract &a)const{
        CGAL_error_msg("not implemented");
        return a;
};

// Specialized function members for Gmpz integer ring.

/*
template <>
IntegerRingRandIter<Linbox_integer_ring<Gmpz> >::
IntegerRingRandIter(const Linbox_integer_ring<Gmpz> &F,
                    const integer &size,
                    const integer &seed):
_F(F),_size(size){
        gmp_randinit_default((gmp_randstate_t)seed);
};

template <>
IntegerRingRandIter<Linbox_integer_ring<Gmpz> >::
IntegerRingRandIter(const IntegerRingRandIter<Linbox_integer_ring<Gmpz> > &R):
_F(R._F),_size(R._size){
        gmp_randinit_set((gmp_randstate_t)_seed,(gmp_randstate_t)R._seed);
};

template <>
IntegerRingRandIter<Linbox_integer_ring<Gmpz> >&
IntegerRingRandIter<Linbox_integer_ring<Gmpz> >::operator=
                (const IntegerRingRandIter<Linbox_integer_ring<Gmpz> > &R){
        if(this!=&R){
                _F=R._F;
                _size=R._size;
                gmp_randinit_set((gmp_randstate_t)_seed,
                                 (gmp_randstate_t)R._seed);
        }
        return *this;
};

template <>
Gmpz&
IntegerRingRandIter<Linbox_integer_ring<Gmpz> >::random(Gmpz &a)const{
        CGAL_assertion(size>0);
        mpz_urandomm(a.mpz(),(gmp_randstate_t)seed,SpyInteger::get_mpz(size));
        return a;
};

template <>
LinBox::ElementAbstract&
IntegerRingRandIter<Linbox_integer_ring<Gmpz> >::random
                (LinBox::ElementAbstract &a)const{
        CGAL_error_msg("not implemented");
        return a;
};

*/
} // namespace CGAL

#endif // CGAL_LINBOX_LINBOX_INTEGER_RANDITER_H
