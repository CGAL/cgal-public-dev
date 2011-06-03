// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_DENSE_VECTOR_H
#define CGAL_LINBOX_DENSE_VECTOR_H

#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
#include <CGAL/memory.h>
#include <linbox/vector/vector-traits.h>

namespace CGAL{

        template <class _NT,class _AL=CGAL_ALLOCATOR(_NT)>
        class Linbox_dense_vector:
        public std::vector<_NT,_AL>{
                public:
                typedef _NT                             NT;
                typedef _AL                             AL;
                typedef std::vector<NT,AL>              StdVec;
                typedef typename StdVec::iterator       iterator;
                typedef typename StdVec::const_iterator const_iterator;
                typedef LinBox::VectorCategories::DenseVectorTag
                                                        VectorCategory;

                template<typename _Tp1>
                struct rebind{
                        typedef typename AL::template rebind<_Tp1> other;
                };

                private:
                typedef Linbox_dense_vector<NT,AL>      DV;
                public:

                // constructors

                Linbox_dense_vector():StdVec(){}
                Linbox_dense_vector(int d):StdVec(d){}
                Linbox_dense_vector(int d,const NT x):StdVec(d,x){}
                Linbox_dense_vector(StdVec &v1):StdVec(v1){}

                // member functions

                int dimension()const{return StdVec::size();}

                bool is_zero()const{
                        const_iterator i;
                        NT zero(0);
                        for(i=StdVec::begin();i!=StdVec::end();++i){
                                if(*i!=zero)
                                        return false;
                        }
                        return true;
                }

                const NT& operator[](int i)const{
                        CGAL_precondition_msg(0<=i&&i<=dimension(),
                                              "vector index out of bounds");
                        return StdVec::operator[](i);
                }

#define CGAL_LINBOX_VECTOR_OPERATOR(_op,_fun) \
                DV _op(DV &v1)const{ \
                        CGAL_precondition_msg(dimension()==v1.dimension(), \
                                              "vector dimensions differ"); \
                        StdVec v2(dimension()); \
                        std::transform(StdVec::begin(), \
                                       StdVec::end(), \
                                       v1.begin(), \
                                       v2.begin(), \
                                       _fun); \
                        return DV(v2); \
                }
                CGAL_LINBOX_VECTOR_OPERATOR(operator+,std::plus<NT>())
                CGAL_LINBOX_VECTOR_OPERATOR(operator-,std::minus<NT>())
                CGAL_LINBOX_VECTOR_OPERATOR(operator*,std::multiplies<NT>())
#undef CGAL_LINBOX_VECTOR_OPERATOR

                DV operator-()const{
                        StdVec v2(dimension());
                        std::transform(StdVec::begin(),
                                       StdVec::end(),
                                       v2.begin(),
                                       std::negate<NT>());
                        return DV(v2);
                }

#define CGAL_LINBOX_VECTOR_SELF_OPERATOR(_op,_fun) \
                DV& _op(DV &v1){ \
                        CGAL_precondition_msg(dimension()==v1.dimension(), \
                                              "vector dimensions differ"); \
                        std::transform(StdVec::begin(), \
                                       StdVec::end(), \
                                       v1.begin(), \
                                       StdVec::begin(), \
                                       _fun); \
                        return *this; \
                }
                CGAL_LINBOX_VECTOR_SELF_OPERATOR(operator+=,std::plus<NT>());
                CGAL_LINBOX_VECTOR_SELF_OPERATOR(operator-=,std::minus<NT>());
                CGAL_LINBOX_VECTOR_SELF_OPERATOR(operator*=,
                                                 std::multiplies<NT>());
                CGAL_LINBOX_VECTOR_SELF_OPERATOR(operator/=,
                                                 std::divides<NT>());
#undef CGAL_LINBOX_VECTOR_SELF_OPERATOR

                DV operator*(const NT &R)const{
                        StdVec v2(dimension());
                        std::transform(StdVec::begin(),
                                       StdVec::end(),
                                       v2.begin(),
                                       std::bind1st(std::multiplies<NT>(),R));
                        return DV(v2);
                }
        };

        template <class NT,class AL>
        Linbox_dense_vector<NT,AL>
        operator*(const NT &R,const Linbox_dense_vector<NT,AL> &v){
                return v*R;
        };

        template <class _NT,class _AL>
        std::ostream&
        operator<<(std::ostream &o,const Linbox_dense_vector<_NT,_AL> &v){
                typedef _NT                                     NT;
                typedef _AL                                     AL;
                typedef Linbox_dense_vector<NT,AL>              V;
                typedef typename V::const_iterator              CI;
                if(is_pretty(o))
                        o<<"Linbox_dense_vector[";
                else
                        o<<'[';
                CI i=v.begin();
                if(i!=v.end()){
                        o<<(*i);
                        for(++i;i!=v.end();++i)
                                o<<','<<(*i);
                }
                o<<']';
                return o;
        };

} // namespace CGAL

#endif // CGAL_LINBOX_DENSE_VECTOR_H
