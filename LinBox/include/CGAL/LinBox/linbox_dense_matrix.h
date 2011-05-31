// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_LINBOX_DENSE_MATRIX_H
#define CGAL_LINBOX_LINBOX_DENSE_MATRIX_H

#include <CGAL/basic.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
#include <CGAL/LinBox/linbox_dense_vector.h>
#include <linbox/blackbox/dense.h>
#include <linbox/matrix/matrix-domain.h>

namespace CGAL{

        // This class represents a matrix that meets both CGAL and LinBox
        // requirements. This is obtained by inheriting from a LinBox
        // matrix and adding the remaining functions from the CGAL
        // specification. The template parameter _LT is a LinBox type, that
        // is, a type that meets the LinBox archetype.
        template<class _LT>
        class Linbox_dense_matrix:
        public LinBox::DenseMatrix<_LT>{
                private:
                typedef _LT                             LT;
                typedef LinBox::DenseMatrix<LT>         LDM;
                typedef Linbox_dense_matrix<LT>         DM;

                // required CGAL types are made public
                public:
                typedef typename LT::Element            NT;
                typedef typename LDM::RawIterator       iterator;
                typedef typename LDM::RawIterator       row_iterator;
                typedef LinBox::Subiterator<iterator>   column_iterator;
                typedef Linbox_dense_vector<NT>         Vector;
                struct Identity{};
                // useful public types not required by CGAL
                typedef typename LDM::ConstRawIterator  const_iterator;
                typedef typename LDM::ConstRowIterator  const_row_iterator;
                typedef typename LDM::ConstColIterator  const_column_iterator;
                typedef LinBox::MatrixCategories::RowColMatrixTag
                                                        MatrixCategory;
                typedef LinBox::MatrixDomain<LT>        LMD;
                private:
                typedef typename LMD::Permutation       Permutation;

                // constructors
                public:
                Linbox_dense_matrix():
                        _field(LT()),
                        _matrixdomain(_field),
                        LDM(_field){}
                Linbox_dense_matrix(int d):
                        _field(LT()),
                        _matrixdomain(_field),
                        LDM(_field,(size_t)d,(size_t)d){}
                Linbox_dense_matrix(const LT &f,int m,int n):
                        _field(f),
                        _matrixdomain(f),
                        LDM(f,(size_t)m,(size_t)n){}
                Linbox_dense_matrix(int m,int n):
                        _field(LT()),
                        _matrixdomain(_field),
                        LDM(_field,(size_t)m,(size_t)n){}
                Linbox_dense_matrix(std::pair<int,int> p):
                        _field(LT()),
                        _matrixdomain(_field),
                        LDM(_field,(size_t)p.first,(size_t)p.second){}

                Linbox_dense_matrix(int m,int n,const NT &x):
                        _field(LT()),
                        _matrixdomain(_field),
                        LDM(_field,(size_t)m,(size_t)n){
                        for(iterator i=this->begin();i!=this->end();++i)
                                *i=x;
                }

                Linbox_dense_matrix(int n,
                                    const Identity &id,
                                    const NT &x=NT(1)):
                        _field(LT()),
                        _matrixdomain(_field),
                        LDM(_field,(size_t)n,(size_t)n){
                        for(int i=0;i<n;++i)
                                this->setEntry(i,i,x);
                }

                // constructor from an iterator range, containing n column
                // vectors of size m
                template <class Forward_iterator>
                Linbox_dense_matrix(Forward_iterator first,
                                    Forward_iterator last):
                                _field(LT()),
                                _matrixdomain(_field),
                                LDM(_field,
                                    std::distance(first->begin(),
                                                  first->end()),
                                    std::distance(first,last)){
                        typedef std::iterator_traits<Forward_iterator>  IT;
                        typedef typename IT::value_type                 V;
                        typedef typename V::const_iterator              VI;

                        size_t column=0,row;
                        for(Forward_iterator ci=first;ci!=last;++ci){
                                row=0;
                                for(VI ri=ci->begin();ri!=ci->end();++ri){
                                        this->setEntry(row,column,*ri);
                                        ++row;
                                }
                                ++column;
                        }
                }

                // constructor from a vector of n columns of size m
                Linbox_dense_matrix(std::vector<Vector> &A):
                                _field(LT()),
                                _matrixdomain(_field),
                                LDM(_field,
                                    std::distance((A.begin())->begin(),
                                                  (A.begin())->end()),
                                    A.size()){
                        typedef typename std::vector<Vector>            VV;
                        typedef typename VV::const_iterator             VVI;
                        typedef typename Vector::const_iterator         VI;
                        size_t column=0,row;
                        for(VVI ci=A.begin();ci!=A.end();++ci){
                                row=0;
                                for(VI ri=ci->begin();ri!=ci->end();++ri){
                                        this->setEntry(row,column,*ri);
                                        ++row;
                                }
                                ++column;
                        }
                }

                // member functions
                public:

                // not needed by the CGAL concept
                const LT& get_field()const{return _field;}
                const LMD& get_matrixdomain()const{return _matrixdomain;}

                // following functions are needed by the CGAL concept
                // n, number of rows
                int row_dimension()const{return this->rowdim();}
                // m, number of columns
                int column_dimension()const{return this->coldim();}
                // returns (m,n)
                std::pair<int,int> dimension()const{
                        return std::make_pair(this->coldim(),this->rowdim());
                }

                const NT& operator()(int i,int j)const{
                        return this->getEntry(i,j);
                }

                // non-const iterators
                row_iterator row_begin(int i)
                        {return (this->rawBegin())+(i*this->coldim());}
                row_iterator row_end(int i)
                        {return (this->row_begin(i))+(this->coldim());}
                column_iterator column_begin(int i)
                        {return (this->colBegin())[i].begin();}
                column_iterator column_end(int i)
                        {return (this->colBegin())[i].end();}
                iterator begin(){return this->rawBegin();}
                iterator end(){return this->rawEnd();}

                // const iterators
                const_row_iterator row_begin(int i)const
                        {return (this->rawBegin())+(i*this->coldim());}
                const_row_iterator row_end(int i)const
                        {return (this->row_begin(i))+(this->coldim());}
                const_column_iterator column_begin(int i)const
                        {return (this->colBegin())[i].begin();}
                const_column_iterator column_end(int i)const
                        {return (this->colBegin())[i].end();}
                const_iterator begin()const{return this->rawBegin();}
                const_iterator end()const{return this->rawEnd();}

                bool operator==(const DM &M1)const{
                        if(this->dimension()!=M1.dimension())
                                return false;
                        const_iterator j=M1.begin();
                        for(const_iterator i=this->begin();i!=this->end();++i){
                                if(*i!=*j)
                                        return false;
                                ++j;
                        }
                        return true;
                }

                bool operator!=(const DM &M1)const{
                        return !(*this==M1);
                }

                void swap_rows(int i,int j){
                        Permutation perm;
                        perm.push_back(std::make_pair(i,j));
                        _matrixdomain.permuteRows(*this,
                                                  perm.begin(),
                                                  perm.end());
                }

                void swap_columns(int i,int j){
                        Permutation perm;
                        perm.push_back(std::make_pair(i,j));
                        _matrixdomain.permuteColumns(*this,
                                                     perm.begin(),
                                                     perm.end());
                }

                Vector row(int i)const{
                        Vector r;
                        std::back_insert_iterator<Vector> ins(r);
                        const_iterator first=begin()+i*column_dimension();
                        std::copy(first,first+column_dimension(),ins);
                        return r;
                }

                Vector column(int i)const{
                        Vector r;
                        const_iterator index=this->rawBegin()+i;
                        const_iterator
                        last=index+row_dimension()*column_dimension();
                        for(;index!=last;index+=column_dimension())
                                r.push_back(*index);
                        return r;
                }

                DM& operator*(const DM &M1)const{
                        CGAL_precondition(row_dimension()==
                                          M1.column_dimension());
                        DM *C=new DM(column_dimension(),M1.row_dimension());
                        return _matrixdomain.mul(*C,*this,M1);
                }

                DM& operator+(const DM &M1)const{
                        CGAL_precondition(row_dimension()==M1.row_dimension());
                        CGAL_precondition(column_dimension()==
                                          M1.column_dimension());
                        DM *C=new DM(column_dimension(),row_dimension());
                        return _matrixdomain.add(*C,*this,M1);
                }

                DM& operator-(const DM &M1)const{
                        CGAL_precondition(row_dimension()==M1.row_dimension());
                        CGAL_precondition(column_dimension()==
                                          M1.column_dimension());
                        DM *C=new DM(column_dimension(),row_dimension());
                        return _matrixdomain.sub(*C,*this,M1);
                }

                DM& operator-()const{
                        DM *B=new DM(column_dimension(),row_dimension());
                        return _matrixdomain.neg(*B,*this);
                }

                DM& operator*(const NT &x)const{
                        DM *B=new DM(column_dimension(),row_dimension());
                        return _matrixdomain.mul(*B,*this,x);
                }

                Vector& operator*(const Vector &vec)const{
                        CGAL_assertion(column_dimension()==vec.dimension());
                        Vector *result=new Vector(vec.dimension());
                        return _matrixdomain.vectorMul(*result,*this,vec);
                }

                // data members
                private:
                LT _field;
                LMD _matrixdomain;

        }; // class Linbox_dense_matrix

        template <class T,class U>
        Linbox_dense_matrix<U>& operator*(const T &x,
                                          const Linbox_dense_matrix<U> &m){
                return m*x;
        };

        template <class T>
        std::ostream& operator<<(std::ostream &o,
                                 const Linbox_dense_matrix<T> &m){
                if(is_pretty(o))
                        o<<"Linbox_dense_matrix[";
                else
                        o<<'[';
                if(m.row_dimension()>0){
                        o<<m.row(0);
                        for(int i=1;i<m.row_dimension();++i)
                                o<<','<<m.row(i);
                }
                o<<']';
                return o;
        };

} // namespace CGAL

#endif // CGAL_LINBOX_LINBOX_DENSE_MATRIX_H
