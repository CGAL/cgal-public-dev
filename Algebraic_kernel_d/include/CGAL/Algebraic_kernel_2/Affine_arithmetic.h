// Copyright (c) 2009, 2010, 2011, 2012 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://eric@scm.gforge.inria.fr/svn/cgal/branches/features/Symbolic-mpi/Algebraic_kernel_d/include/CGAL/Algebraic_kernel_d/Curve_analysis_2.h $
// $Id: Curve_analysis_2.h 71840 2012-08-30 11:29:52Z eric $
// 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_2_AFFINE_ARITHMETIC_H
#define CGAL_ALGEBRAIC_KERNEL_2_AFFINE_ARITHMETIC_H

/*!\file Affine_arithmetic.h
 * \brief this file defines the class \c Affine_form
 */

#include <CGAL/config.h>

namespace CGAL {

template <class NT> class Affine_form;

template <class NT>
std::ostream & operator <<(std::ostream&, const Affine_form<NT>&);

namespace internal {

// \brief tag type to distinguish a certain constructor of \c Affine_form
class Creation_tag {};

template <class NT_>
class Affine_form_rep
{
    typedef NT_ NT; 
    
    typedef Affine_form_rep<NT> Self;   
    
    typedef std::vector<NT>       Vector;
    typedef std::vector<unsigned> Int_vector;
    
    Affine_form_rep(Creation_tag, const NT& center_, unsigned n_coeffs,
        unsigned n_indexes) : 
        center(center_), coeff(n_coeffs, NT(0)), index(n_indexes, 0)
    { }

    Affine_form_rep(const NT& v0 = 0) : center(v0), coeff(), index()
    { }

    // Create an Affine_form from an array of NTs   
    Affine_form_rep(const NT& v0, const NT * t1, const unsigned * t2,
                unsigned t) 
    {
        /*length=t;
        center=v0;
        coeffs = new NT [length];
        indexes = new unsigned [length];
        
        for(unsigned i = 0; i < length; i++)
        {
            coeffs[i]=t1[i];
            indexes[i]=t2[i];
        }
        if(indexes[length-1] > last) set_default(indexes[length-1]);*/
    }
        
    Affine_form_rep(const Self& p, NT alpha, NT dzeta, 
        NT delta) : center(alpha*(p.center)+dzeta), coeff(), index()
    {
        /*coeffs = new NT [length];
        indexes = new unsigned [length];

        for (unsigned i=0; i < p.length;i++)
        {
            indexes[i]=p.indexes[i];
            coeffs[i]=alpha*(p.coeffs[i]);
        }

        indexes[p.length] = inc_last();  
        coeffs[p.length] = delta;*/
    }
    
    Affine_form_rep(const NT& left, const NT& right)
    {
        unsigned idx = inc_last();
        center = (right+left)/2;
        NT rad = (right-left)/2;
        coeff.push_back(rad);
        index.push_back(idx);
    }
    
    //! sets the highest noise symbol in use
    static void set_default(unsigned val = 0) { last = val; }
    
    //! increment the highest noise symbol in use: create a new noise symbol
    static unsigned inc_last()  
    { 
        last++; 
        return last;
    }
    
    friend class ::CGAL::Affine_form<NT>;

    static unsigned last;  // the highest noise symbol in use

    NT center;             // central value of an affine form
    
    Vector coeff;          // values of noise symbols
    Int_vector index;      // indexes of noise symbols

}; // class Affine_form_rep<NT_>

template <class NT> unsigned Affine_form_rep<NT>::last = 0; 

} // namespace internal

template <class NT_>
class Affine_form : public
     ::CGAL::Handle_with_policy<internal::Affine_form_rep<NT_> >
/*      LiS::Handle_without_union,
        ::std::allocator<internal::Affine_form_rep<NT_> > >*/
{
public:
    //! \name Typedefs 
    //@{ 
    //! coefficient type of this instance 
    typedef NT_ NT; 
    //! representation pointed to by this handle 
    typedef internal::Affine_form_rep<NT> Rep;
    //! base class  
    typedef ::CGAL::Handle_with_policy<Rep> Base;
    //! myself
    typedef Affine_form<NT> Self;
    //! container used to store noise symbols
    typedef typename Rep::Vector Vector;
    //! container used to store noise symbol indices
    typedef typename Rep::Int_vector Int_vector;
    //@}

protected:
    //! \name Protected methods
    //@{
    //! access to the internal noise symbol sequence
    Vector& coeffs() { return this->ptr()->coeff; }
    //! const access to the internal noise symbol sequence
    const Vector& coeffs() const { return this->ptr()->coeff; }
    //! access to the internal noise symbol sequence
    Int_vector& indexes() { return this->ptr()->index; }
    //! const access to the internal noise symbol sequence
    const Int_vector& indexes() const { return this->ptr()->index; }
    
    //! create an empty affine_form with n_coeffs noise terms
    Affine_form(internal::Creation_tag f, const NT& center_, unsigned n_coeffs,
        unsigned n_indexes)
        : Base(Rep(f, center_, n_coeffs, n_indexes))
    { }
    //! non-const access to coefficient of noise symbol \c i
    /*! The polynomial's representation must not be shared between
     *  different handles when this function is called.
     *  This can be ensured by calling \c copy_on_write().
     *
     *  If assertions are enabled, the index \c i is range-checked.
     */
    NT& coeff(unsigned int i) {
        //NiX_precond(!this->is_shared() && i<(this->ptr()->coeff.size()));
        return this->ptr()->coeff[i]; 
    }
    unsigned& index(unsigned int i) {
        return this->ptr()->index[i]; 
    }
    //@}
public:
    //! \name Constructors
    //@{
    
    //! constructs an affine form from a constant
    Affine_form(const NT& v0 = 0) : Base(Rep(v0))
    { }

    // Create an Affine_form from an array of NTs   
    explicit Affine_form(const NT& v0, const NT * t1, const unsigned * t2,
                unsigned t) : Base(Rep(v0, t1, t2, t))
    { }
        
    Affine_form(const Self& p) : Base(static_cast<const Base&>(p))
    { }

    explicit Affine_form(const Self& p, NT alpha, NT dzeta, 
        NT delta) : Base(Rep(p, alpha, dzeta, delta))
    { }
    
    explicit Affine_form(const NT& left, const NT& right) :
            Base(Rep(left, right))
    { }
    
//     explicit Affine_form(const NiX::Interval& i) : 
//         Base(Rep(i.lower(), i.upper()))
//     { }
        
    virtual ~Affine_form()
    {  }
    //@}
public:
    //! \name Public Methods
    //@{
    //! a random access iterator pointing to the first noise symbol
    typename Vector::const_iterator coeff_begin() const 
    { return this->ptr()->coeff.begin(); }
    //! a random access iterator pointing beyond the last noise symbol
    typename Vector::const_iterator coeff_end()   const 
    { return this->ptr()->coeff.end(); }
    
    //! a random access iterator pointing to the first index
    typename Int_vector::const_iterator index_begin() const 
    { return this->ptr()->index.begin(); }
    //! a random access iterator pointing beyond the last index
    typename Int_vector::const_iterator index_end()   const 
    { return this->ptr()->index.end(); }
        
    //Self & operator = (const Self & p);
    //! computes the sum of two affine forms
    Self operator + (const Self& p) const 
    {
        internal::Creation_tag TAG;
        unsigned len = length() + p.length();
        Self res(TAG, center() + p.center(), len, len); 
        //res.ptr()->center = center() + p.center();
        
        typename Int_vector::const_iterator it1 = index_begin(),
            it2 = p.index_begin();
        typename Vector::const_iterator val_it1 = coeff_begin(),
            val_it2 = p.coeff_begin();
                    
            
        bool inc_i1 = false, inc_i2 = false;
        unsigned idx, i1, i2, i = 0;
        NT val;
        
        
        while(1) {
            i1 = (it1 != index_end() ? *it1 : -1u);
            i2 = (it2 != p.index_end() ? *it2 : -1u);
            if(i1 == i2) {
                if(i1 == -1u) // both iterators point to the end - stop
                    break;
                val = (*val_it1)+(*val_it2);
                idx = i1;
                inc_i1 = true, inc_i2 = true;
                // assume that indexes come in an ascending order
            } else if(i1 < i2) { 
                val = *val_it1;
                idx = i1;
                inc_i1 = true;
            } else {
                val = *val_it2;
                idx = i2;
                inc_i2 = true;
            }
            //res.coeffs().push_back(val);
            //res.indexes().push_back(idx);
            res.coeff(i) = val;
            res.index(i) = idx;
            i++;
            if(inc_i1) 
            { it1++, val_it1++, inc_i1 = false; }
            if(inc_i2) 
            { it2++, val_it2++, inc_i2 = false; }
        }
        return res;
    }

    //! computes the difference of two affine forms
    Self operator - (const Self& p) const 
    {
        internal::Creation_tag TAG;
        unsigned len = length() + p.length();
        Self res(TAG, center() - p.center(), len, len); 
        
        //internal::Creation_tag TAG;
        //Self res; 
        //res.ptr()->center = center() - p.center();
        //this->copy_on_write();
        typename Int_vector::const_iterator it1 = index_begin(),
                it2 = p.index_begin();
        typename Vector::const_iterator val_it1 = coeff_begin(),
                val_it2 = p.coeff_begin();
        bool inc_i1 = false, inc_i2 = false;
        unsigned idx, i1, i2, i = 0;
        NT val;
        while(1) {
            i1 = (it1 != index_end() ? *it1 : -1u);
            i2 = (it2 != p.index_end() ? *it2 : -1u);
            if(i1 == i2) {
                if(i1 == -1u) // both iterators point to the end - stop
                    break;
                val = (*val_it1)-(*val_it2);
                idx = i1;
                inc_i1 = true, inc_i2 = true;
                // assume that indexes come in an ascending     order
            } else if(i1 < i2) { 
                val = *val_it1;
                idx = i1;
                inc_i1 = true;
            } else {
                val = -(*val_it2);
                idx = i2;
                inc_i2 = true;
            }
            res.coeff(i) = val;
            res.index(i) = idx;
            i++;
            //res.coeffs().push_back(val);
            //res.indexes().push_back(idx);
            if(inc_i1) 
            { it1++, val_it1++, inc_i1 = false; }
            if(inc_i2) 
            { it2++, val_it2++, inc_i2 = false; }
        }
        return res;
    }

    //! computes the product of two affine forms
    Self operator * (const Self& p) const 
    {
        NT cnt1 = center(), cnt2 = p.center();
        internal::Creation_tag TAG;
        unsigned len = length() + p.length()+1;
        Self res(TAG, cnt1 * cnt2, len, len); 
    
        /*Self res; 
        NT cnt1 = center(), cnt2 = p.center();
        res.ptr()->center = cnt1 * cnt2;*/
        //this->copy_on_write();
        typename Int_vector::const_iterator it1 = index_begin(),
                it2 = p.index_begin();
        typename Vector::const_iterator val_it1 = coeff_begin(),
                val_it2 = p.coeff_begin();
        bool inc_i1 = false, inc_i2 = false;
        unsigned idx, i1, i2, i=0;
        NT val;
        while(1) {
            i1 = (it1 != index_end() ? *it1 : -1u);
            i2 = (it2 != p.index_end() ? *it2 : -1u);
            if(i1 == i2) {
                if(i1 == -1u) // both iterators point to the end - stop
                    break;
                val = cnt1*(*val_it2) + cnt2*(*val_it1);
                idx = i1;
                inc_i1 = true, inc_i2 = true;
                // assume that indexes come in an ascending order
            } else if(i1 < i2) { 
                val = cnt2*(*val_it1);
                idx = i1;
                inc_i1 = true;
            } else {
                val = cnt1*(*val_it2);
                idx = i2;
                inc_i2 = true;
            }
            //res.coeffs().push_back(val);
            //res.indexes().push_back(idx);
            res.coeff(i) = val;
            res.index(i) = idx;
            i++;
            if(inc_i1) 
            { it1++, val_it1++, inc_i1 = false; }
            if(inc_i2) 
            { it2++, val_it2++, inc_i2 = false; }
        }
        // add a approximation of non-affine coeff
        res.coeff(i) = spread()*(p.spread());
        res.index(i) = this->ptr()->inc_last();
        
        //res.indexes().push_back(this->ptr()->inc_last());
        //res.coeffs().push_back(spread()*(p.spread()));
        return res;
    }

    //! unary operator
    Self operator -() const
    {
        Self res(*this);
        // from this point on representations must not be shared
        res.copy_on_write();
        res.ptr()->center = -res.center();
        
        typename Vector::iterator val_it = res.coeffs().begin();
        while(val_it != res.coeffs().end())
            *val_it = -(*val_it++); // does it work ??
        return res;
    }
    
    //! multiply by a constant factor
    Self operator *(const NT& x) const
    {
        Self res(*this);
        // from this point on representations must not be shared
        res.copy_on_write();
        res.ptr()->center = res.center() * x;
        
        typename Vector::iterator val_it = res.coeffs().begin();
        while(val_it != res.coeffs().end())
            *val_it = (*val_it++) * x; // does it work ??
        return res;
    }
    
    //! adds a constant to an affine form
    Self operator +(const NT& x) const
    {
        Self res(*this);
        // from this point on representations must not be shared
        res.copy_on_write();
        res.ptr()->center = res.center() + x;
        return res;
    }
    
    //! subtracts a constant from an affine form
    Self operator -(const NT& x) const
    {
        Self res(*this);
        // from this point on representations must not be shared
        res.copy_on_write();
        res.ptr()->center = res.center() - x;
        return res;
    }
   
    //! an output operator
    friend std::ostream & operator << <>(std::ostream & s, const Self &p);
    
    //! returns the number of noise symbols
    unsigned length() const { return this->ptr()->coeff.size(); }
    
    //! returns the center value of an affine form
    NT center() const { return this->ptr()->center; }
    
    //! converts an affine form to interval representation
    void convert(NT& bottom, NT& top) const 
    {
        NT rad = spread();
        bottom = center() - rad; 
        top = center() + rad;
    }
    
    //! returns the total deviation of an affine form, 
    //! i.e. the sum of all noise symbols (their abs values)
    NT spread() const
    {
        NT sum(0), val;
        typename Vector::const_iterator val_it = coeff_begin();
        while(val_it != coeff_end()) {
            val = *val_it++;
            sum += (val >= NT(0)? val: -val); 
        }
        return sum;
    }
    //@}

}; // class Affine_form<NT_>

//! output operator
template <class NT>
std::ostream& operator << (std::ostream& os, const Affine_form<NT>& p)
{
    os << "length: " << p.length() << "\n";
    os << "center: " << p.center() << "\n";
    for(unsigned i=0; i < p.length() ; i++)
        os << "e" << p.index(i) << " = " << p.coeff(i) << "\n";
    return os;
}

//! multiply by a constant factor (the left case)
template <class NT>
Affine_form<NT> operator * (const NT& x, const Affine_form<NT>& p) 
{
    return (p * x);
}

//! add a constant (the left case)
template <class NT>
Affine_form<NT> operator + (const NT& x, const Affine_form<NT>& p) 
{
    return (p + x);
}

//! subtract a constant (the left case)
template <class NT>
Affine_form<NT> operator - (const NT& x, const Affine_form<NT>& p) 
{
    return (Affine_form<NT>(x) - p);
}

//! compute a square of an affine form
template <class NT>
Affine_form<NT> sqr(const Affine_form<NT>& p) 
{
    return (p * p);
}

//! compute a power of an affine form (only for integer exponent)
template <class NT>
Affine_form<NT> pow(const Affine_form<NT>& p, int exp) 
{
    if(exp == 0) 
        return Affine_form<NT>(NT(1));
    else if (exp > 0) {
        if(exp & 1)
            return sqr(pow(p, exp>>1))*p;
        else
            return sqr(pow(p, exp>>1));
    } else {
        return Affine_form<NT>(NT(0));
    }
}

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_AFFINE_ARITHMETIC_H
// EOF
